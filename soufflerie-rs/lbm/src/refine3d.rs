//! Raffinement local de la soufflerie 3D : grilles D3Q19 imbriquées autour de
//! la maquette, chacune deux fois plus fine en espace et en temps que sa
//! parente (même schéma qu'en 2D, voir refine.rs). Couplage par les moments
//! (ρ, u, partie hors équilibre de Π remise à l'échelle, ν̃) : la grille fine
//! a une enveloppe de mailles fantômes remplies par interpolation cubique
//! (dans le plan de chaque face) et linéaire en temps ; sous la zone
//! raffinée, la grille parente reçoit l'état de la grille fine. Les efforts
//! viennent de la grille la plus fine, l'affichage de la veine.

use std::fmt::Write as _;

use bytemuck::{Pod, Zeroable};

use crate::gpu3d::{f, w3, Gpu3D, C3};

/// Maille fantôme d'une grille fine (type de maille).
pub const GHOST: u32 = 0xffff_fffe;
/// Marge (mailles parentes) entre le bord de la grille fine et la zone recopiée.
const MARGIN: usize = 2;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct LinkParams {
    pnx: u32,
    pny: u32,
    pnz: u32,
    pn: u32,
    cnx: u32,
    cny: u32,
    cnz: u32,
    cn: u32,
    ox: u32,
    oy: u32,
    oz: u32,
    ng: u32,
    rx: u32,
    ry: u32,
    rz: u32,
    margin: u32,
    sa: u32,
    alpha: f32,
    nu_p: f32,
    nu_c: f32,
    pad: [u32; 4],
}

const LINK_WGSL: &str = "
struct L {
  pnx: u32, pny: u32, pnz: u32, pn: u32,
  cnx: u32, cny: u32, cnz: u32, cn: u32,
  ox: u32, oy: u32, oz: u32, ng: u32,
  rx: u32, ry: u32, rz: u32, margin: u32,
  sa: u32, alpha: f32, nu_p: f32, nu_c: f32,
  pad0: u32, pad1: u32, pad2: u32, pad3: u32,
};
@group(0) @binding(0) var<uniform> P: L;
const CV1_3 = 357.911;
fn nu_turb(v: f32, nu: f32) -> f32 {
  if (v <= 0.0) { return 0.0; }
  let chi = v / nu;
  let c3 = chi * chi * chi;
  return v * c3 / (c3 + CV1_3);
}
// a = (ρ, u), b = (Π_xx, Π_yy, Π_zz hors équilibre, ν̃), c = (Π_xy, Π_xz, Π_yz, τ).
struct Mom {
  a: vec4<f32>,
  b: vec4<f32>,
  c: vec4<f32>,
};
fn flat(gid: vec3<u32>, nwg: vec3<u32>) -> u32 {
  return gid.x + gid.y * nwg.x * 64u;
}
";

/// Constantes D3Q19 et fonctions communes (moments, populations régularisées).
fn common() -> String {
    let list = |v: Vec<String>| v.join(", ");
    let cx = list((0..19).map(|i| C3[i][0].to_string()).collect());
    let cy = list((0..19).map(|i| C3[i][1].to_string()).collect());
    let cz = list((0..19).map(|i| C3[i][2].to_string()).collect());
    let w = list((0..19).map(|i| f(w3(i))).collect());
    let mut regul = String::new();
    for i in 0..19 {
        let [x, y, z] = C3[i].map(|v| v as f64);
        let cu: Vec<String> = [(x, "u.x"), (y, "u.y"), (z, "u.z")]
            .iter()
            .filter(|(c, _)| *c != 0.0)
            .map(|(c, s)| if *c > 0.0 { s.to_string() } else { format!("-{s}") })
            .collect();
        let cu = if cu.is_empty() { "0.0".into() } else { cu.join(" + ").replace("+ -", "- ") };
        let qp = format!(
            "{} * bxx + {} * byy + {} * bzz + {} * bxy + {} * bxz + {} * byz",
            f(x * x - 1.0 / 3.0),
            f(y * y - 1.0 / 3.0),
            f(z * z - 1.0 / 3.0),
            f(2.0 * x * y),
            f(2.0 * x * z),
            f(2.0 * y * z)
        );
        let h = |a: f64, b: f64| a * (b * b - 1.0 / 3.0);
        let h3: Vec<String> = [
            (h(y, x), "axxy"),
            (h(x, y), "axyy"),
            (h(z, x), "axxz"),
            (h(x, z), "axzz"),
            (h(z, y), "ayyz"),
            (h(y, z), "ayzz"),
        ]
        .iter()
        .filter(|(c, _)| c.abs() > 1e-12)
        .map(|(c, v)| format!("{} * {v}", f(*c)))
        .collect();
        let h3 = if h3.is_empty() { "0.0".into() } else { h3.join(" + ") };
        let wi = w3(i);
        let _ = writeln!(
            regul,
            "  {{ let cu = {cu}; p[{i}] = {} * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - usq) + {} * ({qp}) + {} * ({h3}); }}",
            f(wi),
            f(4.5 * wi),
            f(13.5 * wi)
        );
    }
    format!(
        "{LINK_WGSL}
const CX = array<i32, 19>({cx});
const CY = array<i32, 19>({cy});
const CZ = array<i32, 19>({cz});
const W = array<f32, 19>({w});

fn moments(f: array<f32, 19>) -> Mom {{
  var rho = 0.0;
  var j = vec3<f32>(0.0);
  var pxx = 0.0; var pyy = 0.0; var pzz = 0.0; var pxy = 0.0; var pxz = 0.0; var pyz = 0.0;
  for (var i = 0u; i < 19u; i++) {{
    let c = vec3<f32>(f32(CX[i]), f32(CY[i]), f32(CZ[i]));
    let v = f[i];
    rho += v;
    j += c * v;
    pxx += c.x * c.x * v; pyy += c.y * c.y * v; pzz += c.z * c.z * v;
    pxy += c.x * c.y * v; pxz += c.x * c.z * v; pyz += c.y * c.z * v;
  }}
  let u = j / rho;
  return Mom(
    vec4<f32>(rho, u),
    vec4<f32>(pxx - rho * (1.0 / 3.0 + u.x * u.x), pyy - rho * (1.0 / 3.0 + u.y * u.y), pzz - rho * (1.0 / 3.0 + u.z * u.z), 0.0),
    vec4<f32>(pxy - rho * u.x * u.y, pxz - rho * u.x * u.z, pyz - rho * u.y * u.z, 0.0));
}}

// Populations régularisées récursives (comme la collision) pour ρ, u et la
// partie hors équilibre (b) d'ordre 2.
fn regul(rho: f32, u: vec3<f32>, bxx: f32, byy: f32, bzz: f32, bxy: f32, bxz: f32, byz: f32) -> array<f32, 19> {{
  let usq = 1.5 * dot(u, u);
  let axxy = rho * u.x * u.x * u.y + u.y * bxx + 2.0 * u.x * bxy;
  let axyy = rho * u.x * u.y * u.y + u.x * byy + 2.0 * u.y * bxy;
  let axxz = rho * u.x * u.x * u.z + u.z * bxx + 2.0 * u.x * bxz;
  let axzz = rho * u.x * u.z * u.z + u.x * bzz + 2.0 * u.z * bxz;
  let ayyz = rho * u.y * u.y * u.z + u.z * byy + 2.0 * u.y * byz;
  let ayzz = rho * u.y * u.z * u.z + u.y * bzz + 2.0 * u.z * byz;
  var p: array<f32, 19>;
{regul}  return p;
}}
"
    )
}

/// Accès aux populations `name` (deux tampons A et B, plans de `n` mailles),
/// stockées g = f − w en simple ou demi-précision ; `ld_` rend f.
fn access(name: &str, ba: u32, bb: u32, fp16: bool, write: bool, n: &str) -> String {
    let mode = if write { "read_write" } else { "read" };
    let mut s = format!(
        "@group(0) @binding({ba}) var<storage, {mode}> {name}A: array<u32>;
@group(0) @binding({bb}) var<storage, {mode}> {name}B: array<u32>;
"
    );
    if fp16 {
        s += &format!(
            "fn ld_{name}(i: u32, c: u32) -> f32 {{
  let w = i >> 1u;
  var v: vec2<f32>;
  if (w < 5u) {{ v = unpack2x16float({name}A[w * {n} + c]); }} else {{ v = unpack2x16float({name}B[(w - 5u) * {n} + c]); }}
  return select(v.x, v.y, (i & 1u) == 1u) + W[i];
}}
"
        );
        if write {
            s += &format!("fn st_{name}(c: u32, p: array<f32, 19>) {{\n");
            for wd in 0..10 {
                let (a, b) = (2 * wd, 2 * wd + 1);
                let second = if b < 19 { format!("p[{b}] - W[{b}]") } else { "0.0".into() };
                let (buf, k) = if wd < 5 { ("A", wd) } else { ("B", wd - 5) };
                s += &format!("  {name}{buf}[{k}u * {n} + c] = pack2x16float(vec2<f32>(p[{a}] - W[{a}], {second}));\n");
            }
            s += "}\n";
        }
    } else {
        s += &format!(
            "fn ld_{name}(i: u32, c: u32) -> f32 {{
  if (i < 10u) {{ return bitcast<f32>({name}A[i * {n} + c]) + W[i]; }}
  return bitcast<f32>({name}B[(i - 10u) * {n} + c]) + W[i];
}}
"
        );
        if write {
            s += &format!("fn st_{name}(c: u32, p: array<f32, 19>) {{\n");
            for i in 0..19 {
                let (buf, k) = if i < 10 { ("A", i) } else { ("B", i - 10) };
                s += &format!("  {name}{buf}[{k}u * {n} + c] = bitcast<u32>(p[{i}] - W[{i}]);\n");
            }
            s += "}\n";
        }
    }
    s
}

const GATHER: &str = "
@group(0) @binding(1) var<storage, read_write> gm: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> gidx: array<u32>;
@group(0) @binding(5) var<storage, read> pnut: array<f32>;

fn pcell(x: i32, y: i32, z: i32) -> u32 {
  let q = clamp(vec3<i32>(x, y, z), vec3<i32>(0), vec3<i32>(i32(P.pnx), i32(P.pny), i32(P.pnz)) - vec3<i32>(1));
  return u32(q.x) + P.pnx * (u32(q.y) + P.pny * u32(q.z));
}

// Moments avant collision au nœud parent, ramenés à la grille fine.
fn pre_at(x: i32, y: i32, z: i32) -> Mom {
  var f: array<f32, 19>;
  for (var i = 0u; i < 19u; i++) {
    f[i] = ld_pf(i, pcell(x - CX[i], y - CY[i], z - CZ[i]));
  }
  var m = moments(f);
  var v = 3.0 * P.nu_p;
  var vt = 0.0;
  if (P.sa != 0u) {
    v = pnut[pcell(x, y, z)];
    vt = nu_turb(v, P.nu_p);
  }
  let tau_p = 0.5 + 3.0 * (P.nu_p + vt);
  let tau_c = 0.5 + 3.0 * (P.nu_c + 2.0 * vt);
  let k = tau_c / (2.0 * tau_p);
  m.b = vec4<f32>(m.b.xyz * k, 2.0 * v);
  m.c = vec4<f32>(m.c.xyz * k, tau_c);
  return m;
}

struct St { n: i32, o: vec4<i32>, w: vec4<f32>, };
// Pochoir le long d'un axe : nœud parent (coordonnée paire), ou interpolation
// cubique sur quatre nœuds (impaire).
fn stencil(x2: i32) -> St {
  if ((x2 & 1) == 0) {
    return St(1, vec4<i32>(x2 / 2, 0, 0, 0), vec4<f32>(1.0, 0.0, 0.0, 0.0));
  }
  let i0 = (x2 - 1) / 2;
  return St(4, vec4<i32>(i0 - 1, i0, i0 + 1, i0 + 2), vec4<f32>(-0.0625, 0.5625, 0.5625, -0.0625));
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let g = flat(gid, nwg);
  if (g >= P.ng) { return; }
  gm[6u * g] = gm[6u * g + 3u];
  gm[6u * g + 1u] = gm[6u * g + 4u];
  gm[6u * g + 2u] = gm[6u * g + 5u];
  let c = gidx[g];
  let x = c % P.cnx;
  let y = (c / P.cnx) % P.cny;
  let z = c / (P.cnx * P.cny);
  let sx = stencil(i32(2u * P.ox + x));
  let sy = stencil(i32(2u * P.oy + y));
  let sz = stencil(i32(2u * P.oz + z));
  var m = Mom(vec4<f32>(0.0), vec4<f32>(0.0), vec4<f32>(0.0));
  for (var k = 0; k < sz.n; k++) {
    for (var j = 0; j < sy.n; j++) {
      for (var i = 0; i < sx.n; i++) {
        let w = sx.w[i] * sy.w[j] * sz.w[k];
        let q = pre_at(sx.o[i], sy.o[j], sz.o[k]);
        m.a += w * q.a;
        m.b += w * q.b;
        m.c += w * q.c;
      }
    }
  }
  gm[6u * g + 3u] = m.a;
  gm[6u * g + 4u] = m.b;
  gm[6u * g + 5u] = m.c;
}
";

const FILL: &str = "
@group(0) @binding(1) var<storage, read> gm: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> gidx: array<u32>;
@group(0) @binding(5) var<storage, read_write> cnut: array<f32>;
@group(0) @binding(6) var<storage, read_write> cmacro: array<vec4<f32>>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let g = flat(gid, nwg);
  if (g >= P.ng) { return; }
  let a = mix(gm[6u * g], gm[6u * g + 3u], P.alpha);
  let b = mix(gm[6u * g + 1u], gm[6u * g + 4u], P.alpha);
  let c = mix(gm[6u * g + 2u], gm[6u * g + 5u], P.alpha);
  let k = 1.0 - 1.0 / c.w;
  let p = regul(a.x, a.yzw, k * b.x, k * b.y, k * b.z, k * c.x, k * c.y, k * c.z);
  let idx = gidx[g];
  st_cf(idx, p);
  if (P.sa != 0u) { cnut[idx] = b.w; }
  cmacro[idx] = vec4<f32>(a.yzw, a.x);
}
";

const RESTRICT: &str = "
@group(0) @binding(5) var<storage, read> cwall: array<u32>;
@group(0) @binding(6) var<storage, read> cnut: array<f32>;
@group(0) @binding(9) var<storage, read_write> pnut: array<f32>;
@group(0) @binding(10) var<storage, read_write> pmacro: array<vec4<f32>>;
@group(0) @binding(11) var<storage, read> pwall: array<u32>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let r = flat(gid, nwg);
  if (r >= P.rx * P.ry * P.rz) { return; }
  let px = P.ox + P.margin + r % P.rx;
  let py = P.oy + P.margin + (r / P.rx) % P.ry;
  let pz = P.oz + P.margin + r / (P.rx * P.ry);
  let pidx = px + P.pnx * (py + P.pny * pz);
  if (pwall[pidx] == 0xffffffffu) { return; }
  let cx = 2u * (px - P.ox);
  let cy = 2u * (py - P.oy);
  let cz = 2u * (pz - P.oz);
  let cidx = cx + P.cnx * (cy + P.cny * cz);
  var v = 0.0;
  var vt = 0.0;
  if (P.sa != 0u) {
    v = cnut[cidx];
    vt = nu_turb(v, P.nu_c);
  }
  var f: array<f32, 19>;
  var p: array<f32, 19>;
  var m: Mom;
  if (cwall[cidx] == 0u) {
    // Populations avant collision au temps t + 1, tirées du demi-pas t + ½.
    for (var i = 0u; i < 19u; i++) {
      let s = vec3<i32>(i32(cx) - CX[i], i32(cy) - CY[i], i32(cz) - CZ[i]);
      f[i] = ld_cs(i, u32(s.x) + P.cnx * (u32(s.y) + P.cny * u32(s.z)));
    }
    m = moments(f);
    let tau_c = 0.5 + 3.0 * (P.nu_c + vt);
    let tau_p = 0.5 + 3.0 * (P.nu_p + 0.5 * vt);
    let k = 2.0 * (tau_p - 1.0) / tau_c;
    p = regul(m.a.x, m.a.yzw, k * m.b.x, k * m.b.y, k * m.b.z, k * m.c.x, k * m.c.y, k * m.c.z);
  } else {
    // Maille pariétale ou solide dans la grille fine : équilibre (affichage).
    for (var i = 0u; i < 19u; i++) {
      f[i] = ld_cc(i, cidx);
    }
    m = moments(f);
    p = regul(m.a.x, m.a.yzw, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
  st_pf(pidx, p);
  if (P.sa != 0u) { pnut[pidx] = 0.5 * v; }
  pmacro[pidx] = vec4<f32>(m.a.yzw, m.a.x);
}
";

/// Couplage d'une grille (parente) et de la suivante (fine).
pub(crate) struct Link3 {
    ox: usize,
    oy: usize,
    oz: usize,
    ng: u32,
    r: [usize; 3],
    params: [wgpu::Buffer; 2],
    gm: wgpu::Buffer,
    gidx: wgpu::Buffer,
    gather: wgpu::ComputePipeline,
    fill: wgpu::ComputePipeline,
    restrict: wgpu::ComputePipeline,
    gather_bind: Vec<wgpu::BindGroup>,
    fill_bind: Vec<wgpu::BindGroup>,
    restrict_bind: Vec<wgpu::BindGroup>,
}

fn disp(n: u32) -> (u32, u32) {
    let groups = n.div_ceil(64).max(1);
    let dx = groups.min(65535);
    (dx, groups.div_ceil(dx))
}

impl Link3 {
    pub(crate) fn new(parent: &Gpu3D, child: &Gpu3D, o: [usize; 3]) -> Link3 {
        let d = &parent.device;
        let fp16 = parent.fp16;
        let (cnx, cny, cnz) = (child.nx, child.ny, child.nz);
        let mut ghosts = Vec::new();
        for z in 0..cnz {
            for y in 0..cny {
                for x in 0..cnx {
                    if x == 0 || y == 0 || z == 0 || x == cnx - 1 || y == cny - 1 || z == cnz - 1 {
                        ghosts.push((x + cnx * (y + cny * z)) as u32);
                    }
                }
            }
        }
        let ng = ghosts.len() as u32;
        let base = common();
        let module = |label: &str, body: String| {
            d.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some(label),
                source: wgpu::ShaderSource::Wgsl(format!("{base}{body}").into()),
            })
        };
        let gather = module("gather3d", format!("{}{GATHER}", access("pf", 3, 4, fp16, false, "P.pn")));
        let fill = module("fill3d", format!("{}{FILL}", access("cf", 3, 4, fp16, true, "P.cn")));
        let restrict = module(
            "restrict3d",
            format!(
                "{}{}{}{RESTRICT}",
                access("cs", 1, 2, fp16, false, "P.cn"),
                access("cc", 3, 4, fp16, false, "P.cn"),
                access("pf", 7, 8, fp16, true, "P.pn")
            ),
        );
        let pipe = |m: &wgpu::ShaderModule, label: &str| {
            d.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(label),
                layout: None,
                module: m,
                entry_point: Some("main"),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        use wgpu::util::DeviceExt;
        let uni = wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST;
        let mk = |label: &str, size: u64, usage| {
            d.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: size.max(16),
                usage,
                mapped_at_creation: false,
            })
        };
        let half = [(cnx - 1) / 2, (cny - 1) / 2, (cnz - 1) / 2];
        Link3 {
            ox: o[0],
            oy: o[1],
            oz: o[2],
            ng,
            r: half.map(|h| h + 1 - 2 * MARGIN),
            params: [
                mk("link3_0", std::mem::size_of::<LinkParams>() as u64, uni),
                mk("link3_1", std::mem::size_of::<LinkParams>() as u64, uni),
            ],
            gm: mk("ghost3d", 96 * ng as u64, wgpu::BufferUsages::STORAGE),
            gidx: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("ghost_idx3d"),
                contents: bytemuck::cast_slice(&ghosts),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            gather: pipe(&gather, "gather3d"),
            fill: pipe(&fill, "fill3d"),
            restrict: pipe(&restrict, "restrict3d"),
            gather_bind: Vec::new(),
            fill_bind: Vec::new(),
            restrict_bind: Vec::new(),
        }
    }

    /// Paramètres et groupes de liaison (à refaire si un tampon change).
    pub(crate) fn prepare(&mut self, p: &Gpu3D, c: &Gpu3D) {
        for (a, buf) in self.params.iter().enumerate() {
            let lp = LinkParams {
                pnx: p.nx as u32,
                pny: p.ny as u32,
                pnz: p.nz as u32,
                pn: p.n as u32,
                cnx: c.nx as u32,
                cny: c.ny as u32,
                cnz: c.nz as u32,
                cn: c.n as u32,
                ox: self.ox as u32,
                oy: self.oy as u32,
                oz: self.oz as u32,
                ng: self.ng,
                rx: self.r[0] as u32,
                ry: self.r[1] as u32,
                rz: self.r[2] as u32,
                margin: MARGIN as u32,
                sa: (p.sa && c.sa) as u32,
                alpha: 0.5 * a as f32,
                nu_p: p.nu as f32,
                nu_c: c.nu as f32,
                pad: [0; 4],
            };
            p.queue.write_buffer(buf, 0, bytemuck::bytes_of(&lp));
        }
        let d = &p.device;
        let group = |pipe: &wgpu::ComputePipeline, res: Vec<(u32, &wgpu::Buffer)>| {
            let entries: Vec<wgpu::BindGroupEntry> = res
                .into_iter()
                .map(|(binding, b)| wgpu::BindGroupEntry {
                    binding,
                    resource: b.as_entire_binding(),
                })
                .collect();
            d.create_bind_group(&wgpu::BindGroupDescriptor {
                label: None,
                layout: &pipe.get_bind_group_layout(0),
                entries: &entries,
            })
        };
        // Après un pas, `cur` désigne le nouvel état ; l'ancien fournit les tirages.
        self.gather_bind = (0..2)
            .map(|cur| {
                group(
                    &self.gather,
                    vec![
                        (0, &self.params[0]),
                        (1, &self.gm),
                        (2, &self.gidx),
                        (3, &p.pops[1 - cur][0]),
                        (4, &p.pops[1 - cur][1]),
                        (5, &p.nut[cur]),
                    ],
                )
            })
            .collect();
        self.fill_bind = (0..4)
            .map(|i| {
                let (cur, a) = (i / 2, i % 2);
                group(
                    &self.fill,
                    vec![
                        (0, &self.params[a]),
                        (1, &self.gm),
                        (2, &self.gidx),
                        (3, &c.pops[cur][0]),
                        (4, &c.pops[cur][1]),
                        (5, &c.nut[cur]),
                        (6, &c.macro_buf),
                    ],
                )
            })
            .collect();
        self.restrict_bind = (0..4)
            .map(|i| {
                let (pc, cc) = (i / 2, i % 2);
                group(
                    &self.restrict,
                    vec![
                        (0, &self.params[0]),
                        (1, &c.pops[1 - cc][0]),
                        (2, &c.pops[1 - cc][1]),
                        (3, &c.pops[cc][0]),
                        (4, &c.pops[cc][1]),
                        (5, &c.wall),
                        (6, &c.nut[cc]),
                        (7, &p.pops[pc][0]),
                        (8, &p.pops[pc][1]),
                        (9, &p.nut[pc]),
                        (10, &p.macro_buf),
                        (11, &p.wall),
                    ],
                )
            })
            .collect();
    }

    fn gather(&self, pass: &mut wgpu::ComputePass, pcur: usize) {
        let (dx, dy) = disp(self.ng);
        pass.set_pipeline(&self.gather);
        pass.set_bind_group(0, &self.gather_bind[pcur], &[]);
        pass.dispatch_workgroups(dx, dy, 1);
    }
}

/// Grilles fines d'une veine 3D et leurs couplages.
pub(crate) struct Refinement {
    /// Grilles 1, 2… (2ᵏ fois plus fines que la veine).
    pub(crate) levels: Vec<Gpu3D>,
    pub(crate) links: Vec<Link3>,
    pub(crate) fresh: bool,
    /// Nœud 0 de chaque grille fine, en mailles de la veine.
    pub(crate) origin: Vec<[f64; 3]>,
}

/// Emprise des grilles fines autour de la maquette, en longueurs de
/// référence de la veine : amont, aval, de part et d'autre.
pub(crate) const EXTENT3: [(f64, f64, f64); 3] = [(0.6, 1.5, 0.6), (0.35, 0.7, 0.4), (0.2, 0.4, 0.25)];
pub(crate) const MARGIN3: usize = MARGIN;

/// Un pas de la grille k, suivi de deux pas de chaque grille plus fine.
fn encode(levels: &mut [&mut Gpu3D], links: &[Link3], pass: &mut wgpu::ComputePass, k: usize, slot: u32) {
    {
        let lev = &mut *levels[k];
        pass.set_pipeline(&lev.step_pipe);
        pass.set_bind_group(0, &lev.bind[lev.cur], &[slot]);
        pass.dispatch_workgroups(lev.disp.0, lev.disp.1, 1);
        if lev.sa {
            pass.set_pipeline(&lev.sa_pipe);
            pass.set_bind_group(0, &lev.sa_bind[lev.cur], &[slot]);
            pass.dispatch_workgroups(lev.disp.0, lev.disp.1, 1);
        }
        lev.cur = 1 - lev.cur;
    }
    let Some(link) = links.get(k) else { return };
    link.gather(pass, levels[k].cur);
    let (dx, dy) = disp(link.ng);
    for a in 0..2 {
        pass.set_pipeline(&link.fill);
        pass.set_bind_group(0, &link.fill_bind[levels[k + 1].cur * 2 + a], &[]);
        pass.dispatch_workgroups(dx, dy, 1);
        encode(levels, links, pass, k + 1, slot);
    }
    let nr = (link.r[0] * link.r[1] * link.r[2]) as u32;
    let (rx, ry) = disp(nr);
    pass.set_pipeline(&link.restrict);
    pass.set_bind_group(0, &link.restrict_bind[levels[k].cur * 2 + levels[k + 1].cur], &[]);
    pass.dispatch_workgroups(rx, ry, 1);
}

/// `n` pas de la veine `top` et de ses grilles fines.
pub(crate) fn run(top: &mut Gpu3D, r: &mut Refinement, n: usize, ramp: &crate::cpu::Ramp) {
    for k in 0..r.links.len() {
        let (p, c) = if k == 0 {
            (&*top, &r.levels[0])
        } else {
            (&r.levels[k - 1], &r.levels[k])
        };
        let link = &mut r.links[k];
        link.prepare(p, c);
    }
    let device = top.device.clone();
    let queue = top.queue.clone();
    let mut all: Vec<&mut Gpu3D> = std::iter::once(top).chain(r.levels.iter_mut()).collect();
    if r.fresh {
        let mut enc = device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            for (k, link) in r.links.iter().enumerate() {
                for _ in 0..2 {
                    link.gather(&mut pass, all[k].cur);
                }
            }
        }
        queue.submit([enc.finish()]);
        r.fresh = false;
    }
    let cost: usize = all.iter().enumerate().map(|(k, g)| g.n << k).sum();
    let per_submit = (30_000_000 / cost.max(1)).clamp(1, crate::gpu3d::BATCH);
    let mut left = n;
    while left > 0 {
        let kc = left.min(per_submit);
        let t0 = all[0].time;
        for (k, g) in all.iter().enumerate() {
            let mut bytes = vec![0u8; crate::gpu3d::SLOT as usize * kc];
            for s in 0..kc {
                let (u, a) = ramp.at(t0 + s as u64);
                let last = k == 0 && left == kc && s == kc - 1;
                let p = g.params_for(u, a / (1u64 << k) as f32, last);
                let o = s * crate::gpu3d::SLOT as usize;
                bytes[o..o + std::mem::size_of_val(&p)].copy_from_slice(bytemuck::bytes_of(&p));
            }
            queue.write_buffer(&g.params, 0, &bytes);
        }
        let mut enc = device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            for s in 0..kc {
                encode(&mut all, &r.links, &mut pass, 0, (s as u64 * crate::gpu3d::SLOT) as u32);
            }
        }
        queue.submit([enc.finish()]);
        for (k, g) in all.iter_mut().enumerate() {
            let steps = (kc as u64) << k;
            g.time += steps;
            g.force_steps += steps;
        }
        left -= kc;
    }
}
