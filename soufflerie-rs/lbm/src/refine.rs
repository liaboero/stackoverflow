//! Raffinement local sur GPU : grilles imbriquées autour de la maquette,
//! chacune deux fois plus fine (en espace et en temps) que sa parente. La
//! veine entière reste à la résolution L ; près du corps, la couche limite,
//! les bords d'attaque et de fuite sont résolus à 2ᴷ·L pour une fraction du
//! coût d'une grille uniforme (couplage : refine.wgsl). Les efforts viennent
//! de la grille la plus fine, l'affichage de la grille de la veine.

use bytemuck::{Pod, Zeroable};

use crate::cpu::{ForceSum, LesMode, Macros, Ramp, SolverConfig};
use crate::geometry::Geometry;
use crate::gpu::{GpuSolver, BATCH, SLOT};
use crate::shapes::Body;

const COMMON: &str = include_str!("refine.wgsl");

const GATHER: &str = "
@group(0) @binding(1) var<storage, read_write> gm: array<vec4<f32>>;
@group(0) @binding(3) var<storage, read> psa: array<vec4<f32>>;

// Moments avant collision au nœud parent (px, py), ramenés à la grille fine.
fn pre_at(px: i32, py: i32) -> Mom {
  var f: array<f32, 9>;
  for (var i = 0u; i < 9u; i++) {
    let sx = clamp(px - CX[i], 0, i32(P.pnx) - 1);
    let sy = clamp(py - CY[i], 0, i32(P.pny) - 1);
    f[i] = ld_pf(i, u32(sy) * P.pnx + u32(sx));
  }
  var m = moments(f);
  var v = 3.0 * P.nu_p;
  var vt = 0.0;
  if (P.sa_on != 0u) {
    v = psa[u32(py) * P.pnx + u32(px)].w;
    vt = nu_turb(v, P.nu_p);
  }
  let tau_p = 0.5 + 3.0 * (P.nu_p + vt);
  let tau_c = 0.5 + 3.0 * (P.nu_c + 2.0 * vt);
  let k = tau_c / (2.0 * tau_p);
  m.a.w = 2.0 * v;
  m.b = vec4<f32>(m.b.xyz * k, tau_c);
  return m;
}

fn mix4(m0: Mom, m1: Mom, m2: Mom, m3: Mom) -> Mom {
  let w0 = -0.0625;
  let w1 = 0.5625;
  return Mom(w0 * (m0.a + m3.a) + w1 * (m1.a + m2.a), w0 * (m0.b + m3.b) + w1 * (m1.b + m2.b));
}

@compute @workgroup_size(64, 1, 1)
fn gather(@builtin(global_invocation_id) gid: vec3<u32>) {
  let g = gid.x;
  if (g >= P.ng) {
    return;
  }
  gm[4u * g] = gm[4u * g + 2u];
  gm[4u * g + 1u] = gm[4u * g + 3u];
  let xy = ghost_xy(g);
  let x2 = i32(2u * P.ox + xy.x);
  let y2 = i32(2u * P.oy + xy.y);
  var m: Mom;
  if ((x2 & 1) == 1) {
    let i0 = (x2 - 1) / 2;
    let y = y2 / 2;
    m = mix4(pre_at(i0 - 1, y), pre_at(i0, y), pre_at(i0 + 1, y), pre_at(i0 + 2, y));
  } else if ((y2 & 1) == 1) {
    let j0 = (y2 - 1) / 2;
    let x = x2 / 2;
    m = mix4(pre_at(x, j0 - 1), pre_at(x, j0), pre_at(x, j0 + 1), pre_at(x, j0 + 2));
  } else {
    m = pre_at(x2 / 2, y2 / 2);
  }
  gm[4u * g + 2u] = m.a;
  gm[4u * g + 3u] = m.b;
}
";

const FILL: &str = "
@group(0) @binding(1) var<storage, read> gm: array<vec4<f32>>;
@group(0) @binding(3) var<storage, read_write> csa: array<vec4<f32>>;

@compute @workgroup_size(64, 1, 1)
fn fill(@builtin(global_invocation_id) gid: vec3<u32>) {
  let g = gid.x;
  if (g >= P.ng) {
    return;
  }
  let a = mix(gm[4u * g], gm[4u * g + 2u], P.alpha);
  let b = mix(gm[4u * g + 1u], gm[4u * g + 3u], P.alpha);
  let k = 1.0 - 1.0 / b.w;
  let p = regularized(a.x, a.y, a.z, k * b.x, k * b.y, k * b.z);
  let xy = ghost_xy(g);
  let c = xy.y * P.cnx + xy.x;
  st_cf(c, p);
  if (P.sa_on != 0u) {
    csa[c] = vec4<f32>(a.y, a.z, 0.0, a.w);
  }
}
";

const RESTRICT: &str = "
@group(0) @binding(4) var<storage, read> ckind: array<u32>;
@group(0) @binding(5) var<storage, read> csa: array<vec4<f32>>;
@group(0) @binding(7) var<storage, read_write> psa: array<vec4<f32>>;
@group(0) @binding(8) var<storage, read> pkind: array<u32>;

// Moments avant collision au nœud fin (x, y).
fn pre_c(x: i32, y: i32) -> Mom {
  var f: array<f32, 9>;
  for (var i = 0u; i < 9u; i++) {
    f[i] = ld_cs(i, u32(y - CY[i]) * P.cnx + u32(x - CX[i]));
  }
  return moments(f);
}

@compute @workgroup_size(64, 1, 1)
fn restrict_down(@builtin(global_invocation_id) gid: vec3<u32>) {
  if (gid.x >= P.nrx || gid.y >= P.nry) {
    return;
  }
  let px = P.ox + P.margin + gid.x;
  let py = P.oy + P.margin + gid.y;
  let pidx = py * P.pnx + px;
  if ((pkind[pidx] & 3u) == 1u) {
    return;
  }
  let cx = 2u * (px - P.ox);
  let cy = 2u * (py - P.oy);
  let cidx = cy * P.cnx + cx;
  var s = vec4<f32>(0.0);
  var vt = 0.0;
  if (P.sa_on != 0u) {
    s = csa[cidx];
    vt = nu_turb(s.w, P.nu_c);
  }
  var p: array<f32, 9>;
  var m: Mom;
  if ((ckind[cidx] & 3u) == 0u) {
    // Moments avant collision au temps t + 1 (tirés du demi-pas t + ½),
    // filtrés sur les 3 × 3 nœuds fins voisins quand ils sont tous fluides
    // (Lagrava et al. 2012) : la grille parente ne reçoit pas les modes
    // qu'elle ne peut pas représenter.
    var all = true;
    for (var j = -1; j <= 1; j++) {
      for (var i = -1; i <= 1; i++) {
        all = all && (ckind[u32(i32(cidx) + j * i32(P.cnx) + i)] & 3u) == 0u;
      }
    }
    if (all) {
      m = Mom(vec4<f32>(0.0), vec4<f32>(0.0));
      for (var j = -1; j <= 1; j++) {
        for (var i = -1; i <= 1; i++) {
          let w = f32((2 - abs(i)) * (2 - abs(j))) / 16.0;
          let q = pre_c(i32(cx) + i, i32(cy) + j);
          m.a += w * q.a;
          m.b += w * q.b;
        }
      }
    } else {
      m = pre_c(i32(cx), i32(cy));
    }
    let tau_c = 0.5 + 3.0 * (P.nu_c + vt);
    let tau_p = 0.5 + 3.0 * (P.nu_p + 0.5 * vt);
    let k = 2.0 * (tau_p - 1.0) / tau_c;
    p = regularized(m.a.x, m.a.y, m.a.z, k * m.b.x, k * m.b.y, k * m.b.z);
  } else {
    // Voisin de paroi ou solide dans la grille fine : équilibre (affichage).
    var f: array<f32, 9>;
    for (var i = 0u; i < 9u; i++) {
      f[i] = ld_cc(i, cidx);
    }
    m = moments(f);
    p = regularized(m.a.x, m.a.y, m.a.z, 0.0, 0.0, 0.0);
  }
  st_pf(pidx, p);
  if (P.sa_on != 0u) {
    psa[pidx] = vec4<f32>(m.a.y, m.a.z, s.z, 0.5 * s.w);
  }
}
";

/// Accès aux populations `name` (liaison `binding`, plans de `n` mailles),
/// en simple ou demi-précision (g = f − w, deux par mot).
fn pop_access(name: &str, binding: u32, half: bool, write: bool, n: &str) -> String {
    let mode = if write { "read_write" } else { "read" };
    let ty = if half { "u32" } else { "f32" };
    let mut s = format!("@group(0) @binding({binding}) var<storage, {mode}> {name}: array<{ty}>;\n");
    if half {
        s += &format!(
            "fn ld_{name}(i: u32, c: u32) -> f32 {{ let h = unpack2x16float({name}[(i >> 1u) * {n} + c]); return select(h.x, h.y, (i & 1u) == 1u) + W[i]; }}\n"
        );
        if write {
            s += &format!("fn st_{name}(c: u32, p: array<f32, 9>) {{\n");
            for k in 0..4 {
                s += &format!(
                    "  {name}[{k}u * {n} + c] = pack2x16float(vec2<f32>(p[{a}] - W[{a}], p[{b}] - W[{b}]));\n",
                    a = 2 * k,
                    b = 2 * k + 1
                );
            }
            s += &format!("  {name}[4u * {n} + c] = pack2x16float(vec2<f32>(p[8] - W[8], 0.0));\n}}\n");
        }
    } else {
        s += &format!("fn ld_{name}(i: u32, c: u32) -> f32 {{ return {name}[i * {n} + c]; }}\n");
        if write {
            s += &format!("fn st_{name}(c: u32, p: array<f32, 9>) {{\n");
            for i in 0..9 {
                s += &format!("  {name}[{i}u * {n} + c] = p[{i}];\n");
            }
            s += "}\n";
        }
    }
    s
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct LinkParams {
    pnx: u32,
    pny: u32,
    pn: u32,
    cnx: u32,
    cny: u32,
    cn: u32,
    ox: u32,
    oy: u32,
    ng: u32,
    margin: u32,
    sa_on: u32,
    nrx: u32,
    alpha: f32,
    nu_p: f32,
    nu_c: f32,
    nry: u32,
}

/// Couplage d'une grille (parente) et de la suivante (fine).
struct Link {
    /// Coin de la grille fine dans la grille parente (nœuds parents).
    ox: usize,
    oy: usize,
    ng: u32,
    margin: usize,
    nrx: u32,
    nry: u32,
    /// Paramètres pour alpha = 0 et alpha = ½.
    params: [wgpu::Buffer; 2],
    /// Par maille fantôme : moments au temps t puis t + 1 (2 × 2 vec4).
    gm: wgpu::Buffer,
    gather: wgpu::ComputePipeline,
    fill: wgpu::ComputePipeline,
    restrict: wgpu::ComputePipeline,
    /// Par indice courant de la grille parente.
    gather_bind: Vec<wgpu::BindGroup>,
    /// Indice courant de la grille fine × 2 + alpha.
    fill_bind: Vec<wgpu::BindGroup>,
    /// Indice courant de la parente × 2 + indice courant de la fine.
    restrict_bind: Vec<wgpu::BindGroup>,
}

/// Emprise des grilles fines (en longueurs de référence autour du point de
/// référence de la maquette) : amont, aval, demi-hauteur. Deux niveaux au
/// plus : à 8 L, la couche limite laminaire (τ ≈ 0,5, deux mailles
/// d'épaisseur) devient instable.
const EXTENT: [(f64, f64, f64); 2] = [(0.9, 1.9, 0.75), (0.65, 1.5, 0.58)];
/// Marge (mailles parentes) entre le bord de la grille fine et la zone recopiée.
const MARGIN: usize = 2;

pub struct RefinedGpu {
    /// 0 : veine entière ; k : grille 2ᵏ fois plus fine.
    pub levels: Vec<GpuSolver>,
    /// Maquette placée dans chaque grille.
    pub bodies: Vec<Body>,
    /// Coin de chaque grille dans la veine (mailles de la veine).
    pub origin: Vec<[f64; 2]>,
    links: Vec<Link>,
    u_in: f32,
    fresh: bool,
}

impl RefinedGpu {
    /// `depth` grilles fines autour du point (x0, y0) de la veine, L mailles
    /// par longueur de référence dans la veine.
    pub fn new(cfg: SolverConfig, depth: usize, l: usize, x0: f64, y0: f64) -> Result<RefinedGpu, String> {
        let depth = depth.clamp(1, EXTENT.len());
        let base = GpuSolver::new(cfg.clone())?;
        let device = base.device.clone();
        let queue = base.queue.clone();
        let name = base.adapter_name().to_string();
        let half = cfg.half;
        let mut levels = vec![base];
        let mut origin = vec![[0.0, 0.0]];
        let mut links = Vec::new();
        let lf = l as f64;
        for k in 1..=depth {
            let (up, down, h) = EXTENT[k - 1];
            let sp = (1usize << (k - 1)) as f64;
            let [gx, gy] = origin[k - 1];
            let parent = &levels[k - 1];
            let (pnx, pny) = (parent.cfg.nx, parent.cfg.ny);
            let to_p = |v: f64, g: f64| (v - g) * sp;
            let px0 = to_p(x0 - up * lf, gx).floor().max(2.0) as usize;
            let px1 = (to_p(x0 + down * lf, gx).ceil() as usize).min(pnx - 3);
            let py0 = to_p(y0 - h * lf, gy).floor().max(2.0) as usize;
            let py1 = (to_p(y0 + h * lf, gy).ceil() as usize).min(pny - 3);
            if px1 <= px0 + 2 * MARGIN + 2 || py1 <= py0 + 2 * MARGIN + 2 {
                return Err("grille fine trop petite".into());
            }
            let (cnx, cny) = (2 * (px1 - px0) + 1, 2 * (py1 - py0) + 1);
            let s = (1usize << k) as f64;
            let ccfg = SolverConfig {
                nx: cnx,
                ny: cny,
                nu: cfg.nu * s,
                les_wall_a: cfg.les_wall_a * s,
                sponge_fraction: 0.0,
                ..cfg.clone()
            };
            let child = GpuSolver::with_device(ccfg, device.clone(), queue.clone(), name.clone(), true);
            origin.push([gx + px0 as f64 / sp, gy + py0 as f64 / sp]);
            links.push(Self::link(&levels[k - 1], &child, px0, py0, half));
            levels.push(child);
        }
        let mut r = RefinedGpu {
            bodies: Vec::new(),
            levels,
            origin,
            links,
            u_in: cfg.u_in,
            fresh: true,
        };
        r.sync_les();
        Ok(r)
    }

    fn link(parent: &GpuSolver, child: &GpuSolver, ox: usize, oy: usize, half: bool) -> Link {
        let d = &parent.device;
        let (cnx, cny) = (child.cfg.nx, child.cfg.ny);
        let ng = 2 * cnx + 2 * (cny - 2);
        let module = |label: &str, body: String| {
            d.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some(label),
                source: wgpu::ShaderSource::Wgsl(format!("{COMMON}{body}").into()),
            })
        };
        let gather_mod = module("gather", format!("{}{GATHER}", pop_access("pf", 2, half, false, "P.pn")));
        let fill_mod = module("fill", format!("{}{FILL}", pop_access("cf", 2, half, true, "P.cn")));
        let restrict_mod = module(
            "restrict",
            format!(
                "{}{}{}{RESTRICT}",
                pop_access("cs", 2, half, false, "P.cn"),
                pop_access("cc", 3, half, false, "P.cn"),
                pop_access("pf", 6, half, true, "P.pn")
            ),
        );
        let pipe = |m: &wgpu::ShaderModule, entry: &str| {
            d.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(entry),
                layout: None,
                module: m,
                entry_point: Some(entry),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        let mk = |label: &str, size: u64, usage: wgpu::BufferUsages| {
            d.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size,
                usage,
                mapped_at_creation: false,
            })
        };
        let uni = wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST;
        let half_cells = (cnx - 1) / 2;
        let rows = (cny - 1) / 2;
        Link {
            ox,
            oy,
            ng: ng as u32,
            margin: MARGIN,
            nrx: (half_cells + 1 - 2 * MARGIN) as u32,
            nry: (rows + 1 - 2 * MARGIN) as u32,
            params: [
                mk("link0", std::mem::size_of::<LinkParams>() as u64, uni),
                mk("link1", std::mem::size_of::<LinkParams>() as u64, uni),
            ],
            gm: mk(
                "ghost_moments",
                (64 * ng) as u64,
                wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            ),
            gather: pipe(&gather_mod, "gather"),
            fill: pipe(&fill_mod, "fill"),
            restrict: pipe(&restrict_mod, "restrict_down"),
            gather_bind: Vec::new(),
            fill_bind: Vec::new(),
            restrict_bind: Vec::new(),
        }
    }

    /// Paramètres et groupes de liaison des couplages (les tampons de
    /// Spalart–Allmaras ou d'efforts peuvent avoir été réalloués).
    fn prepare_links(&mut self) {
        for (k, link) in self.links.iter_mut().enumerate() {
            let (p, c) = (&self.levels[k], &self.levels[k + 1]);
            let sa_on = (p.sa_active() && c.sa_active()) as u32;
            for (a, buf) in link.params.iter().enumerate() {
                let lp = LinkParams {
                    pnx: p.cfg.nx as u32,
                    pny: p.cfg.ny as u32,
                    pn: (p.cfg.nx * p.cfg.ny) as u32,
                    cnx: c.cfg.nx as u32,
                    cny: c.cfg.ny as u32,
                    cn: (c.cfg.nx * c.cfg.ny) as u32,
                    ox: link.ox as u32,
                    oy: link.oy as u32,
                    ng: link.ng,
                    margin: link.margin as u32,
                    sa_on,
                    nrx: link.nrx,
                    alpha: 0.5 * a as f32,
                    nu_p: p.cfg.nu as f32,
                    nu_c: c.cfg.nu as f32,
                    nry: link.nry,
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
            // Après un pas parent, l'indice courant `cur` désigne le nouvel
            // état ; l'ancien (f[1 − cur]) fournit les populations tirées.
            link.gather_bind = (0..2)
                .map(|cur| {
                    group(
                        &link.gather,
                        vec![(0, &link.params[0]), (1, &link.gm), (2, &p.f[1 - cur]), (3, &p.sa[cur])],
                    )
                })
                .collect();
            link.fill_bind = (0..4)
                .map(|i| {
                    let (cur, a) = (i / 2, i % 2);
                    group(
                        &link.fill,
                        vec![(0, &link.params[a]), (1, &link.gm), (2, &c.f[cur]), (3, &c.sa[cur])],
                    )
                })
                .collect();
            link.restrict_bind = (0..4)
                .map(|i| {
                    let (pc, cc) = (i / 2, i % 2);
                    group(
                        &link.restrict,
                        vec![
                            (0, &link.params[0]),
                            (2, &c.f[1 - cc]),
                            (3, &c.f[cc]),
                            (4, &c.kind),
                            (5, &c.sa[cc]),
                            (6, &p.f[pc]),
                            (7, &p.sa[pc]),
                            (8, &p.kind),
                        ],
                    )
                })
                .collect();
        }
    }

    /// Un pas de la grille k, suivi récursivement de deux pas de chaque grille plus fine.
    fn encode(levels: &mut [GpuSolver], links: &[Link], pass: &mut wgpu::ComputePass, k: usize, slot: u32) {
        let lev = &mut levels[k];
        let (nx, ny) = (lev.cfg.nx as u32, lev.cfg.ny as u32);
        pass.set_pipeline(&lev.step_pipe);
        pass.set_bind_group(0, &lev.bind[lev.cur], &[slot]);
        pass.dispatch_workgroups(nx.div_ceil(64), ny, 1);
        if k == 0 {
            pass.set_pipeline(&lev.bound_pipe);
            pass.dispatch_workgroups(ny.div_ceil(64), 1, 1);
        }
        if lev.sa_active() {
            pass.set_pipeline(&lev.sa_pipe);
            pass.set_bind_group(0, &lev.sa_bind[lev.cur], &[slot]);
            pass.dispatch_workgroups(nx.div_ceil(64), ny, 1);
        }
        lev.cur = 1 - lev.cur;
        let Some(link) = links.get(k) else { return };
        pass.set_pipeline(&link.gather);
        pass.set_bind_group(0, &link.gather_bind[levels[k].cur], &[]);
        pass.dispatch_workgroups(link.ng.div_ceil(64), 1, 1);
        for a in 0..2 {
            pass.set_pipeline(&link.fill);
            pass.set_bind_group(0, &link.fill_bind[levels[k + 1].cur * 2 + a], &[]);
            pass.dispatch_workgroups(link.ng.div_ceil(64), 1, 1);
            Self::encode(levels, links, pass, k + 1, slot);
        }
        pass.set_pipeline(&link.restrict);
        pass.set_bind_group(0, &link.restrict_bind[levels[k].cur * 2 + levels[k + 1].cur], &[]);
        pass.dispatch_workgroups(link.nrx.div_ceil(64), link.nry, 1);
    }

    /// `n` pas de la veine (2ᵏ pas de la grille k).
    pub fn run(&mut self, n: usize, ramp: &Ramp) {
        self.prepare_links();
        let cost: usize = self
            .levels
            .iter()
            .enumerate()
            .map(|(k, l)| (l.cfg.nx * l.cfg.ny) << k)
            .sum();
        let per_submit = (30_000_000 / cost.max(1)).clamp(1, BATCH);
        let device = self.levels[0].device.clone();
        let queue = self.levels[0].queue.clone();
        if self.fresh {
            // Moments des mailles fantômes aux temps t et t + 1 = état initial.
            let mut enc = device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                for (k, link) in self.links.iter().enumerate() {
                    for _ in 0..2 {
                        pass.set_pipeline(&link.gather);
                        pass.set_bind_group(0, &link.gather_bind[self.levels[k].cur], &[]);
                        pass.dispatch_workgroups(link.ng.div_ceil(64), 1, 1);
                    }
                }
            }
            queue.submit([enc.finish()]);
            self.fresh = false;
        }
        let mut left = n;
        while left > 0 {
            let kc = left.min(per_submit);
            let t0 = self.levels[0].time;
            for (k, lev) in self.levels.iter().enumerate() {
                let mut bytes = vec![0u8; SLOT as usize * kc];
                let scale = (1u64 << k) as f32;
                for s in 0..kc {
                    let (u, a) = ramp.at(t0 + s as u64);
                    let p = lev.step_params(u, a / scale);
                    let o = s * SLOT as usize;
                    bytes[o..o + std::mem::size_of_val(&p)].copy_from_slice(bytemuck::bytes_of(&p));
                }
                queue.write_buffer(&lev.params, 0, &bytes);
            }
            let mut enc = device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                for s in 0..kc {
                    Self::encode(&mut self.levels, &self.links, &mut pass, 0, (s as u64 * SLOT) as u32);
                }
            }
            queue.submit([enc.finish()]);
            for (k, lev) in self.levels.iter_mut().enumerate() {
                let steps = (kc as u64) << k;
                lev.time += steps;
                lev.force_steps += steps;
            }
            left -= kc;
        }
    }

    fn depth(&self) -> usize {
        self.levels.len() - 1
    }

    /// Efforts de la grille la plus fine, ramenés aux unités de la veine :
    /// force par pas ÷ 2ᴷ (longueur de référence 2ᴷ fois plus grande),
    /// moment ÷ 4ᴷ, pas comptés en pas de la veine.
    fn scale_forces(&self, t: u64, f: ForceSum) -> (u64, ForceSum) {
        let k = self.depth() as i32;
        let s = 2f64.powi(k);
        (
            t >> k,
            ForceSum {
                fx: f.fx / (s * s),
                fy: f.fy / (s * s),
                mz: f.mz / (s * s * s),
                steps: f.steps >> k,
            },
        )
    }

    pub fn take_forces(&mut self) -> ForceSum {
        let fine = self.levels.last_mut().unwrap();
        let t = fine.time;
        let f = fine.take_forces();
        self.scale_forces(t, f).1
    }

    pub fn poll_forces(&mut self) -> Vec<(u64, ForceSum)> {
        let v = self.levels.last_mut().unwrap().poll_forces();
        v.into_iter().map(|(t, f)| self.scale_forces(t, f)).collect()
    }

    pub fn macros(&mut self, out: &mut Macros) {
        self.levels[0].macros(out);
    }

    /// ρ, u de la grille k.
    pub fn level_macros(&mut self, k: usize, out: &mut Macros) {
        self.levels[k].macros(out);
    }

    pub fn time(&self) -> u64 {
        self.levels[0].time
    }

    pub fn tau0(&self) -> f64 {
        self.levels[0].tau0()
    }

    pub fn les_active(&self) -> bool {
        self.levels[0].les_active()
    }

    pub fn is_finite(&mut self) -> bool {
        self.levels.iter_mut().all(|l| l.is_finite())
    }

    /// Sous-maille active ou non partout comme dans la veine.
    fn sync_les(&mut self) {
        let on = self.levels[0].les_active();
        for k in 1..self.levels.len() {
            let nu = self.levels[k].cfg.nu;
            self.levels[k].cfg.les = if on { LesMode::On } else { LesMode::Off };
            self.levels[k].set_viscosity(nu);
        }
    }

    pub fn set_viscosity(&mut self, nu: f64) {
        self.levels[0].set_viscosity(nu);
        for k in 1..self.levels.len() {
            self.levels[k].cfg.nu = nu * (1u64 << k) as f64;
        }
        self.sync_les();
    }

    /// Longueur d'amortissement pariétal de la sous-maille de chaque grille.
    pub fn set_les_wall_lengths(&mut self, a: impl Fn(usize) -> f64) {
        for (k, l) in self.levels.iter_mut().enumerate() {
            l.set_les_wall_length(a(k));
        }
    }

    pub fn set_turbulence(&mut self, on: bool) {
        for l in &mut self.levels {
            l.set_turbulence(on);
        }
    }

    pub fn sa_active(&self) -> bool {
        self.levels[0].sa_active()
    }

    pub fn sa_state(&mut self) -> Vec<[f32; 4]> {
        self.levels[0].sa_state()
    }

    pub fn update_laminar_zone(&mut self) {
        for l in &mut self.levels {
            l.update_laminar_zone();
        }
    }

    pub fn geometry(&self) -> &Geometry {
        &self.levels[0].geo
    }

    pub fn geometry_mut(&mut self) -> &mut Geometry {
        &mut self.levels[0].geo
    }

    /// Place la maquette dans chaque grille (même point de référence).
    pub fn set_body(&mut self, body: &Body) {
        self.bodies.clear();
        for k in 0..self.levels.len() {
            let b = if k == 0 {
                body.clone()
            } else {
                let s = (1u64 << k) as f64;
                let [gx, gy] = self.origin[k];
                Body::new(
                    body.kind,
                    body.params(),
                    body.l * s,
                    (body.ref_x - gx) * s,
                    (body.ref_y - gy) * s,
                    self.u_in as f64,
                )
            };
            self.levels[k].set_body(&b);
            self.bodies.push(b);
        }
    }

    pub fn reset(&mut self, u0: f32) {
        for l in &mut self.levels {
            l.reset(u0);
        }
        self.fresh = true;
    }

    pub fn adapter_name(&self) -> String {
        let l = self.levels.last().unwrap();
        format!(
            "{} · {} grille{} fine{} ×{} ({}×{})",
            self.levels[0].adapter_name(),
            self.depth(),
            if self.depth() > 1 { "s" } else { "" },
            if self.depth() > 1 { "s" } else { "" },
            1 << self.depth(),
            l.cfg.nx,
            l.cfg.ny
        )
    }
}
