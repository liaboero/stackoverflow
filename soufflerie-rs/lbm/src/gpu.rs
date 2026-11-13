//! Moteur GPU natif (wgpu : Vulkan, Metal, DirectX 12) : mêmes équations que
//! le moteur processeur. Un pas = une passe de calcul (propagation « pull » et
//! collision fusionnées) puis une passe pour l'entrée et la sortie. Les
//! paramètres variables d'un pas à l'autre (mise en vitesse) passent par un
//! décalage dynamique du tampon uniforme : des centaines de pas partent en une
//! seule soumission. Les efforts s'accumulent par maille pariétale et sont
//! relus en une fois. Avec Spalart–Allmaras, une troisième passe transporte ν̃.

use std::sync::mpsc;

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::cpu::{Collision, ForceSum, LesMode, Macros, Ramp, SolverConfig};
use crate::geometry::Geometry;
use crate::shapes::Body;

const STEP_WGSL: &str = include_str!("d2q9.wgsl");
const MACRO_WGSL: &str = include_str!("macros.wgsl");
const SA_WGSL: &str = include_str!("sa.wgsl");
/// Pas par soumission (emplacements de paramètres).
const BATCH: usize = 256;
const SLOT: u64 = 256;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Params {
    nx: u32,
    ny: u32,
    n: u32,
    flags: u32,
    u_in: f32,
    accel: f32,
    c18: f32,
    spin: f32,
    rot_x: f32,
    rot_y: f32,
    ref_x: f32,
    ref_y: f32,
    nu: f32,
    pad: [f32; 3],
}

pub struct GpuSolver {
    pub cfg: SolverConfig,
    pub geo: Geometry,
    device: wgpu::Device,
    queue: wgpu::Queue,
    adapter_name: String,
    f: [wgpu::Buffer; 2],
    kind: wgpu::Buffer,
    wq: wgpu::Buffer,
    aux: wgpu::Buffer,
    /// Par maille : (amortissement de la sous-maille, distance, Cf laminaire, angle de la normale).
    cell: wgpu::Buffer,
    /// Spalart–Allmaras : (u_x, u_y, u_τ, ν̃), deux tampons alternés.
    sa: [wgpu::Buffer; 2],
    force: wgpu::Buffer,
    force_read: wgpu::Buffer,
    params: wgpu::Buffer,
    macro_buf: wgpu::Buffer,
    macro_read: wgpu::Buffer,
    _macro_params: wgpu::Buffer,
    layout: wgpu::BindGroupLayout,
    sa_layout: wgpu::BindGroupLayout,
    step_pipe: wgpu::ComputePipeline,
    bound_pipe: wgpu::ComputePipeline,
    sa_pipe: wgpu::ComputePipeline,
    macro_pipe: wgpu::ComputePipeline,
    bind: Vec<wgpu::BindGroup>,
    sa_bind: Vec<wgpu::BindGroup>,
    macro_bind: [wgpu::BindGroup; 2],
    sa_active: bool,
    cur: usize,
    pub time: u64,
    tau0: f64,
    les_active: bool,
    force_steps: u64,
    wall_count: usize,
    finite: bool,
    ref_x: f64,
    ref_y: f64,
    spin: f64,
}

/// Entrées d'un groupe de liaison, numérotées dans l'ordre.
fn entries(res: Vec<wgpu::BindingResource<'_>>) -> Vec<wgpu::BindGroupEntry<'_>> {
    res.into_iter()
        .enumerate()
        .map(|(k, resource)| wgpu::BindGroupEntry {
            binding: k as u32,
            resource,
        })
        .collect()
}

fn storage(read_only: bool, binding: u32) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

impl GpuSolver {
    pub fn new(cfg: SolverConfig) -> Result<GpuSolver, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle_from_env());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        }))
        .map_err(|e| format!("aucun adaptateur GPU : {e}"))?;
        let info = adapter.get_info();
        let adapter_name = format!("{} ({:?})", info.name, info.backend);
        let mut limits = wgpu::Limits::default();
        let al = adapter.limits();
        limits.max_storage_buffer_binding_size = al.max_storage_buffer_binding_size;
        limits.max_buffer_size = al.max_buffer_size;
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("soufflerie"),
            required_limits: limits,
            ..Default::default()
        }))
        .map_err(|e| format!("GPU indisponible : {e}"))?;
        let (nx, ny) = (cfg.nx, cfg.ny);
        let n = nx * ny;
        let mk = |label: &str, size: u64, usage: wgpu::BufferUsages| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: size.max(16),
                usage,
                mapped_at_creation: false,
            })
        };
        let st = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC;
        let f = [mk("f0", 36 * n as u64, st), mk("f1", 36 * n as u64, st)];
        let kind = mk("kind", 4 * n as u64, st);
        let wq = mk("wq", 36 * n as u64, st);
        let aux = mk("aux", 4 * (nx + ny) as u64, st);
        let cell = mk("cell", 16 * n as u64, st);
        // Alloués à la première activation de Spalart–Allmaras.
        let sa = [mk("sa0", 16, st), mk("sa1", 16, st)];
        let force = mk("force", 16, st);
        let force_read = mk(
            "force_read",
            16,
            wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        );
        let params = mk(
            "params",
            SLOT * BATCH as u64,
            wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        );
        let macro_buf = mk("macros", 12 * n as u64, st);
        let macro_read = mk(
            "macros_read",
            12 * n as u64,
            wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        );
        let macro_params = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("macro_params"),
            contents: bytemuck::cast_slice(&[nx as u32, ny as u32, n as u32, 0u32]),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("d2q9"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: true,
                        min_binding_size: wgpu::BufferSize::new(std::mem::size_of::<Params>() as u64),
                    },
                    count: None,
                },
                storage(true, 1),
                storage(false, 2),
                storage(true, 3),
                storage(true, 4),
                storage(true, 5),
                storage(true, 6),
                storage(false, 7),
                storage(false, 8),
            ],
        });
        let sa_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("sa"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: true,
                        min_binding_size: wgpu::BufferSize::new(std::mem::size_of::<Params>() as u64),
                    },
                    count: None,
                },
                storage(true, 1),
                storage(true, 2),
                storage(true, 3),
                storage(true, 4),
                storage(false, 5),
            ],
        });
        let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("d2q9"),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("d2q9"),
            source: wgpu::ShaderSource::Wgsl(STEP_WGSL.into()),
        });
        let pipe = |entry: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(entry),
                layout: Some(&pl),
                module: &module,
                entry_point: Some(entry),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        let step_pipe = pipe("step");
        let bound_pipe = pipe("boundary");
        let sa_pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("sa"),
            bind_group_layouts: &[Some(&sa_layout)],
            immediate_size: 0,
        });
        let sa_mod = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("sa"),
            source: wgpu::ShaderSource::Wgsl(SA_WGSL.into()),
        });
        let sa_pipe = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("sa_step"),
            layout: Some(&sa_pl),
            module: &sa_mod,
            entry_point: Some("sa_step"),
            compilation_options: Default::default(),
            cache: None,
        });
        let mmod = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("macros"),
            source: wgpu::ShaderSource::Wgsl(MACRO_WGSL.into()),
        });
        let macro_pipe = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("macros"),
            layout: None,
            module: &mmod,
            entry_point: Some("macros"),
            compilation_options: Default::default(),
            cache: None,
        });
        let geo = Geometry::new(nx, ny, cfg.les_wall_a);
        let mut s = GpuSolver {
            bind: Vec::new(),
            sa_bind: Vec::new(),
            sa_active: false,
            macro_bind: [
                Self::macro_group(&device, &macro_pipe, &macro_params, &f[0], &kind, &macro_buf),
                Self::macro_group(&device, &macro_pipe, &macro_params, &f[1], &kind, &macro_buf),
            ],
            cfg,
            geo,
            device,
            queue,
            adapter_name,
            f,
            kind,
            wq,
            aux,
            cell,
            sa,
            force,
            force_read,
            params,
            macro_buf,
            macro_read,
            _macro_params: macro_params,
            layout,
            sa_layout,
            step_pipe,
            bound_pipe,
            sa_pipe,
            macro_pipe,
            cur: 0,
            time: 0,
            tau0: 0.5,
            les_active: false,
            force_steps: 0,
            wall_count: 0,
            finite: true,
            ref_x: nx as f64 / 4.0,
            ref_y: (ny as f64 - 1.0) / 2.0,
            spin: 0.0,
        };
        s.set_viscosity(s.cfg.nu);
        s.upload_geometry();
        s.rebuild_binds();
        s.reset(s.cfg.u_in);
        Ok(s)
    }

    /// Groupes de liaison des passes (à refaire quand un tampon change).
    fn rebuild_binds(&mut self) {
        let d = &self.device;
        let params = wgpu::BindingResource::Buffer(wgpu::BufferBinding {
            buffer: &self.params,
            offset: 0,
            size: wgpu::BufferSize::new(std::mem::size_of::<Params>() as u64),
        });
        self.bind = (0..2)
            .map(|c| {
                d.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("d2q9"),
                    layout: &self.layout,
                    entries: &entries(vec![
                        params.clone(),
                        self.f[c].as_entire_binding(),
                        self.f[1 - c].as_entire_binding(),
                        self.kind.as_entire_binding(),
                        self.wq.as_entire_binding(),
                        self.aux.as_entire_binding(),
                        self.cell.as_entire_binding(),
                        self.force.as_entire_binding(),
                        self.sa[c].as_entire_binding(),
                    ]),
                })
            })
            .collect();
        self.sa_bind = (0..2)
            .map(|c| {
                d.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("sa"),
                    layout: &self.sa_layout,
                    entries: &entries(vec![
                        params.clone(),
                        self.kind.as_entire_binding(),
                        self.aux.as_entire_binding(),
                        self.cell.as_entire_binding(),
                        self.sa[c].as_entire_binding(),
                        self.sa[1 - c].as_entire_binding(),
                    ]),
                })
            })
            .collect();
    }

    fn macro_group(
        device: &wgpu::Device,
        pipe: &wgpu::ComputePipeline,
        params: &wgpu::Buffer,
        f: &wgpu::Buffer,
        kind: &wgpu::Buffer,
        out: &wgpu::Buffer,
    ) -> wgpu::BindGroup {
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("macros"),
            layout: &pipe.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: f.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: kind.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: out.as_entire_binding(),
                },
            ],
        })
    }

    pub fn adapter_name(&self) -> &str {
        &self.adapter_name
    }

    pub fn tau0(&self) -> f64 {
        self.tau0
    }

    pub fn les_active(&self) -> bool {
        self.les_active
    }

    pub fn set_viscosity(&mut self, nu: f64) {
        self.cfg.nu = nu;
        self.tau0 = 3.0 * nu + 0.5;
        self.les_active = match self.cfg.les {
            LesMode::On => true,
            LesMode::Off => false,
            LesMode::Auto => self.tau0 < 0.51,
        };
        self.upload_aux();
    }

    pub fn set_les_wall_length(&mut self, a: f64) {
        self.geo.set_les_wall_length(a);
        self.upload_cell();
    }

    /// Active ou coupe Spalart–Allmaras avec loi de paroi.
    pub fn set_turbulence(&mut self, on: bool) {
        let was = self.sa_active;
        self.sa_active = on;
        let size = 16 * (self.cfg.nx * self.cfg.ny) as u64;
        if on && self.sa[0].size() < size {
            let st = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC;
            let mk = |label| {
                self.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some(label),
                    size,
                    usage: st,
                    mapped_at_creation: false,
                })
            };
            self.sa = [mk("sa0"), mk("sa1")];
            self.rebuild_binds();
        }
        if on && !was {
            self.init_sa();
        }
    }

    pub fn sa_active(&self) -> bool {
        self.sa_active
    }

    /// Zone laminaire (`geo.lam`) modifiée : renvoie les données par maille.
    pub fn update_laminar_zone(&mut self) {
        self.upload_cell();
    }

    /// ν̃ = 3ν au large, nul dans le corps ; vitesses et u_τ nuls.
    fn init_sa(&mut self) {
        let n = self.cfg.nx * self.cfg.ny;
        if self.sa[0].size() < 16 * n as u64 {
            return;
        }
        let v = crate::sa::nu_tilde_inf(self.cfg.nu);
        let init: Vec<[f32; 4]> = self
            .geo
            .kind
            .iter()
            .map(|&k| [0.0, 0.0, 0.0, if k == crate::geometry::SOLID { 0.0 } else { v }])
            .collect();
        for b in &self.sa {
            self.queue.write_buffer(b, 0, bytemuck::cast_slice(&init));
        }
    }

    /// ν̃ et u_τ (lecture synchrone ; vide si Spalart–Allmaras est inactif).
    pub fn sa_state(&mut self) -> Vec<[f32; 4]> {
        let n = self.cfg.nx * self.cfg.ny;
        if !self.sa_active || self.sa[0].size() < 16 * n as u64 {
            return Vec::new();
        }
        let size = 16 * n as u64;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("sa_read"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let v = self.read_buffer(&self.sa[self.cur], &staging, size);
        v.chunks_exact(4).map(|c| [c[0], c[1], c[2], c[3]]).collect()
    }

    /// Amortissement de la sous-maille, distance, Cf laminaire et angle de la normale.
    fn upload_cell(&mut self) {
        let g = &self.geo;
        let cell: Vec<[f32; 4]> = (0..g.n)
            .map(|i| {
                [
                    g.les_damp[i],
                    g.dist[i],
                    g.lam[i],
                    g.normal[2 * i + 1].atan2(g.normal[2 * i]),
                ]
            })
            .collect();
        self.queue.write_buffer(&self.cell, 0, bytemuck::cast_slice(&cell));
    }

    /// τ par colonne (éponge) puis coefficient des bandes absorbantes par ligne.
    fn upload_aux(&mut self) {
        let (nx, ny) = (self.cfg.nx, self.cfg.ny);
        let s0 = (nx as f64 * (1.0 - self.cfg.sponge_fraction)).round() as usize;
        let tau_sponge = self.tau0.max(0.9);
        let len = (nx - 1).saturating_sub(s0).max(1) as f64;
        let absorb = (ny as f64 / 8.0).round();
        let mut aux: Vec<f32> = (0..nx)
            .map(|x| {
                let mut tau = self.tau0;
                if x > s0 {
                    let s = (x - s0) as f64 / len;
                    tau += (tau_sponge - self.tau0) * s * s;
                }
                tau as f32
            })
            .collect();
        aux.extend((0..ny).map(|y| {
            let d = y.min(ny - 1 - y) as f64;
            let s = if absorb > 0.0 { 1.0 - d / absorb } else { 0.0 };
            if s > 0.0 {
                (0.05 * s * s) as f32
            } else {
                0.0
            }
        }));
        self.queue.write_buffer(&self.aux, 0, bytemuck::cast_slice(&aux));
    }

    /// Envoie la géométrie : type de maille et indice de maille pariétale, q, amortissement.
    fn upload_geometry(&mut self) {
        let mut kind: Vec<u32> = self.geo.kind.iter().map(|&k| k as u32).collect();
        for (w, &idx) in self.geo.wall_cells.iter().enumerate() {
            kind[idx] |= ((w + 1) as u32) << 2;
        }
        self.queue.write_buffer(&self.kind, 0, bytemuck::cast_slice(&kind));
        self.queue
            .write_buffer(&self.wq, 0, bytemuck::cast_slice(&self.geo.wall_q));
        self.upload_cell();
        let wc = self.geo.wall_cells.len();
        if wc != self.wall_count || self.force.size() < (16 * wc.max(1)) as u64 {
            let size = (16 * wc.max(1)) as u64;
            self.force = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("force"),
                size,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            self.force_read = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("force_read"),
                size,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            self.wall_count = wc;
            if !self.bind.is_empty() {
                self.rebuild_binds();
            }
        }
        self.wall_count = wc;
    }

    pub fn set_body(&mut self, body: &Body) {
        self.ref_x = body.ref_x;
        self.ref_y = body.ref_y;
        self.spin = body.spin;
        self.geo.build(body);
        self.upload_geometry();
    }

    /// Fluide à la vitesse u0, avec la petite perturbation du moteur processeur.
    pub fn reset(&mut self, u0: f32) {
        let (nx, ny) = (self.cfg.nx, self.cfg.ny);
        let n = nx * ny;
        let kick = 0.15 * u0.max(self.cfg.u_in).max(0.05);
        let r = ny as f64 / 8.0;
        let px = self.ref_x + 1.5 * r;
        let py = self.ref_y + 0.5 * r;
        let mut f = vec![0f32; 9 * n];
        for y in 0..ny {
            for x in 0..nx {
                let idx = y * nx + x;
                let eq = if self.geo.kind[idx] == crate::geometry::SOLID {
                    crate::d2q9::equilibrium(1.0, 0.0, 0.0)
                } else {
                    let d2 = (x as f64 - px).powi(2) + (y as f64 - py).powi(2);
                    crate::d2q9::equilibrium(1.0, u0, (kick as f64 * (-d2 / (r * r)).exp()) as f32)
                };
                for i in 0..9 {
                    f[i * n + idx] = eq[i];
                }
            }
        }
        for b in &self.f {
            self.queue.write_buffer(b, 0, bytemuck::cast_slice(&f));
        }
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.clear_buffer(&self.force, 0, None);
        self.queue.submit([enc.finish()]);
        self.cur = 0;
        self.time = 0;
        self.force_steps = 0;
        self.finite = true;
        self.init_sa();
    }

    fn flags(&self) -> u32 {
        (self.cfg.collision == Collision::Rr) as u32
            | ((self.les_active as u32) << 1)
            | ((self.cfg.bouzidi as u32) << 2)
            | ((self.sa_active as u32) << 3)
    }

    /// `n` pas ; la vitesse d'entrée suit `ramp`.
    pub fn run(&mut self, n: usize, ramp: &Ramp) {
        let (nx, ny) = (self.cfg.nx as u32, self.cfg.ny as u32);
        let mut left = n;
        while left > 0 {
            let k = left.min(BATCH);
            let mut bytes = vec![0u8; SLOT as usize * k];
            for s in 0..k {
                let (u, a) = ramp.at(self.time + s as u64);
                let p = Params {
                    nx,
                    ny,
                    n: nx * ny,
                    flags: self.flags(),
                    u_in: u,
                    accel: a,
                    c18: 18.0 * self.cfg.cs * self.cfg.cs,
                    spin: self.spin as f32,
                    rot_x: self.ref_x as f32,
                    rot_y: self.ref_y as f32,
                    ref_x: self.ref_x as f32,
                    ref_y: self.ref_y as f32,
                    nu: self.cfg.nu as f32,
                    pad: [0.0; 3],
                };
                let o = s * SLOT as usize;
                bytes[o..o + std::mem::size_of::<Params>()].copy_from_slice(bytemuck::bytes_of(&p));
            }
            self.queue.write_buffer(&self.params, 0, &bytes);
            let mut enc = self.device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                for s in 0..k {
                    let off = (s as u64 * SLOT) as u32;
                    pass.set_pipeline(&self.step_pipe);
                    pass.set_bind_group(0, &self.bind[self.cur], &[off]);
                    pass.dispatch_workgroups(nx.div_ceil(64), ny, 1);
                    pass.set_pipeline(&self.bound_pipe);
                    pass.dispatch_workgroups(ny.div_ceil(64), 1, 1);
                    if self.sa_active {
                        pass.set_pipeline(&self.sa_pipe);
                        pass.set_bind_group(0, &self.sa_bind[self.cur], &[off]);
                        pass.dispatch_workgroups(nx.div_ceil(64), ny, 1);
                    }
                    self.cur = 1 - self.cur;
                }
            }
            self.queue.submit([enc.finish()]);
            self.time += k as u64;
            self.force_steps += k as u64;
            left -= k;
        }
    }

    fn wait(&self) {
        let _ = self.device.poll(wgpu::PollType::wait_indefinitely());
    }

    fn read_buffer(&self, src: &wgpu::Buffer, staging: &wgpu::Buffer, size: u64) -> Vec<f32> {
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(src, 0, staging, 0, size);
        self.queue.submit([enc.finish()]);
        let (tx, rx) = mpsc::channel();
        staging.map_async(wgpu::MapMode::Read, 0..size, move |r| {
            let _ = tx.send(r);
        });
        self.wait();
        let ok = rx.recv().map(|r| r.is_ok()).unwrap_or(false);
        let out = if ok {
            let view = staging.get_mapped_range(0..size).expect("lecture GPU");
            bytemuck::cast_slice::<u8, f32>(&view).to_vec()
        } else {
            Vec::new()
        };
        staging.unmap();
        out
    }

    /// Efforts cumulés depuis le dernier appel (lecture synchrone), puis remise à zéro.
    pub fn take_forces(&mut self) -> ForceSum {
        if self.force_steps == 0 || self.wall_count == 0 {
            return ForceSum {
                steps: std::mem::take(&mut self.force_steps),
                ..Default::default()
            };
        }
        let size = 16 * self.wall_count as u64;
        let v = self.read_buffer(&self.force, &self.force_read, size);
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.clear_buffer(&self.force, 0, None);
        self.queue.submit([enc.finish()]);
        let (mut fx, mut fy, mut mz) = (0.0f64, 0.0f64, 0.0f64);
        for c in v.chunks_exact(4) {
            fx += c[0] as f64;
            fy += c[1] as f64;
            mz += c[2] as f64;
        }
        self.finite = fx.is_finite() && fy.is_finite();
        ForceSum {
            fx,
            fy,
            mz,
            steps: std::mem::take(&mut self.force_steps),
        }
    }

    /// ρ, u (lecture synchrone).
    pub fn macros(&mut self, out: &mut Macros) {
        let (nx, ny) = (self.cfg.nx, self.cfg.ny);
        let n = nx * ny;
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.macro_pipe);
            pass.set_bind_group(0, &self.macro_bind[self.cur], &[]);
            pass.dispatch_workgroups((n as u32).div_ceil(64), 1, 1);
        }
        self.queue.submit([enc.finish()]);
        let v = self.read_buffer(&self.macro_buf, &self.macro_read, 12 * n as u64);
        out.nx = nx;
        out.ny = ny;
        if v.len() == 3 * n {
            out.rho = v[..n].to_vec();
            out.ux = v[n..2 * n].to_vec();
            out.uy = v[2 * n..].to_vec();
            let mid = out.rho[n / 2];
            if !(mid > 0.3 && mid < 3.0) {
                self.finite = false;
            }
        }
    }

    pub fn is_finite(&mut self) -> bool {
        self.finite
    }
}
