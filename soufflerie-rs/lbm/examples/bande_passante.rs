//! Bande passante mémoire atteignable sur le GPU (référence des passes de
//! `banc3d --profil 1`) : copie d'un tableau de 256 Mo par un shader (lecture
//! et écriture contiguës, comme la collision), et lecture seule.
//!
//!   cargo run --release --example bande_passante
use std::time::Instant;

const N: u32 = 64 << 20; // mots de 4 octets : 256 Mo

fn main() {
    let instance = wgpu::Instance::new({
        let mut d = wgpu::InstanceDescriptor::new_without_display_handle();
        d.backends = wgpu::Backends::VULKAN | wgpu::Backends::METAL;
        d
    });
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        ..Default::default()
    }))
    .expect("adaptateur");
    println!("{}", adapter.get_info().name);
    let mut limits = wgpu::Limits::default();
    let al = adapter.limits();
    limits.max_storage_buffer_binding_size = al.max_storage_buffer_binding_size;
    limits.max_buffer_size = al.max_buffer_size;
    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
        required_limits: limits,
        ..Default::default()
    }))
    .expect("périphérique");
    let mk = |usage| {
        device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: N as u64 * 4,
            usage: wgpu::BufferUsages::STORAGE | usage,
            mapped_at_creation: false,
        })
    };
    let a = mk(wgpu::BufferUsages::COPY_DST);
    let b = mk(wgpu::BufferUsages::COPY_SRC);
    for (name, body, bytes) in [
        ("copie (lecture + écriture)", "dst[i] = src[i];", 8.0),
        ("copie vec4", "", 8.0),
        ("lecture seule", "if (src[i] == 12345u) { dst[0] = 1u; }", 4.0),
    ] {
        let vec4 = body.is_empty();
        let ty = if vec4 { "vec4<u32>" } else { "u32" };
        let body = if vec4 { "dst[i] = src[i];" } else { body };
        let n = if vec4 { N / 4 } else { N };
        let src = format!(
            "@group(0) @binding(0) var<storage, read> src: array<{ty}>;
@group(0) @binding(1) var<storage, read_write> dst: array<{ty}>;
@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) g: vec3<u32>, @builtin(num_workgroups) nw: vec3<u32>) {{
  let i = g.x + g.y * nw.x * 256u;
  if (i >= {n}u) {{ return; }}
  {body}
}}"
        );
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: None,
            source: wgpu::ShaderSource::Wgsl(src.into()),
        });
        let pipe = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: None,
            layout: None,
            module: &module,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });
        let bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &pipe.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: a.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: b.as_entire_binding() },
            ],
        });
        let groups = n.div_ceil(256);
        let (gx, gy) = (groups.min(65535), groups.div_ceil(65535));
        let run = |reps: usize| {
            let _ = device.poll(wgpu::PollType::wait_indefinitely());
            let t0 = Instant::now();
            let mut enc = device.create_command_encoder(&Default::default());
            {
                let mut p = enc.begin_compute_pass(&Default::default());
                p.set_pipeline(&pipe);
                p.set_bind_group(0, &bind, &[]);
                for _ in 0..reps {
                    p.dispatch_workgroups(gx, gy, 1);
                }
            }
            queue.submit([enc.finish()]);
            let _ = device.poll(wgpu::PollType::wait_indefinitely());
            t0.elapsed().as_secs_f64() / reps as f64
        };
        run(3);
        let t = (0..3).map(|_| run(20)).fold(f64::INFINITY, f64::min);
        println!("  {name:28} {:6.1} Go/s", N as f64 * bytes / t / 1e9);
    }
}
