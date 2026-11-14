// ρ, ux, uy par maille (trois plans), pour l'affichage.
struct M {
  nx: u32, ny: u32, n: u32, pad: u32,
};
@group(0) @binding(0) var<uniform> P: M;
// Binding 1 (populations f) et fonction ldm : préambule ajouté par gpu.rs.
@group(0) @binding(2) var<storage, read> kind: array<u32>;
@group(0) @binding(3) var<storage, read_write> outm: array<f32>;

@compute @workgroup_size(64, 1, 1)
fn macros(@builtin(global_invocation_id) gid: vec3<u32>) {
  // Une ligne de la veine par groupe y (au plus 65 535 groupes par dimension).
  if (gid.x >= P.nx || gid.y >= P.ny) {
    return;
  }
  let idx = gid.y * P.nx + gid.x;
  if ((kind[idx] & 3u) == 1u) {
    outm[idx] = 1.0;
    outm[P.n + idx] = 0.0;
    outm[2u * P.n + idx] = 0.0;
    return;
  }
  var v: array<f32, 9>;
  for (var i = 0u; i < 9u; i++) {
    v[i] = ldm(i, idx);
  }
  let rho = v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7] + v[8];
  outm[idx] = rho;
  outm[P.n + idx] = (v[1] - v[3] + v[5] - v[6] - v[7] + v[8]) / rho;
  outm[2u * P.n + idx] = (v[2] - v[4] + v[5] + v[6] - v[7] - v[8]) / rho;
}
