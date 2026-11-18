// Geometry + pose metrics for the docking toy. No DOM access, so it can be tested with node
// (see tests/test_metrics.js).
(function (root) {
  "use strict";

  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm = (a) => Math.sqrt(dot(a, a));
  const normalize = (a) => { const n = norm(a); return n > 0 ? scale(a, 1 / n) : [0, 0, 0]; };
  const dist2 = (a, b) => { const d = sub(a, b); return dot(d, d); };

  function centroid(xyz) {
    const c = [0, 0, 0];
    xyz.forEach((p) => { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; });
    return scale(c, 1 / xyz.length);
  }

  // 3x3 matrices as arrays of rows.
  const identity = () => [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const matVec = (m, v) => [dot(m[0], v), dot(m[1], v), dot(m[2], v)];
  function matMul(a, b) {
    const bt = [[b[0][0], b[1][0], b[2][0]], [b[0][1], b[1][1], b[2][1]], [b[0][2], b[1][2], b[2][2]]];
    return a.map((row) => bt.map((col) => dot(row, col)));
  }

  // Rotation matrix for angle (radians) about a unit axis (Rodrigues).
  function axisAngle(axis, angle) {
    const [x, y, z] = normalize(axis);
    const c = Math.cos(angle), s = Math.sin(angle), t = 1 - c;
    return [
      [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
      [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
      [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ];
  }

  // Uniformly random rotation (Shoemake 1992, random unit quaternion). `rand` defaults to Math.random.
  function randomRotation(rand = Math.random) {
    const u1 = rand(), u2 = rand(), u3 = rand();
    const a = Math.sqrt(1 - u1), b = Math.sqrt(u1);
    const [w, x, y, z] = [a * Math.sin(2 * Math.PI * u2), a * Math.cos(2 * Math.PI * u2),
      b * Math.sin(2 * Math.PI * u3), b * Math.cos(2 * Math.PI * u3)];
    return [
      [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
      [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
      [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ];
  }

  // Re-orthonormalise a rotation matrix (accumulated float error after many drags).
  function orthonormalize(m) {
    const r0 = normalize(m[0]);
    const r1 = normalize(sub(m[1], scale(r0, dot(r0, m[1]))));
    return [r0, r1, cross(r0, r1)];
  }

  // Pose = rotation R about the reference centroid c, then translation t: x' = R (x - c) + c + t.
  function applyPose(refXyz, pose, c = centroid(refXyz)) {
    return refXyz.map((p) => add(add(matVec(pose.R, sub(p, c)), c), pose.t));
  }

  // Plain RMSD with fixed atom correspondence i <-> i (no superposition: docking RMSD is
  // measured in the receptor frame).
  function rmsd(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += dist2(a[i], b[i]);
    return Math.sqrt(s / a.length);
  }

  // Symmetry-corrected RMSD: minimum over graph automorphisms `perms` (perm[i] = atom matched to i).
  function symmetryRmsd(pred, ref, perms) {
    let best = Infinity, bestPerm = null;
    (perms && perms.length ? perms : [pred.map((_, i) => i)]).forEach((perm) => {
      let s = 0;
      for (let i = 0; i < pred.length; i++) s += dist2(pred[i], ref[perm[i]]);
      const r = Math.sqrt(s / pred.length);
      if (r < best) { best = r; bestPerm = perm; }
    });
    return { rmsd: best, perm: bestPerm };
  }

  // Number of ligand-receptor heavy-atom pairs closer than `cutoff` (crude clash count).
  function countClashes(ligXyz, recXyz, cutoff) {
    const c2 = cutoff * cutoff;
    let n = 0;
    for (const p of ligXyz) for (const q of recXyz) if (dist2(p, q) < c2) n++;
    return n;
  }

  const api = {
    sub, add, scale, dot, cross, norm, normalize, centroid, identity, matVec, matMul,
    axisAngle, randomRotation, orthonormalize, applyPose, rmsd, symmetryRmsd, countClashes,
    dist: (a, b) => Math.sqrt(dist2(a, b)),
  };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Metrics = api;
})(typeof window !== "undefined" ? window : this);
