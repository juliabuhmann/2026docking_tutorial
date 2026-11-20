// Checks docs/js/diffusion.js: AF3's schedule endpoints, proper rotations, and that the sampler with
// the oracle denoiser ends on the crystal structure.   node tests/test_diffusion.js
const D = require("../docs/js/diffusion.js"), C = require("../docs/js/compare.js"), { load } = require("./load_data.js");
let fail = 0;
const check = (ok, msg) => { if (!ok) fail++; console.log(`${ok ? "ok  " : "FAIL"} ${msg}`); };
const close = (a, b, rel) => Math.abs(a - b) <= rel * Math.abs(b);

const L = D.noiseLevels();
check(L.length === D.AF3.STEPS + 1, `${L.length} noise levels for ${D.AF3.STEPS} steps`);
check(close(L[0], D.AF3.SIGMA_DATA * D.AF3.S_MAX, 1e-9), `first level ${L[0].toFixed(1)} A = sigma_data * s_max`);
check(close(L[L.length - 1], D.AF3.SIGMA_DATA * D.AF3.S_MIN, 1e-9), `last level ${L[L.length - 1].toFixed(4)} A = sigma_data * s_min`);
check(L.every((x, i) => i === 0 || x < L[i - 1]), "levels decrease monotonically");

const rng = D.mulberry32(7);
let worst = 0;
for (let k = 0; k < 200; k++) {
  const R = D.randomRotation(rng);
  for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) worst = Math.max(worst, Math.abs(R[a][0] * R[b][0] + R[a][1] * R[b][1] + R[a][2] * R[b][2] - (a === b ? 1 : 0)));
  const det = R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
  if (Math.abs(det - 1) > 1e-9) worst = Infinity;
}
check(worst < 1e-9, `200 random rotations are orthonormal with det +1 (max error ${worst.toExponential(1)})`);

for (const t of load("index.js").cofold) {
  const prot = load(`prot_${t.id}.js`);
  const xyz = [];
  for (const l of prot.receptor_pdb.split("\n")) if (l.startsWith("ATOM")) xyz.push(+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54));
  prot.ligand.atoms.forEach((a) => xyz.push(...a.xyz));
  const ref = Float64Array.from(xyz), t0 = Date.now();
  const run = D.sample(ref, C, { seed: 1 });
  const ms = Date.now() - t0, last = run.frames[run.frames.length - 1];
  let s = 0;
  for (let i = 0; i < ref.length; i++) s += (last[i] - ref[i]) ** 2;
  const rmsd = Math.sqrt(s / run.n), start = D.spread(run.frames[0]);
  check(rmsd < 0.05, `${t.id}: ${run.n} atoms, final RMSD to the crystal ${rmsd.toFixed(4)} A (${ms} ms)`);
  check(close(start, L[0] * Math.sqrt(3), 0.05), `${t.id}: start cloud spread ${start.toFixed(0)} A ~ sigma_0 * sqrt(3)`);
  check(run.info.every((f, i) => i === 0 || (f.gamma > 0) === (f.sigma > D.AF3.GAMMA_MIN)), `${t.id}: extra noise (gamma) only while sigma > ${D.AF3.GAMMA_MIN} A`);
}
console.log(fail ? `${fail} failed` : "all passed");
process.exit(fail ? 1 : 0);
