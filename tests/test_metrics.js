// Sanity checks for docs/js/metrics.js (ported from the docking_game prototype). Run:
//   node tests/test_metrics.js && python tests/check_rmsd_rdkit.py
// Writes random poses to tests/out/<case>_random_poses.json so check_rmsd_rdkit.py can cross-check
// the symmetry-corrected RMSD against RDKit, for every case of the docking chapter.
const fs = require("fs");
const path = require("path");
const M = require("../docs/js/metrics.js");
const { load } = require("./load_data.js");
const N_RANDOM_POSES = 200;
const SEED = 42;

let state = SEED; // small LCG so the random poses are reproducible
const rand = () => ((state = (state * 1664525 + 1013904223) % 4294967296) / 4294967296);

let failures = 0;
const check = (cond, msg) => { console.log(`${cond ? "ok  " : "FAIL"} ${msg}`); if (!cond) failures++; };

const index = load("index.js");
for (const dc of index.docking.cases) {
  const cs = load(`prot_${dc.level}.js`);
  const ref = cs.ligand.atoms.map((a) => a.xyz);
  const c = M.centroid(ref);
  const id = { R: M.identity(), t: [0, 0, 0] };
  check(M.rmsd(M.applyPose(ref, id), ref) < 1e-9, `${dc.level}: identity pose has RMSD 0`);
  const shifted = M.applyPose(ref, { R: M.identity(), t: [2, 0, 0] });
  check(Math.abs(M.rmsd(shifted, ref) - 2) < 1e-9, `${dc.level}: 2 A translation gives RMSD 2`);
  const R = M.randomRotation(rand);
  check(Math.abs(M.dot(M.cross(R[0], R[1]), R[2]) - 1) < 1e-9, `${dc.level}: random rotation has det +1`);
  const poses = [];
  for (let k = 0; k < N_RANDOM_POSES; k++) {
    const pose = { R: M.randomRotation(rand), t: [(rand() - 0.5) * 6, (rand() - 0.5) * 6, (rand() - 0.5) * 6] };
    const xyz = M.applyPose(ref, pose, c);
    poses.push({ xyz, rmsd: M.rmsd(xyz, ref), sym_rmsd: M.symmetryRmsd(xyz, ref, cs.ligand.automorphisms).rmsd });
  }
  check(poses.every((p) => p.sym_rmsd <= p.rmsd + 1e-12), `${dc.level}: symmetry RMSD <= plain RMSD`);
  const out = path.join(__dirname, "out");
  fs.mkdirSync(out, { recursive: true });
  fs.writeFileSync(path.join(out, `${dc.level}_random_poses.json`), JSON.stringify(poses));
}
process.exit(failures ? 1 : 0);
