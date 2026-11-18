// Pocket algorithms on the ten Pocket-hunt proteins. Run: node tests/test_pockets.js
//  1. Writes tests/out/pocket_ranks.md: rank of the true pocket (DCC <= 4 A) for each protein under
//     the default settings of both methods and the scoring variants. The README table is this file.
//  2. Checks that every Chapter 2 challenge can be completed, with the settings named below.
//  3. Checks the morph endpoints: no pocket at the ligand site at the start, one at the end.
const fs = require("fs"), path = require("path");
const P = require("../docs/js/pockets.js"), { load } = require("./load_data.js");

let failures = 0;
const check = (cond, msg) => { console.log(`${cond ? "ok  " : "FAIL"} ${msg}`); if (!cond) failures++; };
const index = load("index.js");
const prot = {};
for (const L of index.levels) {
  const d = load(`prot_${L.id}.js`);
  prot[L.id] = { d, title: L.title, atoms: P.atomsFromPdb(d.receptor_pdb), lig: d.ligand.atoms.map((a) => a.xyz) };
}
const rank = (id, method, params) => {
  const p = prot[id];
  const res = method === "grid" ? P.gridPockets(p.atoms, params) : P.alphaPockets(p.atoms, p.d.alpha, params);
  return res.pockets.findIndex((pk) => P.evaluate(pk, p.lig).dcc <= P.DCC_SUCCESS_A) + 1;
};

// 1. table
const variants = [
  ["Grid, volume", "grid", P.GRID_DEFAULTS],
  ["Grid, volume × burial", "grid", { ...P.GRID_DEFAULTS, score: "buriedness" }],
  ["Alpha, size", "alpha", P.ALPHA_DEFAULTS],
  ["Alpha, hydrophobicity", "alpha", { ...P.ALPHA_DEFAULTS, score: "apolar" }],
  ["Alpha, burial", "alpha", { ...P.ALPHA_DEFAULTS, score: "buriedness" }],
];
const rows = Object.keys(prot).map((id) => [prot[id].title, ...variants.map(([, m, prm]) => rank(id, m, prm))]);
const fmt = (r) => (r === 0 ? "–" : `#${r}`);
let md = `| Protein | ${variants.map((v) => v[0]).join(" | ")} |\n|---|${variants.map(() => "---").join("|")}|\n`;
rows.forEach((r) => { md += `| ${r[0]} | ${r.slice(1).map(fmt).join(" | ")} |\n`; });
md += `| **Rank 1** | ${variants.map((_, j) => `${rows.filter((r) => r[j + 1] === 1).length}/10`).join(" | ")} |\n`;
fs.mkdirSync(path.join(__dirname, "out"), { recursive: true });
fs.writeFileSync(path.join(__dirname, "out", "pocket_ranks.md"), md);
console.log(md);

// 2. challenges (ids as in docs/js/ch_algo.js)
check(rank("streptavidin", "grid", P.GRID_DEFAULTS) === 1, "a1: grid ranks biotin's pocket #1 on streptavidin");
check(rank("t4_lysozyme", "grid", P.GRID_DEFAULTS) > 1, "a3 premise: grid by volume does not rank the T4L cavity #1");
check(rank("t4_lysozyme", "alpha", { ...P.ALPHA_DEFAULTS, score: "apolar" }) === 1, "a3: alpha + hydrophobicity ranks the T4L cavity #1");
check(rank("carbonic_anhydrase", "alpha", P.ALPHA_DEFAULTS) > 3, "a4 premise: alpha defaults miss the CA zinc pocket");
const caFix = rank("carbonic_anhydrase", "alpha", { ...P.ALPHA_DEFAULTS, rMin: 3.0 });
check(caFix >= 1 && caFix <= 3, `a4: alpha with rMin 3.0 puts the CA pocket at #${caFix}`);
check(rank("ptp1b", "alpha", P.ALPHA_DEFAULTS) === 1 && rank("ptp1b", "grid", P.GRID_DEFAULTS) === 0, "a5: alpha finds the PTP1B allosteric site at #1, grid does not find it");
check(rank("trypsin", "alpha", { ...P.ALPHA_DEFAULTS, score: "apolar" }) > rank("trypsin", "alpha", P.ALPHA_DEFAULTS), "a6: hydrophobicity score drops trypsin's S1 pocket");
const best = { ...P.ALPHA_DEFAULTS, rMin: 3.2, score: "apolar" };
const top1 = Object.keys(prot).filter((id) => rank(id, "alpha", best) === 1).length;
check(top1 >= 8, `a7: alpha, hydrophobicity, rMin 3.2 reaches ${top1}/10 at rank 1`);

// 3. morphs
for (const m of index.morphs) {
  const d = load(`morph_${m.id}.js`);
  const base = P.atomsFromPdb(d.pdb), la = P.atomsFromPdb(d.ligand_pdb), lig = [];
  for (let i = 0; i < la.n; i++) lig.push([la.xyz[3 * i], la.xyz[3 * i + 1], la.xyz[3 * i + 2]]);
  const r = (xyz) => P.gridPockets({ ...base, xyz: Float64Array.from(xyz) }, P.GRID_DEFAULTS).pockets
    .findIndex((pk) => P.evaluate(pk, lig).dcc <= P.DCC_SUCCESS_A) + 1;
  const [r0, r1] = [r(d.start_xyz), r(d.end_xyz)];
  check(r0 === 0 && r1 === 1, `morph ${m.id}: pocket at the ligand site start ${fmt(r0)}, end ${fmt(r1)}`);
}
process.exit(failures ? 1 : 0);
