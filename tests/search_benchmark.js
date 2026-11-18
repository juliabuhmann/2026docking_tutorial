// Search simple parameter grids for the best top-1 count (DCC <= 4 A) over the ten proteins.
// Used to set a reachable target for the Chapter 2 benchmark challenge.
const P = require("../docs/js/pockets.js"), { load } = require("./load_data.js");
const index = load("index.js");
const prots = index.levels.map((L) => { const d = load(`prot_${L.id}.js`); return { id: L.id, d, atoms: P.atomsFromPdb(d.receptor_pdb), lig: d.ligand.atoms.map((a) => a.xyz) }; });
const top1 = (fn) => prots.filter((p) => { const r = fn(p); return r.pockets.length && P.evaluate(r.pockets[0], p.lig).dcc <= P.DCC_SUCCESS_A; }).map((p) => p.id);
let best = [];
for (const score of ["size", "apolar", "buriedness"]) for (const rMin of [2.8, 3.0, 3.2, 3.4, 3.6]) for (const rMax of [5.5, 6.2, 7, 8])
  for (const linkDist of [1.6, 2.0, 2.5, 3.0]) for (const minSpheres of [5, 10, 20]) {
    const prm = { rMin, rMax, linkDist, minSpheres, score };
    const ids = top1((p) => P.alphaPockets(p.atoms, p.d.alpha, prm));
    best.push({ n: ids.length, method: "alpha", prm, ids });
  }
for (const score of ["size", "buriedness"]) for (const probe of [0.8, 1.2, 1.6]) for (const minBuried of [4, 5, 6, 7]) for (const maxScan of [6, 8, 10, 12]) {
  const prm = { ...P.GRID_DEFAULTS, probe, minBuried, maxScan, score };
  const ids = top1((p) => P.gridPockets(p.atoms, prm));
  best.push({ n: ids.length, method: "grid", prm, ids });
}
best.sort((a, b) => b.n - a.n);
for (const b of best.slice(0, 8)) console.log(b.n, b.method, JSON.stringify(b.prm), "missed:", index.levels.map((L) => L.id).filter((i) => !b.ids.includes(i)).join(","));
