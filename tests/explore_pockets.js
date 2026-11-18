// Parameter exploration: rank of the true pocket for each protein under a few settings.
const P = require("../docs/js/pockets.js"), { load } = require("./load_data.js");
const index = load("index.js");
const settings = JSON.parse(process.argv[2] || "{}");
const method = settings.method || "grid";
let hits = 0;
for (const L of index.levels) {
  const d = load(`prot_${L.id}.js`);
  const atoms = P.atomsFromPdb(d.receptor_pdb);
  const lig = d.ligand.atoms.map((a) => a.xyz);
  const t0 = Date.now();
  const res = method === "grid" ? P.gridPockets(atoms, settings) : P.alphaPockets(atoms, d.alpha, settings);
  const ev = res.pockets.slice(0, 10).map((p) => P.evaluate(p, lig));
  const rDcc = ev.findIndex((e) => e.dcc <= 4) + 1;
  const rCov = ev.findIndex((e) => e.coverage >= 0.5) + 1;
  if (rCov === 1) hits++;
  const top = res.pockets.slice(0, 3).map((p, i) => `${(p.volume || p.n).toFixed(0)}/dcc${ev[i].dcc.toFixed(1)}/cov${ev[i].coverage.toFixed(2)}`).join("  ");
  console.log(`${L.id.padEnd(22)} n=${String(res.pockets.length).padStart(3)} rankDCC=${rDcc || "-"} rankCov=${rCov || "-"} ${Date.now() - t0}ms | ${top}`);
}
console.log("top-1 by coverage:", hits, "/", index.levels.length);
