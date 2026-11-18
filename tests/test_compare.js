// Checks docs/js/compare.js against the fixtures from tests/make_compare_fixtures.py.
//   python tests/make_compare_fixtures.py && node tests/test_compare.js
const fs = require("fs"), path = require("path");
const C = require("../docs/js/compare.js"), { load } = require("./load_data.js");
const OUT = path.join(__dirname, "out");
const expect = JSON.parse(fs.readFileSync(path.join(OUT, "compare_expect.json"), "utf8"));
let fail = 0;
for (const [id, cases] of Object.entries(expect)) {
  const cof = load(`cofold_${id}.js`), prot = load(`prot_${id}.js`);
  const ref = { ca_by_seq: cof.ca_by_seq, pocket_seq: cof.pocket_seq, ligand: prot.ligand };
  for (const [file, want] of Object.entries(cases)) {
    const res = C.compareToCrystal(C.parseStructure(fs.readFileSync(path.join(OUT, `${id}_${file}`), "utf8")), ref);
    let ok;
    if (want === null) ok = res.ligandRmsd === null && res.warnings.some((w) => w.includes("atom names"));
    else ok = Math.abs(res.ligandRmsd - want) < 0.01 && res.caRmsd < 0.01;
    if (!ok) fail++;
    console.log(`${ok ? "ok  " : "FAIL"} ${id.padEnd(14)} ${file.padEnd(12)} residues ${res.nResidues} pocket ${res.nPocket} ` +
      `CA ${res.caRmsd.toFixed(3)} ligand ${res.ligandRmsd === null ? "n/a" : res.ligandRmsd.toFixed(3)} (want ${want === null ? "warning" : want}) ${res.warnings.join("; ").slice(0, 60)}`);
  }
}
console.log(fail ? `${fail} FAILED` : "all passed");
process.exit(fail ? 1 : 0);
