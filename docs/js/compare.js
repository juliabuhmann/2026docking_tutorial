// Compare a predicted protein-ligand complex with the crystal structure. No DOM access, so
// tests/test_compare.js can run it with node.
//
//   parseStructure(text)    PDB or mmCIF -> { protein: Map(resNumber -> CA xyz), ligand: [{name, elem, xyz}] }
//   superpose(P, Q)         least-squares rotation + translation taking P onto Q (Horn 1987, quaternions)
//   compareToCrystal(pred, ref)
//       ref: { ca_by_seq: {seqIndex: xyz}, pocket_seq: [seqIndex], ligand: {atoms, automorphisms} }
//       Superposes the prediction onto the crystal on the CA atoms of the pocket residues (as in
//       PoseBusters-style "pocket-aligned ligand RMSD") and reports the ligand RMSD, matching
//       atoms by name and taking the best symmetry-equivalent matching. Also reports the CA RMSD
//       of the whole chain after a global superposition.
(function (root) {
  "use strict";

  // ------------------------------------------------------------------ parsing
  function parsePdb(text) {
    const protein = new Map(), ligand = [];
    for (const l of text.split(/\r?\n/)) {
      if (l.startsWith("ENDMDL")) break;
      const rec = l.slice(0, 6).trim();
      if (rec !== "ATOM" && rec !== "HETATM") continue;
      const name = l.slice(12, 16).trim(), resn = l.slice(17, 20).trim();
      const xyz = [+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)];
      const elem = (l.slice(76, 78).trim() || name.replace(/[^A-Za-z]/g, "").slice(0, 1)).toUpperCase();
      if (resn === "HOH" || elem === "H") continue;
      if (rec === "ATOM") { if (name === "CA") protein.set(+l.slice(22, 26), xyz); }
      else ligand.push({ name, elem, resn, xyz });
    }
    return { protein, ligand };
  }

  // mmCIF tokens: whitespace separated, with '...' or "..." quoting.
  function cifTokens(line) {
    const out = [];
    const re = /'([^']*)'|"([^"]*)"|(\S+)/g;
    let m;
    while ((m = re.exec(line))) out.push(m[1] !== undefined ? m[1] : m[2] !== undefined ? m[2] : m[3]);
    return out;
  }
  function parseCif(text) {
    const lines = text.split(/\r?\n/);
    const protein = new Map(), ligand = [];
    for (let i = 0; i < lines.length; i++) {
      if (lines[i].trim() !== "loop_" || !(lines[i + 1] || "").startsWith("_atom_site.")) continue;
      const cols = [];
      let j = i + 1;
      while (j < lines.length && lines[j].startsWith("_atom_site.")) cols.push(lines[j++].trim().slice(11));
      const c = (n) => cols.indexOf(n);
      const iGroup = c("group_PDB"), iAtom = c("label_atom_id") >= 0 ? c("label_atom_id") : c("auth_atom_id");
      const iResn = c("label_comp_id") >= 0 ? c("label_comp_id") : c("auth_comp_id");
      const iSeq = c("label_seq_id"), iAuth = c("auth_seq_id"), iEl = c("type_symbol");
      const iX = c("Cartn_x"), iY = c("Cartn_y"), iZ = c("Cartn_z"), iModel = c("pdbx_PDB_model_num");
      let firstModel = null;
      for (; j < lines.length; j++) {
        const s = lines[j].trim();
        if (!s || s.startsWith("#") || s.startsWith("loop_") || s.startsWith("_")) break;
        const t = cifTokens(s);
        if (iModel >= 0) { if (firstModel === null) firstModel = t[iModel]; else if (t[iModel] !== firstModel) break; }
        const name = t[iAtom], resn = t[iResn], elem = (t[iEl] || "").toUpperCase();
        if (resn === "HOH" || elem === "H") continue;
        const xyz = [+t[iX], +t[iY], +t[iZ]];
        if (t[iGroup] === "ATOM") {
          const seq = iSeq >= 0 && t[iSeq] !== "." ? +t[iSeq] : +t[iAuth];
          if (name === "CA") protein.set(seq, xyz);
        } else ligand.push({ name, elem, resn, xyz });
      }
      break;
    }
    return { protein, ligand };
  }

  function parseStructure(text) {
    return /^\s*(data_|loop_|_atom_site)/m.test(text.slice(0, 20000)) && text.includes("_atom_site.") ? parseCif(text) : parsePdb(text);
  }

  // ------------------------------------------------------------------ superposition (Horn's quaternion method)
  function centroid(P) {
    const c = [0, 0, 0];
    P.forEach((p) => { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; });
    return c.map((v) => v / P.length);
  }
  // Largest eigenpair of a symmetric 4x4 matrix by Jacobi rotations.
  function eigSym4(A) {
    const a = A.map((r) => r.slice()), V = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    for (let sweep = 0; sweep < 60; sweep++) {
      let off = 0;
      for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) off += a[p][q] * a[p][q];
      if (off < 1e-20) break;
      for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) {
        if (Math.abs(a[p][q]) < 1e-30) continue;
        const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const cs = 1 / Math.sqrt(t * t + 1), sn = t * cs;
        for (let k = 0; k < 4; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = cs * x - sn * y; a[k][q] = sn * x + cs * y; }
        for (let k = 0; k < 4; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = cs * x - sn * y; a[q][k] = sn * x + cs * y; }
        for (let k = 0; k < 4; k++) { const x = V[k][p], y = V[k][q]; V[k][p] = cs * x - sn * y; V[k][q] = sn * x + cs * y; }
      }
    }
    let best = 0;
    for (let i = 1; i < 4; i++) if (a[i][i] > a[best][best]) best = i;
    return [V[0][best], V[1][best], V[2][best], V[3][best]];
  }
  // Returns {R (3x3 rows), t} so that R p + t best matches q, for paired point lists P, Q.
  function superpose(P, Q) {
    const cp = centroid(P), cq = centroid(Q);
    const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < P.length; i++) {
      const p = [P[i][0] - cp[0], P[i][1] - cp[1], P[i][2] - cp[2]], q = [Q[i][0] - cq[0], Q[i][1] - cq[1], Q[i][2] - cq[2]];
      for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) S[a][b] += p[a] * q[b];
    }
    const [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = S;
    const N = [
      [xx + yy + zz, yz - zy, zx - xz, xy - yx],
      [yz - zy, xx - yy - zz, xy + yx, zx + xz],
      [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
      [xy - yx, zx + xz, yz + zy, -xx - yy + zz],
    ];
    const [w, x, y, z] = eigSym4(N);
    const R = [
      [w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
      [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
      [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z],
    ];
    const Rcp = apply(R, [0, 0, 0], cp);
    return { R, t: [cq[0] - Rcp[0], cq[1] - Rcp[1], cq[2] - Rcp[2]] };
  }
  function apply(R, t, p) {
    return [R[0][0] * p[0] + R[0][1] * p[1] + R[0][2] * p[2] + t[0],
      R[1][0] * p[0] + R[1][1] * p[1] + R[1][2] * p[2] + t[1],
      R[2][0] * p[0] + R[2][1] * p[1] + R[2][2] * p[2] + t[2]];
  }
  function rmsd(P, Q) {
    let s = 0;
    for (let i = 0; i < P.length; i++) s += (P[i][0] - Q[i][0]) ** 2 + (P[i][1] - Q[i][1]) ** 2 + (P[i][2] - Q[i][2]) ** 2;
    return Math.sqrt(s / P.length);
  }

  // ------------------------------------------------------------------ the comparison
  function compareToCrystal(pred, ref) {
    const out = { warnings: [] };
    // CA pairs by sequence position
    const pairs = [];
    for (const [k, xyz] of Object.entries(ref.ca_by_seq)) if (pred.protein.has(+k)) pairs.push([+k, pred.protein.get(+k), xyz]);
    out.nResidues = pairs.length;
    if (pairs.length < 10) throw new Error("fewer than 10 residues could be matched to the crystal structure. Was this protein folded from the sequence given on this page?");
    const g = superpose(pairs.map((p) => p[1]), pairs.map((p) => p[2]));
    out.caRmsd = rmsd(pairs.map((p) => apply(g.R, g.t, p[1])), pairs.map((p) => p[2]));

    const pocket = new Set(ref.pocket_seq);
    const pp = pairs.filter((p) => pocket.has(p[0]));
    out.nPocket = pp.length;
    const s = pp.length >= 3 ? superpose(pp.map((p) => p[1]), pp.map((p) => p[2])) : g;
    if (pp.length < 3) out.warnings.push("too few pocket residues matched; used the global superposition");
    out.pocketCaRmsd = pp.length ? rmsd(pp.map((p) => apply(s.R, s.t, p[1])), pp.map((p) => p[2])) : NaN;
    out.transform = s;

    // ligand atoms matched by name, in crystal order
    const byName = new Map(pred.ligand.map((a) => [a.name, a]));
    const refAtoms = ref.ligand.atoms;
    const matched = refAtoms.map((a) => byName.get(a.name));
    const missing = refAtoms.filter((a, i) => !matched[i]).map((a) => a.name);
    out.nLigand = refAtoms.length - missing.length;
    if (!pred.ligand.length) { out.ligandRmsd = null; out.warnings.push("no ligand (HETATM) atoms in the file"); return out; }
    if (missing.length) {
      out.ligandRmsd = null;
      out.warnings.push(`ligand atom names do not match the crystal (${missing.slice(0, 5).join(", ")}${missing.length > 5 ? ", …" : ""} missing). Fold with the CCD code, not a SMILES, so the atoms keep their names.`);
      return out;
    }
    const predXyz = matched.map((a) => apply(s.R, s.t, a.xyz));
    const refXyz = refAtoms.map((a) => a.xyz);
    let best = Infinity;
    for (const perm of (ref.ligand.automorphisms && ref.ligand.automorphisms.length ? ref.ligand.automorphisms : [refXyz.map((_, i) => i)])) {
      let sum = 0;
      for (let i = 0; i < predXyz.length; i++) { const q = refXyz[perm[i]]; sum += (predXyz[i][0] - q[0]) ** 2 + (predXyz[i][1] - q[1]) ** 2 + (predXyz[i][2] - q[2]) ** 2; }
      best = Math.min(best, Math.sqrt(sum / predXyz.length));
    }
    out.ligandRmsd = best;
    out.predLigandXyz = predXyz;
    return out;
  }

  const api = { parsePdb, parseCif, parseStructure, superpose, apply, rmsd, compareToCrystal };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Compare = api;
})(typeof window !== "undefined" ? window : this);
