// Pocket-finding algorithms from the lecture (Protein Bioinformatics II, slides 22-29), written
// to be read by students. No DOM access, so tests/test_pockets.js can run them with node.
//
//   gridPockets   grid-based, in the spirit of LIGSITE (Hendlich et al. 1997): a grid point
//                 outside the atoms is "buried" in a direction if protein lies on both sides of
//                 it along that line (a protein-solvent-protein event). Count the events over 7
//                 directions, keep well-buried points, cluster, rank.
//   alphaPockets  grid-free, in the spirit of fpocket (Le Guilloux et al. 2009): alpha spheres
//                 (empty spheres through 4 atom centres, from a Delaunay triangulation) are
//                 filtered by radius, clustered, ranked.
//
// Both are simplified teaching versions, not the published programs.
(function (root) {
  "use strict";

  // van der Waals radii (A), Bondi 1964; anything else 1.8.
  const VDW = { C: 1.7, N: 1.55, O: 1.52, S: 1.8, P: 1.8, F: 1.47, Cl: 1.75, Br: 1.85, I: 1.98, Zn: 1.39, Ca: 1.97, Mg: 1.73, Se: 1.9 };
  const vdw = (e) => VDW[e] || 1.8;

  // The 7 scan directions of LIGSITE: the 3 axes and the 4 cube diagonals. Each is written with
  // its first non-zero component positive, which lets one ascending x/y/z loop visit every
  // point after its predecessor along the line.
  const DIRS = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1]];

  // ------------------------------------------------------------------ parsing helpers
  function atomsFromPdb(pdb) {
    const xyz = [], elem = [], resi = [], chain = [], resn = [], name = [];
    for (const l of pdb.split("\n")) {
      if (!l.startsWith("ATOM") && !l.startsWith("HETATM")) continue;
      xyz.push(+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54));
      const e = (l.slice(76, 78).trim() || l.slice(12, 14).trim());
      elem.push(e.charAt(0).toUpperCase() + e.slice(1).toLowerCase());
      resi.push(+l.slice(22, 26)); chain.push(l[21]); resn.push(l.slice(17, 20).trim()); name.push(l.slice(12, 16).trim());
    }
    return { xyz: Float64Array.from(xyz), elem, resi, chain, resn, name, n: elem.length };
  }

  function centroid(points) {
    const c = [0, 0, 0];
    for (const p of points) { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; }
    return c.map((v) => v / points.length);
  }
  const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);

  // Spatial hash for fixed-radius neighbour queries.
  function hashPoints(pts, cell) {
    const map = new Map();
    const key = (i, j, k) => (i * 73856093) ^ (j * 19349663) ^ (k * 83492791);
    pts.forEach((p, idx) => {
      const kk = key(Math.floor(p[0] / cell), Math.floor(p[1] / cell), Math.floor(p[2] / cell));
      if (!map.has(kk)) map.set(kk, []);
      map.get(kk).push(idx);
    });
    return {
      near(p, r) {
        const out = [];
        const i0 = Math.floor((p[0] - r) / cell), i1 = Math.floor((p[0] + r) / cell);
        const j0 = Math.floor((p[1] - r) / cell), j1 = Math.floor((p[1] + r) / cell);
        const k0 = Math.floor((p[2] - r) / cell), k1 = Math.floor((p[2] + r) / cell);
        for (let i = i0; i <= i1; i++) for (let j = j0; j <= j1; j++) for (let k = k0; k <= k1; k++) {
          const b = map.get(key(i, j, k));
          if (b) for (const q of b) if (dist(pts[q], p) <= r) out.push(q);
        }
        return out;
      },
    };
  }

  // ------------------------------------------------------------------ grid-based (LIGSITE-like)
  // params: spacing (A), minBuried (events out of 7), maxScan (A, how far to look for protein
  // along a line; Infinity = to the edge of the box), minVolume (A^3), probe (A added to every
  // atom radius, so narrow gaps between atoms count as protein).
  function gridPockets(atoms, params) {
    const h = params.spacing, pad = 4;
    let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (let a = 0; a < atoms.n; a++) for (let d = 0; d < 3; d++) {
      lo[d] = Math.min(lo[d], atoms.xyz[3 * a + d]); hi[d] = Math.max(hi[d], atoms.xyz[3 * a + d]);
    }
    lo = lo.map((v) => v - pad); hi = hi.map((v) => v + pad);
    const nx = Math.ceil((hi[0] - lo[0]) / h) + 1, ny = Math.ceil((hi[1] - lo[1]) / h) + 1, nz = Math.ceil((hi[2] - lo[2]) / h) + 1;
    const N = nx * ny * nz, idx = (i, j, k) => (i * ny + j) * nz + k;

    // Steps 1-2: place the grid, mark points inside atoms as protein.
    const occ = new Uint8Array(N);
    for (let a = 0; a < atoms.n; a++) {
      const r = vdw(atoms.elem[a]) + (params.probe || 0), r2 = r * r;
      const x = atoms.xyz[3 * a], y = atoms.xyz[3 * a + 1], z = atoms.xyz[3 * a + 2];
      const i0 = Math.max(0, Math.floor((x - r - lo[0]) / h)), i1 = Math.min(nx - 1, Math.ceil((x + r - lo[0]) / h));
      const j0 = Math.max(0, Math.floor((y - r - lo[1]) / h)), j1 = Math.min(ny - 1, Math.ceil((y + r - lo[1]) / h));
      const k0 = Math.max(0, Math.floor((z - r - lo[2]) / h)), k1 = Math.min(nz - 1, Math.ceil((z + r - lo[2]) / h));
      for (let i = i0; i <= i1; i++) { const dx = lo[0] + i * h - x;
        for (let j = j0; j <= j1; j++) { const dy = lo[1] + j * h - y;
          for (let k = k0; k <= k1; k++) { const dz = lo[2] + k * h - z;
            if (dx * dx + dy * dy + dz * dz <= r2) occ[idx(i, j, k)] = 1; } } }
    }

    // Step 3: count protein-solvent-protein events. For each direction, `before` holds how many
    // steps back the nearest protein point is (65535 = none), `after` the same looking forward.
    // Because every direction's first non-zero component is +1, the ascending loop always reaches
    // a point after the point one step behind it, and the descending loop after the one ahead.
    const maxSteps = Number.isFinite(params.maxScan) ? Math.round(params.maxScan / h) : 65534;
    const buried = new Uint8Array(N);
    const before = new Uint16Array(N), after = new Uint16Array(N);
    for (const [dx, dy, dz] of DIRS) {
      const step = (dx * ny + dy) * nz + dz;
      const inside = (i, j, k) => i >= 0 && j >= 0 && k >= 0 && i < nx && j < ny && k < nz;
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) for (let k = 0; k < nz; k++) {
        const p = idx(i, j, k);
        if (!inside(i - dx, j - dy, k - dz)) { before[p] = 65535; continue; }
        const q = p - step;
        before[p] = occ[q] ? 1 : (before[q] >= 65534 ? 65535 : before[q] + 1);
      }
      for (let i = nx - 1; i >= 0; i--) for (let j = ny - 1; j >= 0; j--) for (let k = nz - 1; k >= 0; k--) {
        const p = idx(i, j, k);
        if (!inside(i + dx, j + dy, k + dz)) { after[p] = 65535; continue; }
        const q = p + step;
        after[p] = occ[q] ? 1 : (after[q] >= 65534 ? 65535 : after[q] + 1);
      }
      for (let p = 0; p < N; p++) if (!occ[p] && before[p] <= maxSteps && after[p] <= maxSteps) buried[p]++;
    }

    // Optional: every solvent point buried in at least `showBuried` directions, for display.
    let buriedPoints = null;
    if (params.showBuried) {
      buriedPoints = [];
      for (let p = 0; p < N; p++) if (!occ[p] && buried[p] >= params.showBuried) {
        const i = Math.floor(p / (ny * nz)), j = Math.floor(p / nz) % ny, k = p % nz;
        buriedPoints.push([lo[0] + i * h, lo[1] + j * h, lo[2] + k * h, buried[p]]);
      }
    }

    // Step 4: keep buried points; step 5: cluster (26-neighbour connectivity), score, rank.
    const keep = new Uint8Array(N);
    for (let p = 0; p < N; p++) if (!occ[p] && buried[p] >= params.minBuried) keep[p] = 1;
    const label = new Int32Array(N).fill(-1);
    const pockets = [];
    const stack = [];
    for (let s = 0; s < N; s++) {
      if (!keep[s] || label[s] >= 0) continue;
      const pts = [], bur = [];
      label[s] = pockets.length; stack.push(s);
      while (stack.length) {
        const p = stack.pop();
        const i = Math.floor(p / (ny * nz)), j = Math.floor(p / nz) % ny, k = p % nz;
        pts.push([lo[0] + i * h, lo[1] + j * h, lo[2] + k * h]); bur.push(buried[p]);
        for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let c = -1; c <= 1; c++) {
          const ii = i + a, jj = j + b, kk = k + c;
          if (ii < 0 || jj < 0 || kk < 0 || ii >= nx || jj >= ny || kk >= nz) continue;
          const q = idx(ii, jj, kk);
          if (keep[q] && label[q] < 0) { label[q] = label[s]; stack.push(q); }
        }
      }
      const volume = pts.length * h * h * h;
      if (volume < params.minVolume) { pockets.push(null); continue; }
      const meanBuried = bur.reduce((x, y) => x + y, 0) / bur.length;
      pockets.push({ points: pts, volume, meanBuried, centre: centroid(pts) });
    }
    const found = pockets.filter(Boolean);
    const score = params.score === "buriedness" ? (p) => p.volume * (p.meanBuried - params.minBuried + 1) : (p) => p.volume;
    found.forEach((p) => { p.score = score(p); });
    found.sort((a, b) => b.score - a.score);
    return { pockets: found, buriedPoints, grid: { nx, ny, nz, h, points: N } };
  }

  // ------------------------------------------------------------------ alpha spheres (fpocket-like)
  // alpha: {x,y,z,r,v} as written by build_data.py (v = 4 atom indices per sphere).
  // params: rMin, rMax (A), linkDist (A, single-linkage distance between sphere centres),
  // minSpheres, score: "size" | "apolar" | "buriedness".
  function alphaPockets(atoms, alpha, params) {
    const kept = [];
    for (let s = 0; s < alpha.r.length; s++) {
      const r = alpha.r[s];
      if (r < params.rMin || r > params.rMax) continue;
      let carbons = 0;
      for (let t = 0; t < 4; t++) if (atoms.elem[alpha.v[4 * s + t]] === "C") carbons++;
      kept.push({ c: [alpha.x[s], alpha.y[s], alpha.z[s]], r, apolar: carbons >= 3 });
    }
    // single-linkage clustering of the sphere centres
    const centres = kept.map((k) => k.c);
    const hash = hashPoints(centres, params.linkDist);
    const label = new Int32Array(kept.length).fill(-1);
    const clusters = [];
    for (let s = 0; s < kept.length; s++) {
      if (label[s] >= 0) continue;
      const members = [s]; label[s] = clusters.length;
      for (let m = 0; m < members.length; m++) {
        for (const q of hash.near(centres[members[m]], params.linkDist)) {
          if (label[q] < 0) { label[q] = label[s]; members.push(q); }
        }
      }
      clusters.push(members);
    }
    // burial proxy: protein atoms within 8 A of the pocket centre
    const atomPts = [];
    for (let a = 0; a < atoms.n; a++) atomPts.push([atoms.xyz[3 * a], atoms.xyz[3 * a + 1], atoms.xyz[3 * a + 2]]);
    const atomHash = hashPoints(atomPts, 8);
    const pockets = clusters.filter((m) => m.length >= params.minSpheres).map((m) => {
      const spheres = m.map((i) => kept[i]);
      const centre = centroid(spheres.map((s) => s.c));
      const nApolar = spheres.filter((s) => s.apolar).length;
      return { spheres, centre, n: spheres.length, apolarFrac: nApolar / spheres.length,
               atomsNear: atomHash.near(centre, 8).length };
    });
    const score = {
      size: (p) => p.n,
      apolar: (p) => p.n * p.apolarFrac,
      buriedness: (p) => p.n * p.atomsNear / 100,
    }[params.score || "size"];
    pockets.forEach((p) => { p.score = score(p); });
    pockets.sort((a, b) => b.score - a.score);
    return { pockets, kept: kept.length };
  }

  // ------------------------------------------------------------------ evaluation against the ligand
  // DCC: distance from the pocket centre to the ligand centre (success if <= 4 A, the usual
  // criterion in pocket-finding benchmarks). Coverage: fraction of ligand atoms that have a pocket
  // point (grid point or alpha-sphere centre) within `covR` A.
  function evaluate(pocket, ligandXyz, covR = 3) {
    const pts = pocket.points || pocket.spheres.map((s) => s.c);
    const dcc = dist(pocket.centre, centroid(ligandXyz));
    let covered = 0;
    for (const l of ligandXyz) if (pts.some((p) => dist(p, l) <= covR)) covered++;
    return { dcc, coverage: covered / ligandXyz.length };
  }

  // Defaults chosen with tests/explore_pockets.js on the ten Pocket-hunt proteins (see README).
  const GRID_DEFAULTS = { spacing: 1.0, probe: 1.2, maxScan: 8, minBuried: 6, minVolume: 30, score: "size" };
  // rMin/rMax are fpocket's default alpha-sphere radius limits (3.4 and 6.2 A).
  const ALPHA_DEFAULTS = { rMin: 3.4, rMax: 6.2, linkDist: 2.0, minSpheres: 10, score: "size" };
  // A pocket counts as found when its centre is within 4 A of the ligand centre (DCC).
  const DCC_SUCCESS_A = 4;

  const api = { GRID_DEFAULTS, ALPHA_DEFAULTS, DCC_SUCCESS_A, VDW, vdw, DIRS, atomsFromPdb, centroid, dist, hashPoints, gridPockets, alphaPockets, evaluate };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Pockets = api;
})(typeof window !== "undefined" ? window : this);
