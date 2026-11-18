// Chapter 2: the two pocket-finding algorithms from the lecture, live, on the same ten proteins.
(function () {
  "use strict";
  const { $, esc, el, toast, segmented, slider, checkbox, vec } = PW;
  const RANK_COLORS = ["#ffd166", "#5eead4", "#b69cff", "#ff8fd0", "#60a5fa", "#ffb454", "#7ee787", "#ff8f8f"];

  const TASKS = [
    { id: "a1", text: "Run the <b>grid</b> method on <b>streptavidin</b>. Is biotin's pocket ranked #1?",
      check: (r) => r.method === "grid" && r.protein === "streptavidin" && r.trueRank === 1,
      answer: "Yes. Biotin's pocket is deep and enclosed, so its grid points are buried in many directions, and it is the biggest pocket on this small protein." },
    { id: "a2", text: "Grid method on <b>HIV protease</b>: lower <i>buried directions</i> to 3 or less. What happens to the pockets?",
      check: (r) => r.method === "grid" && r.protein === "hiv_protease" && r.params.minBuried <= 3,
      answer: "Shallow dents on the surface now pass the burial test as well, so pockets grow, merge and multiply. The threshold trades sensitivity for specificity, like the score cut-off of a PSSM." },
    { id: "a3", text: "<b>T4 lysozyme</b>: by size, the active-site cleft wins. Find a method or score that ranks the buried <b>benzene cavity #1</b>.",
      check: (r) => r.protein === "t4_lysozyme" && r.trueRank === 1,
      answer: "The engineered cavity is small but lined almost entirely by carbon atoms. Alpha-sphere scores that reward hydrophobicity or burial prefer it over the larger, more open and more polar cleft. Size is not everything." },
    { id: "a4", text: "<b>Carbonic anhydrase</b>: with fpocket's default radii the zinc pocket ranks around #20. Get it into the <b>top 3</b> with alpha spheres.",
      check: (r) => r.method === "alpha" && r.protein === "carbonic_anhydrase" && r.trueRank >= 1 && r.trueRank <= 3,
      answer: "The bottom of the cone, next to the zinc, is narrow: its alpha spheres are smaller than 3.4 Å, so the default throws them away and the pocket falls apart into small pieces. A minimum radius of about 3.0 Å keeps them (slide 26)." },
    { id: "a5", text: "<b>PTP1B</b>: one method ranks the allosteric site #1, the other does not find it at all. Which is which?",
      check: (r) => r.protein === "ptp1b" && r.method === "alpha" && r.trueRank === 1,
      answer: "Alpha spheres find it; the grid with these settings does not. The allosteric site is a shallow groove between two helices: not buried enough for the grid's burial threshold, but lined by enough atoms to hold a cluster of medium-sized alpha spheres." },
    { id: "a6", text: "<b>Trypsin</b>: switch the alpha-sphere score to <b>hydrophobicity</b>. The S1 pocket drops. Why?",
      check: (r) => r.protein === "trypsin" && r.method === "alpha" && r.params.score === "apolar",
      answer: "S1 is made to bind a positive arginine: Asp189, Ser190 and backbone oxygens line it. A hydrophobicity score assumes drug pockets are greasy, which is often true but not here. Every descriptor carries an assumption." },
    { id: "a7", text: "<b>Benchmark</b> your settings on all ten proteins (button below the table). Can you reach <b>8 of 10 at rank 1</b>?",
      check: (r) => r.benchmarkTop1 >= 8,
      answer: "One way: alpha spheres, hydrophobicity score, smallest sphere about 3.2 Å. No setting here reaches 10/10: a score that wins on greasy pockets loses on polar ones such as trypsin's S1. And because you tuned on these same ten proteins, 8/10 is an optimistic estimate. A fair test needs proteins you did not tune on (slide 29: positives, negatives and a held-out set)." },
  ];

  let S = null;

  function defaultsFor(method) { return { ...(method === "grid" ? Pockets.GRID_DEFAULTS : Pockets.ALPHA_DEFAULTS) }; }

  // ------------------------------------------------------------------ left panel
  function renderLeft() {
    const left = $("#left");
    left.innerHTML = `<h3>Protein</h3><select id="algo-prot" style="width:100%;padding:7px;border-radius:10px;background:var(--panel2);color:var(--text);border:1px solid var(--line)"></select>
      <h3>Method</h3><div id="algo-method"></div><p class="small muted" id="algo-method-why"></p>
      <h3>Parameters</h3><div id="algo-params"></div>
      <h3>Score used for ranking</h3><div id="algo-score"></div>
      <h3>Show</h3><div id="algo-show"></div>
      <div class="row"><button id="algo-defaults" class="ghost">Reset to defaults</button></div>`;
    const sel = $("#algo-prot");
    S.index.levels.forEach((L) => sel.appendChild(el(`<option value="${L.id}">${esc(L.title)}</option>`)));
    sel.value = S.protId;
    sel.addEventListener("change", () => openProtein(sel.value));
    $("#algo-method").appendChild(segmented([{ value: "grid", label: "Grid", title: "LIGSITE-like" }, { value: "alpha", label: "Alpha spheres", title: "fpocket-like" }], S.method, (m) => {
      S.method = m; S.params = defaultsFor(m); renderParams(); run();
    }));
    $("#algo-defaults").addEventListener("click", () => { S.params = defaultsFor(S.method); renderParams(); run(); });
    renderParams();
  }

  function renderParams() {
    const box = $("#algo-params"), P = S.params;
    box.innerHTML = "";
    $("#algo-method-why").innerHTML = S.method === "grid"
      ? "Lay a grid over the protein. A point outside the atoms is <i>buried</i> in a direction if protein lies on both sides of it along that line. Keep well-buried points, cluster them, rank the clusters (slide 23)."
      : "Find every empty sphere touching four atom centres (Delaunay triangulation). Keep the medium-sized ones, cluster them, rank the clusters (slides 24 to 27).";
    const add = (o) => box.appendChild(slider({ ...o, onChange: (v) => { P[o.key] = v; run(); } }));
    if (S.method === "grid") {
      add({ key: "spacing", label: "Grid spacing", min: 0.8, max: 2, step: 0.1, value: P.spacing, unit: " Å", why: "Finer grids see more detail and cost much more (slide 23)." });
      add({ key: "probe", label: "Probe radius", min: 0, max: 2, step: 0.1, value: P.probe, unit: " Å", why: "Added to every atom radius, so gaps narrower than a small molecule count as protein." });
      add({ key: "minBuried", label: "Buried directions needed", min: 1, max: 7, step: 1, value: P.minBuried, format: (v) => `${v} of 7`, why: "Lines checked: the 3 axes and 4 cube diagonals. Protein on both sides counts as one." });
      add({ key: "maxScan", label: "Look for protein within", min: 4, max: 20, step: 1, value: P.maxScan, unit: " Å", why: "How far along each line to look. Very far, and the whole outside of a concave protein looks buried." });
      add({ key: "minVolume", label: "Smallest pocket kept", min: 10, max: 200, step: 10, value: P.minVolume, unit: " Å³" });
    } else {
      add({ key: "rMin", label: "Smallest alpha sphere", min: 2.8, max: 5, step: 0.1, value: P.rMin, unit: " Å", why: "Smaller spheres sit between packed atoms inside the protein (slide 26). fpocket: 3.4 Å." });
      add({ key: "rMax", label: "Largest alpha sphere", min: 4, max: 8, step: 0.1, value: P.rMax, unit: " Å", why: "Larger spheres lie outside the protein or over flat surface. fpocket: 6.2 Å." });
      add({ key: "linkDist", label: "Cluster spheres closer than", min: 1, max: 4, step: 0.1, value: P.linkDist, unit: " Å", why: "Sphere centres this close end up in the same pocket (single linkage)." });
      add({ key: "minSpheres", label: "Smallest pocket kept", min: 3, max: 60, step: 1, value: P.minSpheres, format: (v) => `${v} spheres` });
    }
    const scoreBox = $("#algo-score");
    scoreBox.innerHTML = "";
    const scores = S.method === "grid"
      ? [{ value: "size", label: "Volume" }, { value: "buriedness", label: "Volume × burial" }]
      : [{ value: "size", label: "Size" }, { value: "apolar", label: "Hydrophobicity", title: "spheres touching ≥3 carbon atoms" }, { value: "buriedness", label: "Burial" }];
    scoreBox.appendChild(segmented(scores, P.score, (v) => { P.score = v; run(); }));
    scoreBox.appendChild(el(`<p class="small muted">${S.method === "grid" ? "Burial = how many of the 7 directions are blocked, on average." :
      "Hydrophobicity = number of apolar spheres (touching at least 3 carbon atoms). Burial = protein atoms within 8 Å of the pocket centre."}</p>`));
    const show = $("#algo-show");
    show.innerHTML = "";
    show.appendChild(checkbox("True ligand (green)", S.showLigand, (v) => { S.showLigand = v; draw(); }));
    if (S.method === "grid") show.appendChild(checkbox("All buried grid points, by burial (slow)", S.showRaw, (v) => { S.showRaw = v; run(); }));
    else show.appendChild(checkbox("Every alpha sphere, by size class (slide 26)", S.showRaw, (v) => { S.showRaw = v; draw(); }));
    show.appendChild(checkbox("Spheres at full size (top pocket)", S.fullSpheres, (v) => { S.fullSpheres = v; draw(); }));
  }

  // ------------------------------------------------------------------ run + draw
  function run() {
    const t0 = performance.now();
    const p = { ...S.params, showBuried: S.method === "grid" && S.showRaw ? 3 : 0 };
    S.result = S.method === "grid" ? Pockets.gridPockets(S.atoms, p) : Pockets.alphaPockets(S.atoms, S.prot.alpha, p);
    S.ms = performance.now() - t0;
    S.evals = S.result.pockets.map((pk) => Pockets.evaluate(pk, S.ligXyz));
    const k = S.evals.findIndex((e) => e.dcc <= Pockets.DCC_SUCCESS_A);
    S.trueRank = k >= 0 ? k + 1 : 0;
    record({ protein: S.protId, method: S.method, params: { ...S.params }, trueRank: S.trueRank });
    draw();
    renderRight();
  }

  function draw() {
    const v = S.viewer;
    v.removeAllShapes(); v.removeAllLabels();
    (S.dotModels || []).forEach((m) => v.removeModel(m));
    S.dotModels = [];
    const dots = (pts, col, r, op) => { if (pts.length) S.dotModels.push(PW.dots(v, pts, col, r, op)); };
    S.model.setStyle({}, { cartoon: { color: "#b8c4dc", opacity: 0.72 } });
    S.model.setStyle({ hetflag: true }, { sphere: { color: "#c0c7d4", scale: 0.8 } });
    if (S.ligModel) { v.removeModel(S.ligModel); S.ligModel = null; }
    if (S.showLigand) {
      S.ligModel = v.addModel(S.prot.ligand.molblock, "sdf");
      S.ligModel.setStyle({}, { stick: { radius: 0.24, colorscheme: "greenCarbon" } });
    }
    if (S.showRaw && S.method === "grid" && S.result.buriedPoints) {
      const col = { 3: "#2f4675", 4: "#3b5f9e", 5: "#5585cc", 6: "#8bb6f2", 7: "#d6e8ff" };
      for (let b = 3; b <= 7; b++) dots(S.result.buriedPoints.filter((q) => q[3] === b), col[b], 0.16, 0.9);
    }
    if (S.showRaw && S.method === "alpha") {
      const a = S.prot.alpha, P = S.params, small = [], mid = [], large = [];
      for (let i = 0; i < a.r.length; i++) (a.r[i] < P.rMin ? small : a.r[i] > P.rMax ? large : mid).push([a.x[i], a.y[i], a.z[i]]);
      dots(small, "#4d7cff", 0.14, 0.9); dots(mid, "#ffd166", 0.14, 0.9); dots(large, "#ff6b6b", 0.14, 0.9);
    }
    S.result.pockets.slice(0, RANK_COLORS.length).forEach((pk, i) => {
      const col = RANK_COLORS[i];
      dots(S.method === "grid" ? pk.points : pk.spheres.map((sp) => sp.c), col, S.method === "grid" ? 0.3 : 0.38);
      v.addLabel(`#${i + 1}`, { position: vec(pk.centre), fontSize: 13, fontColor: "#06131c", backgroundColor: col, backgroundOpacity: 0.95, inFront: true, borderThickness: 0 });
    });
    if (S.fullSpheres && S.result.pockets[0]) {
      const pk = S.result.pockets[0];
      const big = v.addShape({ color: RANK_COLORS[0], alpha: S.method === "alpha" ? 0.12 : 0.3 });
      if (S.method === "alpha") pk.spheres.forEach((sp) => big.addSphere({ center: vec(sp.c), radius: sp.r, color: RANK_COLORS[0] }));
      else pk.points.forEach((q) => big.addSphere({ center: vec(q), radius: S.params.spacing / 2, color: RANK_COLORS[0] }));
    }
    v.render();
  }

  // ------------------------------------------------------------------ right panel
  function renderRight() {
    const right = $("#right");
    const n = S.result.pockets.length;
    const verdict = S.trueRank === 1 ? `<span class="yes">Found: the true pocket is #1.</span>`
      : S.trueRank > 0 ? `<span style="color:var(--warn)">The true pocket is #${S.trueRank}.</span>`
      : `<span class="no">The true pocket was not found.</span>`;
    const rows = S.result.pockets.slice(0, 8).map((pk, i) => {
      const e = S.evals[i], hit = e.dcc <= Pockets.DCC_SUCCESS_A;
      const size = S.method === "grid" ? `${pk.volume.toFixed(0)} Å³` : `${pk.n}`;
      return `<tr class="${hit ? "hit" : ""}"><td><span class="swatch" style="background:${RANK_COLORS[i]}"></span>#${i + 1}</td><td class="num">${size}</td>
        <td class="num">${pk.score.toFixed(S.method === "grid" ? 0 : 1)}</td><td class="num">${e.dcc.toFixed(1)}</td><td>${hit ? "✓" : ""}</td></tr>`;
    }).join("");
    right.innerHTML = `<h3>Result</h3>
      <div class="box ${S.trueRank === 1 ? "gold" : ""}"><b>${n} pocket${n === 1 ? "" : "s"}</b> in ${S.ms.toFixed(0)} ms. ${verdict}
        <div class="small muted">A pocket counts as correct when its centre is within ${Pockets.DCC_SUCCESS_A} Å of the ligand's centre (DCC, the usual benchmark criterion).</div></div>
      <table><thead><tr><th>Rank</th><th class="num">${S.method === "grid" ? "Volume" : "Spheres"}</th><th class="num">Score</th><th class="num">DCC Å</th><th></th></tr></thead><tbody>${rows || `<tr><td colspan="5" class="muted">No pockets. Loosen a parameter.</td></tr>`}</tbody></table>
      <div class="row" style="margin-top:8px"><button id="algo-bench">Benchmark on all 10 proteins</button></div>
      <div id="algo-bench-out"></div>
      <h3>Challenges</h3><div id="algo-tasks"></div>`;
    $("#algo-bench").addEventListener("click", benchmark);
    if (S.bench) renderBench();
    renderTasks();
  }

  async function benchmark() {
    const out = $("#algo-bench-out");
    out.innerHTML = `<p class="small muted">Running…</p>`;
    await new Promise((r) => setTimeout(r, 30));
    const rows = [];
    for (const L of S.index.levels) {
      const d = await PW.get("prot", L.id);
      const atoms = Pockets.atomsFromPdb(d.receptor_pdb);
      const lig = d.ligand.atoms.map((a) => a.xyz);
      const res = S.method === "grid" ? Pockets.gridPockets(atoms, S.params) : Pockets.alphaPockets(atoms, d.alpha, S.params);
      const k = res.pockets.findIndex((pk) => Pockets.evaluate(pk, lig).dcc <= Pockets.DCC_SUCCESS_A);
      rows.push({ title: L.title, rank: k >= 0 ? k + 1 : 0 });
      await new Promise((r) => setTimeout(r, 0));
    }
    const top1 = rows.filter((r) => r.rank === 1).length, top3 = rows.filter((r) => r.rank >= 1 && r.rank <= 3).length;
    S.bench = { rows, top1, top3, method: S.method };
    record({ protein: S.protId, method: S.method, params: { ...S.params }, trueRank: S.trueRank, benchmarkTop1: top1 });
    renderBench();
    renderTasks();
  }
  function renderBench() {
    const b = S.bench;
    $("#algo-bench-out").innerHTML = `<div class="box"><b>${b.top1}/10</b> at rank 1, <b>${b.top3}/10</b> in the top 3 (${b.method === "grid" ? "grid" : "alpha spheres"}, these settings).
      <table style="margin-top:6px">${b.rows.map((r) => `<tr><td>${esc(r.title)}</td><td class="num ${r.rank === 1 ? "yes" : r.rank ? "" : "no"}">${r.rank ? `#${r.rank}` : "–"}</td></tr>`).join("")}</table></div>`;
  }

  function record(r) {
    S.last = r;
    TASKS.forEach((t) => {
      if (!PW.progress.tasks[t.id] && t.check(r)) {
        PW.progress.tasks[t.id] = true;
        PW.save();
        toast(`Challenge done: ${t.text.replace(/<[^>]+>/g, "").slice(0, 70)}…`, "good");
      }
    });
  }
  function renderTasks() {
    const box = $("#algo-tasks");
    if (!box) return;
    box.innerHTML = TASKS.map((t) => {
      const done = !!PW.progress.tasks[t.id];
      return `<div class="task ${done ? "done" : ""}"><span class="tick"></span><div><p>${t.text}</p>${done ? `<div class="ans">${t.answer}</div>` : ""}</div></div>`;
    }).join("");
  }

  // ------------------------------------------------------------------ open
  async function openProtein(id) {
    S.protId = id;
    PW.setArg(id);
    const L = S.index.levels.find((x) => x.id === id);
    PW.setHud("Pocket algorithms", `${L.title} · PDB ${L.pdb_id}`);
    S.prot = await PW.get("prot", id);
    S.viewer = PW.resetView();
    S.atoms = Pockets.atomsFromPdb(S.prot.receptor_pdb);
    S.ligXyz = S.prot.ligand.atoms.map((a) => a.xyz);
    S.model = S.viewer.addModel(S.prot.receptor_pdb, "pdb");
    S.ligModel = null;
    S.viewer.zoomTo();
    S.bench = null;
    run();
    $("#stage-overlay").innerHTML = "";
    $("#stage-overlay").appendChild(el(`<div class="legend-card">Pockets are coloured by rank: <span style="color:${RANK_COLORS[0]}">#1</span>
      <span style="color:${RANK_COLORS[1]}">#2</span> <span style="color:${RANK_COLORS[2]}">#3</span> … Tick <b>True ligand</b> to check.</div>`));
  }

  PW.chapter({
    id: "algo", title: "Pocket algorithms", color: "#b69cff", num: "Chapter 2",
    blurb: "Run the grid method and alpha spheres from the lecture on the same proteins. Turn the knobs and watch the ranking change.",
    learn: "LIGSITE grid · fpocket alpha spheres · descriptors · benchmarking",
    progress() {
      const done = TASKS.filter((t) => PW.progress.tasks[t.id]).length;
      return { text: `${done}/${TASKS.length}`, done: done === TASKS.length };
    },
    async enter(arg) {
      const index = await PW.get("index", "index");
      S = { index, method: "grid", params: defaultsFor("grid"), showLigand: false, showRaw: false, fullSpheres: false };
      S.protId = index.levels.some((L) => L.id === arg) ? arg : index.levels[0].id;
      renderLeft();
      await openProtein(S.protId);
    },
    leave() { S = null; },
    recentre() { if (S) { S.viewer.zoomTo(); S.viewer.render(); } },
  });
})();
