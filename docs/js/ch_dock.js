// Chapter 3: dock it yourself. Ported from the docking_game prototype ("Pose the Ligand"):
// move and rotate the crystal ligand back into its pocket and watch the RMSD.
// All pose math is in js/metrics.js (tested in tests/test_metrics.js).
(function () {
  "use strict";
  const { $, $$, esc, el, stars, citation, toast, segmented, checkbox } = PW;
  const M = window.Metrics;

  const GAUGE_MAX_A = 10; // RMSD at which the gauge is full (feel only)

  const MODES = {
    easy: { show: { ghost: true, contacts: false, pocket: true, surface: false, live: true }, stars: 1,
      hint: "The green ghost is the crystal pose. Put the pink ligand on top of it and watch the RMSD drop." },
    medium: { show: { ghost: false, contacts: true, pocket: true, surface: false, live: true }, stars: 2,
      hint: "No ghost. Dashed lines show the crystal interactions; they turn green when you are within contact distance." },
    hard: { show: { ghost: false, contacts: false, pocket: false, surface: true, live: false }, stars: 3,
      hint: "Only the protein surface. Find the pocket, place the ligand, then press Submit to see how close you were." },
  };

  let S = null;

  // ------------------------------------------------------------------ geometry helpers (from the prototype)
  function molblockAt(xyz) {
    const lines = S.molLines.slice();
    xyz.forEach((p, i) => { lines[4 + i] = p.map((v) => v.toFixed(4).padStart(10)).join("") + lines[4 + i].slice(30); });
    return lines.join("\n");
  }
  // ------------------------------------------------------------------ drawing
  function applyStaticStyles() {
    const v = S.viewer, show = S.show, rec = S.recModel;
    rec.setStyle({}, { cartoon: { color: "#9aa3b2", opacity: 0.85 } });
    rec.setStyle({ hetflag: true }, { sphere: { color: "#c0c7d4", scale: 0.8 } });
    if (show.pocket) S.pocketSel.forEach((sel) => rec.setStyle(sel, { stick: { radius: 0.13, colorscheme: "grayCarbon" } }, true));
    v.removeAllLabels();
    if (show.contacts) {
      const seen = new Set();
      S.prot.ligand.contacts.forEach((c) => {
        if (c.kind !== "water") rec.setStyle({ chain: c.partner_chain, resi: c.partner_resi }, { stick: { radius: 0.17, colorscheme: "orangeCarbon" } }, true);
        const key = c.partner_resn + c.partner_resi;
        if (seen.has(key)) return;
        seen.add(key);
        v.addLabel(c.kind === "water" ? "water" : `${c.partner_resn} ${c.partner_resi}`, {
          position: { x: c.partner_xyz[0], y: c.partner_xyz[1], z: c.partner_xyz[2] },
          fontSize: 12, fontColor: "white", backgroundColor: "black", backgroundOpacity: 0.55, inFront: true });
      });
      S.waterModel.setStyle({}, { sphere: { radius: 0.35, color: "#ff5c5c" } });
    } else S.waterModel.setStyle({}, {});
    S.truthModel.setStyle({}, (show.ghost || S.revealed) ? { stick: { radius: 0.2, colorscheme: "greenCarbon", opacity: 0.55 } } : {});
    v.removeAllSurfaces();
    if (show.surface) v.addSurface($3Dmol.SurfaceType.MS, { opacity: 0.6, color: "#e8ecf2" }, { model: rec });
  }

  function drawLigand(xyz) {
    const v = S.viewer;
    if (S.ligModel) v.removeModel(S.ligModel);
    S.ligModel = v.addModel(molblockAt(xyz), "sdf");
    S.ligModel.setStyle({}, { stick: { radius: 0.24, colorscheme: "magentaCarbon" } });
    v.removeAllShapes();
    const piv = S.controls && S.controls.pivot; // the atom the ligand is being turned around
    if (piv) v.addSphere({ center: { x: piv[0], y: piv[1], z: piv[2] }, radius: 0.45, color: "#22d3ee", opacity: 0.9 });
    if (S.show.contacts) {
      S.prot.ligand.contacts.forEach((c) => {
        const p = xyz[c.lig_index], ok = M.dist(p, c.partner_xyz) <= S.cfg.polar_contact_cutoff_A;
        v.addCylinder({ start: { x: p[0], y: p[1], z: p[2] }, end: { x: c.partner_xyz[0], y: c.partner_xyz[1], z: c.partner_xyz[2] },
          radius: 0.06, dashed: true, dashLength: 0.25, gapLength: 0.2, color: ok ? "#2ecc71" : "#ffb020", fromCap: 1, toCap: 1 });
      });
    }
  }

  function evaluate(xyz) {
    const L = S.prot.ligand;
    return {
      sym: M.symmetryRmsd(xyz, S.refXyz, L.automorphisms).rmsd,
      naive: M.rmsd(xyz, S.refXyz),
      made: L.contacts.filter((c) => M.dist(xyz[c.lig_index], c.partner_xyz) <= S.cfg.polar_contact_cutoff_A).length,
      centroid: M.dist(M.centroid(xyz), S.refCentroid),
      clashes: M.countClashes(xyz, S.recXyz, S.cfg.clash_cutoff_A),
    };
  }

  function updateReadouts(xyz) {
    const e = evaluate(xyz), visible = S.show.live || S.revealed, ok = S.cfg.success_rmsd_A;
    const q = (id) => document.getElementById(id);
    q("dk-rmsd").textContent = visible ? `${e.sym.toFixed(2)} Å` : "hidden";
    q("dk-naive").textContent = visible ? `${e.naive.toFixed(2)} Å` : "–";
    q("dk-centroid").textContent = visible ? `${e.centroid.toFixed(2)} Å` : "–";
    q("dk-clashes").textContent = e.clashes;
    q("dk-contacts").textContent = S.show.contacts || S.revealed ? `${e.made} / ${S.prot.ligand.contacts.length}` : "hidden";
    q("dk-gauge").style.width = visible ? `${Math.min(e.sym / GAUGE_MAX_A, 1) * 100}%` : "0";
    const b = q("dk-badge");
    if (!visible) { b.textContent = ""; b.className = "chip"; }
    else if (e.sym <= ok) { b.textContent = "success"; b.className = "chip yes"; }
    else if (e.sym <= 2 * ok) { b.textContent = "close"; b.className = "chip"; b.style.color = "var(--warn)"; }
    else { b.textContent = "off"; b.className = "chip no"; }
    const card = q("dk-contacts-card");
    card.style.display = (S.show.contacts || S.revealed) && S.prot.ligand.contacts.length ? "" : "none";
    q("dk-contacts-body").innerHTML = S.prot.ligand.contacts.map((c) => {
      const d = M.dist(xyz[c.lig_index], c.partner_xyz), hit = d <= S.cfg.polar_contact_cutoff_A;
      return `<tr><td>${esc(c.lig_atom)}</td><td>${esc(c.partner)}<br><span class="muted">${esc(c.kind)}</span></td><td class="num">${c.distance.toFixed(2)}</td><td class="num ${hit ? "yes" : "no"}">${d.toFixed(2)}</td></tr>`;
    }).join("");
  }

  function refresh() {
    const xyz = S.curXyz = M.applyPose(S.refXyz, S.pose, S.refCentroid);
    drawLigand(xyz);
    updateReadouts(xyz);
    S.viewer.render();
  }

  function newStart() {
    S.revealed = false;
    S.pose = { R: M.randomRotation(), t: M.scale(S.outward, S.cfg.start_offset_A) };
    applyStaticStyles();
    refresh();
  }
  function recentre() {
    // Fit the receptor atoms around the pocket. Picked by index: 3Dmol's "within" does not reach across
    // models, and with it the view fitted only the 9-atom ghost and ended up inside the protein.
    const r = S.cfg.start_offset_A + 9, index = [];
    S.recXyz.forEach((q, i) => { if (M.dist(q, S.refCentroid) <= r) index.push(i); });
    S.viewer.zoomTo({ model: S.recModel, index });
    S.viewer.render();
  }

  // ------------------------------------------------------------------ input (js/controls.js)
  function bindInput() {
    S.controls = window.LigandControls.attach({
      viewer: PW.view(), el: $("#stage"), M, center: S.refCentroid || [0, 0, 0],
      getPose: () => S.pose, getXyz: () => S.curXyz, refresh, enabled: () => !!S.pose,
    });
  }
  function unbindInput() {
    if (S && S.controls) { S.controls.detach(); S.controls = null; }
  }

  // ------------------------------------------------------------------ panels
  function renderLeft() {
    const left = $("#left");
    left.innerHTML = `<h3>Dock it yourself</h3>
      <p class="small">Docking = <b>search</b> + <b>scoring</b> (slide 32). Here <i>you</i> are the search. Put the ligand back where the crystal structure has it.</p>
      <div class="levels" id="dk-cases"></div>
      <h3>Difficulty</h3><div id="dk-mode"></div><p class="small muted" id="dk-mode-hint"></p>
      <div class="row"><button id="dk-start">New start</button></div>
      <h3>Show</h3><div id="dk-show"></div>
      <details class="small"><summary><b>Colour legend</b></summary>
        <p style="margin:4px 0"><b>Atoms</b> (sticks): <span style="color:#6d8df0">N</span> nitrogen and <span style="color:#ee7b7b">O</span> oxygen, the hydrogen-bond donors and acceptors; <span style="color:#f0cd62">S</span> sulfur. Hydrogens are not shown: X-ray structures at this resolution do not resolve them.</p>
        <p style="margin:4px 0"><b>Carbon colour = which molecule</b>: <span style="color:#ff8fd0">pink</span> your ligand, <span style="color:#7ee787">green</span> the crystal pose, grey the pocket residues, <span style="color:#ffb454">orange</span> residues making the crystal contacts. <span style="color:#ff5c5c">Red balls</span>: waters in the contacts. Grey spheres: metal ions.</p>
        <p style="margin:4px 0"><b>Dashed lines</b>: a crystal contact, <span style="color:#2ecc71">green</span> when your atom is within ${S.cfg.polar_contact_cutoff_A} Å of its partner, <span style="color:#ffb020">amber</span> when not yet.</p></details>
      <h3>Controls</h3>
      <ul class="small ctl-help">
        <li><b>Camera:</b> drag empty space to rotate, scroll to zoom, middle-drag to pan.</li>
        <li><b>Move the ligand:</b> grab it and drag (the cursor turns into a hand). It moves in the screen plane; <kbd>Shift</kbd>+scroll pushes or pulls it. <kbd>Ctrl</kbd>/<kbd>⌘</kbd>+drag moves it from anywhere.</li>
        <li><b>Rotate the ligand:</b> <kbd>Shift</kbd>+drag, or right-drag on it. Start on an atom and the ligand turns around that atom (cyan ball), so an end you have already placed stays put; from empty space it turns around its centre. <kbd>Alt</kbd>+scroll spins it.</li>
        <li><b>Keys:</b> <kbd>W</kbd><kbd>A</kbd><kbd>S</kbd><kbd>D</kbd> move, <kbd>R</kbd>/<kbd>F</kbd> toward/away, arrows rotate, <kbd>Q</kbd>/<kbd>E</kbd> spin; <kbd>Shift</kbd> for small steps.</li></ul>`;
    const list = $("#dk-cases");
    S.cases.forEach((c, i) => {
      const L = S.index.levels.find((x) => x.id === c.level);
      const s = PW.progress.docking[c.level];
      const b = el(`<button class="lvl ${i === S.caseIdx ? "on" : ""}"><span class="n">${i + 1}</span>
        <span class="t">${esc(L.ligand_label)}<small><span class="diff ${c.difficulty}">${c.difficulty}</span> ${esc(L.title.split(" + ")[0])} · ${L.n_ligand_atoms} atoms</small></span>
        ${s ? stars(s) : `<span class="stars"><span class="off">★★★</span></span>`}</button>`);
      b.addEventListener("click", () => openCase(i));
      list.appendChild(b);
    });
    $("#dk-mode").appendChild(segmented([{ value: "easy", label: "Easy" }, { value: "medium", label: "Medium" }, { value: "hard", label: "Hard" }], S.mode, setMode));
    $("#dk-start").addEventListener("click", newStart);
    renderShow();
  }
  function renderShow() {
    const box = $("#dk-show");
    box.innerHTML = "";
    [["ghost", "Crystal pose (ghost)"], ["contacts", "Key interactions"], ["pocket", "Pocket residues"], ["surface", "Protein surface"], ["live", "Live RMSD"]].forEach(([k, label]) => {
      box.appendChild(checkbox(label, S.show[k], (v) => { S.show[k] = v; applyStaticStyles(); refresh(); }));
    });
    $("#dk-mode-hint").textContent = MODES[S.mode].hint;
  }
  function setMode(m) {
    S.mode = m;
    S.show = { ...MODES[m].show };
    renderShow();
    newStart();
  }

  function renderRight() {
    const c = S.cases[S.caseIdx], L = S.index.levels.find((x) => x.id === c.level);
    $("#right").innerHTML = `<h3>Ligand</h3><div class="big" style="font-size:22px">${esc(L.ligand_label)}</div>
      <p class="small muted" style="margin:2px 0">${esc(L.ligand_name.toLowerCase())} (${L.ligand_id})</p>
      <p class="small muted">${esc(L.title)} · PDB <a href="https://www.rcsb.org/structure/${L.pdb_id}" target="_blank" rel="noopener">${L.pdb_id}</a> · ${esc(c.note)}</p>
      <div class="box acc">
        <div class="row" style="justify-content:space-between"><span class="small muted">Symmetry-corrected RMSD</span><span id="dk-badge" class="chip"></span></div>
        <div class="big" id="dk-rmsd">–</div>
        <div class="gauge"><i id="dk-gauge"></i></div>
        <div class="small muted">0 · success ≤ ${S.cfg.success_rmsd_A} Å · ${GAUGE_MAX_A} Å</div>
        <div class="readout" style="margin-top:8px">
          <span title="Atom i matched to atom i, ignoring symmetry">Naive RMSD</span><b id="dk-naive">–</b>
          <span>Centroid distance</span><b id="dk-centroid">–</b>
          <span title="Ligand-protein heavy atoms closer than ${S.cfg.clash_cutoff_A} Å">Clashes</span><b id="dk-clashes">–</b>
          <span>Polar contacts made</span><b id="dk-contacts">–</b></div>
        <div class="row" style="margin-top:10px"><button class="primary" id="dk-submit">Submit pose</button></div>
      </div>
      <div class="box" id="dk-contacts-card"><b class="small">Interactions in the crystal</b>
        <table><thead><tr><th>Ligand</th><th>Partner</th><th class="num">Crystal</th><th class="num">Now</th></tr></thead><tbody id="dk-contacts-body"></tbody></table></div>
      <h3>Your attempts</h3><table><thead><tr><th>#</th><th>Mode</th><th class="num">RMSD</th><th class="num">Naive</th></tr></thead><tbody id="dk-history"></tbody></table>
      <details style="margin-top:12px"><summary class="small">What is the RMSD here?</summary>
        <p class="small">The RMSD to the crystal pose is a <b>metric</b> for judging a pose, not a score that predicts binding. In plain words: the average distance between each ligand atom and the same atom in the crystal structure.</p>
        <p class="small">RMSD = √(Σ|x<sub>i</sub> − x<sub>i</sub><sup>crystal</sup>|² / N) over the ${L.n_ligand_atoms} heavy atoms. Nothing is superimposed: both poses stay in the protein's frame, so it measures both where the ligand sits and how it is turned.
        Symmetric molecules have several valid atom matchings (${S.prot.ligand.automorphisms.length} here); the headline number takes the best one, the naive one does not.</p>
        <p class="small">Under ${S.cfg.success_rmsd_A} Å usually counts as a docking success (slide 35). RMSD alone is not enough, though: it does not check clashes, bond geometry or the key interactions, and many poses under ${S.cfg.success_rmsd_A} Å are physically implausible (<a href="https://doi.org/10.1039/D3SC04185A" target="_blank" rel="noopener">PoseBusters: Buttenschoen, Morris &amp; Deane, <i>Chem. Sci.</i> 2024</a>). That is why clashes and contacts are counted too.</p></details>`;
    $("#dk-submit").addEventListener("click", submit);
    renderHistory();
  }
  function renderHistory() {
    const h = S.history.filter((x) => x.case === S.cases[S.caseIdx].level);
    $("#dk-history").innerHTML = h.length ? h.map((x, i) => `<tr><td>${i + 1}</td><td>${x.mode}</td><td class="num ${x.sym <= S.cfg.success_rmsd_A ? "yes" : ""}">${x.sym.toFixed(2)}</td><td class="num">${x.naive.toFixed(2)}</td></tr>`).reverse().join("")
      : `<tr><td colspan="4" class="muted">Nothing submitted yet.</td></tr>`;
  }
  // Stars follow the aids that were switched on, not the mode's name: ticking the ghost in Hard mode
  // is Easy mode.
  function earnedStars() {
    if (S.show.ghost) return MODES.easy.stars;
    if (S.show.live || S.show.contacts || S.show.pocket) return MODES.medium.stars;
    return MODES.hard.stars;
  }
  function submit() {
    const e = evaluate(M.applyPose(S.refXyz, S.pose, S.refCentroid));
    const c = S.cases[S.caseIdx];
    S.history.push({ case: c.level, mode: S.mode, sym: e.sym, naive: e.naive });
    if (e.sym <= S.cfg.success_rmsd_A && S.revealed) {
      toast(`${e.sym.toFixed(2)} Å, but the crystal pose was already shown: no stars. Press <b>New start</b> to try again.`, "warn");
    } else if (e.sym <= S.cfg.success_rmsd_A) {
      const n = earnedStars();
      PW.progress.docking[c.level] = Math.max(n, PW.progress.docking[c.level] || 0);
      PW.save();
      toast(`Docked! ${e.sym.toFixed(2)} Å. ${stars(n)}${n < 3 ? " Fewer aids switched on (Show) earn more stars." : ""}`, "good");
      renderLeft();
    } else toast(`${e.sym.toFixed(2)} Å: not yet under ${S.cfg.success_rmsd_A} Å. The crystal pose is now shown in green.`, "warn");
    S.revealed = true;
    applyStaticStyles(); refresh(); renderHistory();
  }

  // ------------------------------------------------------------------ open
  async function openCase(i) {
    S.caseIdx = i;
    const c = S.cases[i], L = S.index.levels.find((x) => x.id === c.level);
    PW.setArg(i + 1);
    PW.setHud("Dock it yourself", `${L.ligand_label} into ${L.title.split(" + ")[0]}`);
    S.prot = await PW.get("prot", c.level);
    const v = S.viewer = PW.resetView();
    S.refXyz = S.prot.ligand.atoms.map((a) => a.xyz);
    S.refCentroid = M.centroid(S.refXyz);
    unbindInput(); bindInput(); // the controls rotate about this case's ligand centre
    S.recXyz = S.prot.receptor_pdb.split("\n").filter((l) => l.startsWith("ATOM") || l.startsWith("HETATM")).map((l) => [+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)]);
    S.molLines = S.prot.ligand.molblock.split("\n");
    S.outward = M.normalize(M.sub(S.refCentroid, M.centroid(S.recXyz)));
    S.pocketSel = S.prot.pocket_residues.map((r) => ({ chain: r.chain, resi: r.resi }));
    S.recModel = v.addModel(S.prot.receptor_pdb, "pdb");
    S.waterModel = v.addModel(S.prot.water_pdb, "pdb");
    S.truthModel = v.addModel(S.prot.ligand.molblock, "sdf");
    S.ligModel = null;
    renderLeft();
    renderRight();
    setMode(S.mode);
    recentre();
    const hint = $("#stage-hint");
    hint.innerHTML = `<span id="dk-hint"><b>Drag ligand</b> move · <kbd>Shift</kbd>+drag rotate (start on an atom to pivot there) · <b>drag background</b> camera · <kbd>Shift</kbd>+scroll depth · <kbd>W</kbd><kbd>A</kbd><kbd>S</kbd><kbd>D</kbd>/arrows</span>`;
    hint.classList.remove("hidden");
  }

  PW.chapter({
    id: "dock", title: "Dock it yourself", color: "#ff8fd0", num: "Chapter 3",
    blurb: "You are the search algorithm. Push and turn each ligand back into its pocket until the RMSD drops under 2 Å.",
    learn: "Docking = search + scoring · RMSD · symmetry · clashes and contacts",
    progress(index) {
      const done = index.docking.cases.filter((c) => PW.progress.docking[c.level]).length;
      return { text: `${done}/${index.docking.cases.length}`, done: done === index.docking.cases.length };
    },
    async enter(arg) {
      const index = await PW.get("index", "index");
      S = { index, cfg: index.docking, cases: index.docking.cases, mode: "easy", history: [], show: { ...MODES.easy.show } };
      const i = Math.max(0, Math.min(S.cases.length - 1, (parseInt(arg, 10) || 1) - 1));
      await openCase(i);
    },
    tour: () => [
      { sel: "#dk-cases", title: "Five ligands", text: "From a rigid fragment to imatinib. Each goes back into its own crystal pocket." },
      { sel: "#dk-mode", title: "Difficulty", text: "Easy shows the crystal pose as a green ghost, medium the key contacts, hard only the protein surface. Stars depend on the aids that are switched on." },
      { sel: "#stage", title: "Move and turn the ligand", text: "<b>Drag the ligand</b> to move it, <b>Shift+drag</b> (or right-drag on it) to rotate it, drag empty space to turn the camera. Start a rotation on an atom to keep that atom in place." },
      { sel: "#dk-hint", title: "Controls at a glance", text: "The short version stays here; the full list (keys too) is in the left panel." },
      { sel: "#dk-submit", title: "Submit", text: "Checks your pose against the crystal structure. After a miss the crystal pose is shown; <b>New start</b> gives a fresh try." },
    ],
    leave() { unbindInput(); S = null; },
    recentre() { if (S) recentre(); },
  });
})();
