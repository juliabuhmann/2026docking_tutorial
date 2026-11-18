// Chapter 3: dock it yourself. Ported from the docking_game prototype ("Pose the Ligand"):
// move and rotate the crystal ligand back into its pocket and watch the RMSD.
// All pose math is in js/metrics.js (tested in tests/test_metrics.js).
(function () {
  "use strict";
  const { $, $$, esc, el, stars, citation, toast, segmented, checkbox } = PW;
  const M = window.Metrics;

  // UI feel only, no scientific meaning.
  const KEY_STEP_A = 0.25, KEY_STEP_FINE_A = 0.05, KEY_ROT_DEG = 5, KEY_ROT_FINE_DEG = 1;
  const DRAG_ROT_RAD_PER_PX = 0.01, WHEEL_MOVE_A_PER_UNIT = 0.01, WHEEL_ROT_RAD_PER_UNIT = 0.003;
  const GAUGE_MAX_A = 10;

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
  function screenBasis() {
    const c = M.add(S.refCentroid, S.pose.t), v = S.viewer;
    const toP = (p) => ({ x: p[0], y: p[1], z: p[2] });
    const s0 = v.modelToScreen(toP(c));
    const cols = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map((e) => v.modelToScreen(toP(M.add(c, e))));
    const Jx = cols.map((s) => s.x - s0.x), Jy = cols.map((s) => s.y - s0.y);
    const a = M.dot(Jx, Jx), b = M.dot(Jx, Jy), d = M.dot(Jy, Jy), det = a * d - b * b;
    const pinv = (sx, sy) => M.add(M.scale(Jx, (d * sx - b * sy) / det), M.scale(Jy, (-b * sx + a * sy) / det));
    const right = pinv(1, 0), down = pinv(0, 1), up = M.scale(down, -1);
    return { right, down, up, toward: M.normalize(M.cross(right, up)) };
  }
  const translateScreen = (dx, dy) => { const B = screenBasis(); S.pose.t = M.add(S.pose.t, M.add(M.scale(B.right, dx), M.scale(B.down, dy))); };
  const translateDepth = (a) => { S.pose.t = M.add(S.pose.t, M.scale(screenBasis().toward, a)); };
  function rotate(axis, angle) {
    if (M.norm(axis) === 0 || angle === 0) return;
    S.pose.R = M.orthonormalize(M.matMul(M.axisAngle(axis, angle), S.pose.R));
  }
  function rotateScreen(dx, dy) {
    const B = screenBasis();
    rotate(M.add(M.scale(M.normalize(B.up), dx), M.scale(M.normalize(B.right), dy)), Math.hypot(dx, dy) * DRAG_ROT_RAD_PER_PX);
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
    q("dk-contacts").textContent = `${e.made} / ${S.prot.ligand.contacts.length}`;
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
    const xyz = M.applyPose(S.refXyz, S.pose, S.refCentroid);
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
    S.viewer.zoomTo({ within: { distance: S.cfg.start_offset_A + 9, sel: { model: S.truthModel } } });
    S.viewer.render();
  }

  // ------------------------------------------------------------------ input
  function bindInput() {
    const stage = $("#stage");
    let drag = null;
    const actionFor = (ev) => (S.tool === "move" || ev.shiftKey ? "move" : S.tool === "rotate" || ev.altKey ? "rotate" : null);
    const pt = (ev) => (ev.touches ? ev.touches[0] : ev);
    const h = S.handlers = {
      down(ev) {
        if (ev.touches && ev.touches.length > 1) return;
        const a = actionFor(ev);
        if (!a) return;
        ev.stopPropagation(); ev.preventDefault();
        drag = { a, x: pt(ev).clientX, y: pt(ev).clientY };
      },
      move(ev) {
        if (!drag) return;
        ev.stopPropagation(); ev.preventDefault();
        const p = pt(ev), dx = p.clientX - drag.x, dy = p.clientY - drag.y;
        drag.x = p.clientX; drag.y = p.clientY;
        if (drag.a === "move") translateScreen(dx, dy); else rotateScreen(dx, dy);
        refresh();
      },
      up() { drag = null; },
      wheel(ev) {
        const a = actionFor(ev);
        if (!a) return;
        ev.stopPropagation(); ev.preventDefault();
        const delta = ev.deltaY || ev.deltaX;
        if (a === "move") translateDepth(-delta * WHEEL_MOVE_A_PER_UNIT); else rotate(screenBasis().toward, delta * WHEEL_ROT_RAD_PER_UNIT);
        refresh();
      },
      key(ev) {
        if (ev.target.closest && ev.target.closest("input, textarea, select")) return;
        const step = ev.shiftKey ? KEY_STEP_FINE_A : KEY_STEP_A;
        const rot = ((ev.shiftKey ? KEY_ROT_FINE_DEG : KEY_ROT_DEG) * Math.PI) / 180;
        const B = screenBasis(), px = step / Math.max(M.norm(B.right), 1e-9);
        const moves = {
          w: () => translateScreen(0, -px), s: () => translateScreen(0, px), a: () => translateScreen(-px, 0), d: () => translateScreen(px, 0),
          r: () => translateDepth(step), f: () => translateDepth(-step),
          arrowleft: () => rotate(B.up, -rot), arrowright: () => rotate(B.up, rot), arrowup: () => rotate(B.right, -rot), arrowdown: () => rotate(B.right, rot),
          q: () => rotate(B.toward, rot), e: () => rotate(B.toward, -rot),
        };
        const fn = moves[ev.key.toLowerCase()];
        if (!fn) return;
        ev.preventDefault(); fn(); refresh();
      },
    };
    stage.addEventListener("mousedown", h.down, true);
    stage.addEventListener("touchstart", h.down, { capture: true, passive: false });
    window.addEventListener("mousemove", h.move, true);
    window.addEventListener("touchmove", h.move, { capture: true, passive: false });
    window.addEventListener("mouseup", h.up, true);
    window.addEventListener("touchend", h.up, true);
    stage.addEventListener("wheel", h.wheel, { capture: true, passive: false });
    window.addEventListener("keydown", h.key);
  }
  function unbindInput() {
    if (!S || !S.handlers) return;
    const h = S.handlers, stage = $("#stage");
    stage.removeEventListener("mousedown", h.down, true);
    stage.removeEventListener("touchstart", h.down, true);
    window.removeEventListener("mousemove", h.move, true);
    window.removeEventListener("touchmove", h.move, true);
    window.removeEventListener("mouseup", h.up, true);
    window.removeEventListener("touchend", h.up, true);
    stage.removeEventListener("wheel", h.wheel, true);
    window.removeEventListener("keydown", h.key);
    S.handlers = null;
  }

  // ------------------------------------------------------------------ panels
  function renderLeft() {
    const left = $("#left");
    left.innerHTML = `<h3>Dock it yourself</h3>
      <p class="small">Docking = <b>search</b> + <b>scoring</b> (slide 32). Here <i>you</i> are the search. Put the ligand back where the crystal structure has it.</p>
      <div class="levels" id="dk-cases"></div>
      <h3>Difficulty</h3><div id="dk-mode"></div><p class="small muted" id="dk-mode-hint"></p>
      <h3>Mouse does</h3><div id="dk-tool"></div>
      <div class="row"><button id="dk-start">New start</button></div>
      <h3>Show</h3><div id="dk-show"></div>
      <h3>Controls</h3>
      <p class="small"><b>Move ligand</b>: pick the tool or hold <kbd>Shift</kbd> and drag; scroll = towards/away.<br>
      <b>Rotate ligand</b>: pick the tool or hold <kbd>Alt</kbd> and drag; scroll = spin.<br>
      Keys: <kbd>W</kbd><kbd>A</kbd><kbd>S</kbd><kbd>D</kbd> move, <kbd>R</kbd>/<kbd>F</kbd> depth, arrows rotate, <kbd>Q</kbd>/<kbd>E</kbd> spin, <kbd>Shift</kbd> = fine.</p>`;
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
    $("#dk-tool").appendChild(segmented([{ value: "camera", label: "Camera" }, { value: "move", label: "Move ligand" }, { value: "rotate", label: "Rotate ligand" }], S.tool, (t) => { S.tool = t; }));
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
        <p class="small">RMSD = √(Σ|x<sub>i</sub> − x<sub>i</sub><sup>crystal</sup>|² / N) over the ${L.n_ligand_atoms} heavy atoms. Nothing is superimposed: both poses stay in the protein's frame.
        Symmetric molecules have several valid atom matchings (${S.prot.ligand.automorphisms.length} here); the headline number takes the best one, the naive one does not.
        Under ${S.cfg.success_rmsd_A} Å usually counts as a docking success (slide 35).</p></details>`;
    $("#dk-submit").addEventListener("click", submit);
    renderHistory();
  }
  function renderHistory() {
    const h = S.history.filter((x) => x.case === S.cases[S.caseIdx].level);
    $("#dk-history").innerHTML = h.length ? h.map((x, i) => `<tr><td>${i + 1}</td><td>${x.mode}</td><td class="num ${x.sym <= S.cfg.success_rmsd_A ? "yes" : ""}">${x.sym.toFixed(2)}</td><td class="num">${x.naive.toFixed(2)}</td></tr>`).reverse().join("")
      : `<tr><td colspan="4" class="muted">Nothing submitted yet.</td></tr>`;
  }
  function submit() {
    const e = evaluate(M.applyPose(S.refXyz, S.pose, S.refCentroid));
    const c = S.cases[S.caseIdx];
    S.history.push({ case: c.level, mode: S.mode, sym: e.sym, naive: e.naive });
    if (e.sym <= S.cfg.success_rmsd_A) {
      const n = MODES[S.mode].stars;
      PW.progress.docking[c.level] = Math.max(n, PW.progress.docking[c.level] || 0);
      PW.save();
      toast(`Docked! ${e.sym.toFixed(2)} Å on ${S.mode}. ${stars(n)}${n < 3 ? " Try a harder mode for more stars." : ""}`, "good");
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
      S = { index, cfg: index.docking, cases: index.docking.cases, mode: "easy", tool: "camera", history: [], show: { ...MODES.easy.show } };
      const i = Math.max(0, Math.min(S.cases.length - 1, (parseInt(arg, 10) || 1) - 1));
      bindInput();
      await openCase(i);
    },
    leave() { unbindInput(); S = null; },
    recentre() { if (S) recentre(); },
  });
})();
