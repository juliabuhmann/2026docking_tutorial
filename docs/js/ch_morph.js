// Chapter 4: apo to holo. Morph between two crystal structures of the same protein and watch the
// pocket appear. The in-between frames are a straight-line interpolation, not a physical pathway.
(function () {
  "use strict";
  const { $, esc, el, citation, slider, checkbox, segmented, vec } = PW;
  const PLAY_SECONDS = 3.5;           // one sweep start -> end
  const DISP_MAX_A = 10;              // colour scale for atom displacement

  let S = null;

  // displacement colour: blue (still) -> yellow -> red (moves a lot)
  function dispColor(d) {
    const t = Math.min(d / DISP_MAX_A, 1);
    const stops = [[0, [110, 140, 200]], [0.5, [255, 209, 102]], [1, [255, 107, 107]]];
    let i = 0;
    while (i < stops.length - 2 && t > stops[i + 1][0]) i++;
    const [t0, c0] = stops[i], [t1, c1] = stops[i + 1], f = (t - t0) / (t1 - t0);
    const c = c0.map((v, k) => Math.round(v + (c1[k] - v) * f));
    return `rgb(${c[0]},${c[1]},${c[2]})`;
  }

  function coordsAt(t) {
    const s = S.data.start_xyz, e = S.data.end_xyz, out = new Float64Array(s.length);
    for (let i = 0; i < s.length; i++) out[i] = s[i] + (e[i] - s[i]) * t;
    return out;
  }

  function setFrame(t) {
    S.t = t;
    const xyz = coordsAt(t);
    S.atoms3d.forEach((a, i) => { a.x = xyz[3 * i]; a.y = xyz[3 * i + 1]; a.z = xyz[3 * i + 2]; });
    S.xyz = xyz;
    styleProtein();
    updateReadout();
  }

  function styleProtein() {
    const m = S.model;
    if (S.colorBy === "motion") m.setStyle({}, { cartoon: { colorfunc: (a) => S.dispColors[a.index] } });
    else m.setStyle({}, { cartoon: { color: "#9fb3d6" } });
    if (S.showSurface) {
      S.viewer.removeAllSurfaces();
      S.viewer.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.35, color: "#dfe6f2" }, { model: m });
    }
    if (S.ligModel) S.ligModel.setStyle({}, S.showLigand ? { stick: { radius: 0.24, colorscheme: "greenCarbon" } } : {});
    // residues the story is about (e.g. Abl's DFG motif and Tyr393), as labelled sticks
    S.viewer.removeAllLabels();
    (S.data.highlight || []).forEach((h) => {
      m.setStyle({ resi: h.resi }, { stick: { radius: 0.3, colorscheme: "orangeCarbon" } }, true);
      const ca = S.atoms3d.filter((a) => h.resi.includes(a.resi) && a.atom === "CA");
      if (ca.length) {
        const c = ca.reduce((acc, a) => [acc[0] + a.x / ca.length, acc[1] + a.y / ca.length, acc[2] + a.z / ca.length], [0, 0, 0]);
        S.viewer.addLabel(h.label, { position: vec(c), fontSize: 12, fontColor: "#06131c", backgroundColor: "#ffb454", backgroundOpacity: 0.95, inFront: true });
      }
    });
    (S.pocketLabels || []).forEach((l) => S.viewer.addLabel(l.text, l.opts));
    S.viewer.render();
  }

  function clashes() {
    const c2 = S.clashA * S.clashA;
    let n = 0;
    for (const l of S.ligXyz) for (let i = 0; i < S.xyz.length / 3; i++) {
      const dx = S.xyz[3 * i] - l[0], dy = S.xyz[3 * i + 1] - l[1], dz = S.xyz[3 * i + 2] - l[2];
      if (dx * dx + dy * dy + dz * dz < c2) n++;
    }
    return n;
  }

  function updateReadout() {
    const q = (id) => document.getElementById(id);
    if (!q("mo-t")) return;
    q("mo-t").textContent = `${Math.round(S.t * 100)}%`;
    if (q("mo-t2")) q("mo-t2").textContent = `${Math.round(S.t * 100)}%`;
    q("mo-clash").textContent = S.ligXyz.length ? clashes() : "–";
    S.slider.set(Math.round(S.t * 100));
  }

  function findPocketsNow() {
    const atoms = { ...S.baseAtoms, xyz: S.xyz };
    const res = Pockets.gridPockets(atoms, Pockets.GRID_DEFAULTS);
    const ev = res.pockets.map((p) => Pockets.evaluate(p, S.ligXyz));
    const k = ev.findIndex((e) => e.dcc <= Pockets.DCC_SUCCESS_A);
    (S.dotModels || []).forEach((mm) => S.viewer.removeModel(mm));
    const cols = ["#ffd166", "#5eead4", "#b69cff"];
    S.pocketLabels = [];
    S.dotModels = res.pockets.slice(0, 3).map((p, i) => {
      S.pocketLabels.push({ text: `#${i + 1}`, opts: { position: vec(p.centre), fontSize: 13, fontColor: "#06131c", backgroundColor: cols[i], backgroundOpacity: 0.95, inFront: true } });
      return PW.dots(S.viewer, p.points, cols[i], 0.3);
    });
    styleProtein();
    const pct = Math.round(S.t * 100), A = Pockets.DCC_SUCCESS_A;
    const near = ev.length ? Math.min(...ev.map((e) => e.dcc)) : null;
    const msg = k >= 0 ? `At ${pct}%: a pocket at the ligand site, rank <b>#${k + 1}</b>, ${res.pockets[k].volume.toFixed(0)} Å³ (centre ${ev[k].dcc.toFixed(1)} Å from the ligand's).`
      : `At ${pct}%: <b>no pocket</b> centred within ${A} Å of the ligand. ${near === null ? "No pockets at all." : `The nearest pocket centre is ${near.toFixed(1)} Å away${near < 2 * A ? ": close, so look at the dots; part of the site may already be a pocket" : ""}.`}`;
    document.getElementById("mo-pocket").innerHTML = msg;
    if (pct === 0) taskDone("q1");
  }
  function clearPockets() {
    (S.dotModels || []).forEach((mm) => S.viewer.removeModel(mm));
    S.dotModels = []; S.pocketLabels = [];
    const msg = document.getElementById("mo-pocket");
    if (msg) msg.innerHTML = "Grid method, Chapter 2 defaults.";
    styleProtein();
  }

  // ------------------------------------------------------------------ two small tasks per morph
  // q1: run the pocket finder on the start frame. q2: find the first frame where the ligand fits.
  const taskKey = (q) => `morph:${S.data.id}:${q}`;
  function taskDone(q) {
    if (PW.progress.tasks[taskKey(q)]) return;
    PW.progress.tasks[taskKey(q)] = true; PW.save();
    PW.toast("Task done.", "good", 2000);
    renderTasks();
  }
  // first frame (in % steps) with at most FIT_CLASHES clashes, computed once per morph
  const FIT_CLASHES = 2;
  function firstFit() {
    for (let pct = 0; pct <= 100; pct++) { S.xyz = coordsAt(pct / 100); if (clashes() <= FIT_CLASHES) return pct; }
    return null;
  }
  function renderTasks() {
    const box = document.getElementById("mo-tasks");
    if (!box || !S.data) return;
    const [startRank] = S.summary, f = S.fit;
    const t = (q, text, ans) => `<div class="task ${PW.progress.tasks[taskKey(q)] ? "done" : ""}"><span class="tick"></span><div><p>${text}</p>${PW.progress.tasks[taskKey(q)] ? `<div class="ans">${ans}</div>` : ""}</div></div>`;
    box.innerHTML = t("q1", "Slide to <b>0%</b> and press <b>Find pockets</b>. Does the pocket finder see a pocket where the ligand will bind?",
      `At the start the grid ranks the ligand site ${startRank === "not found" ? "<b>nowhere</b>: no pocket is centred there" : `<b>${startRank}</b>`}. A pocket finder only sees the shape it is given.`)
      + (S.ligXyz.length ? t("q2", `Find the <b>first frame</b> where the ligand fits: at most ${FIT_CLASHES} clashes.`,
      f === null ? "It never fits along this straight-line path." : f <= 10
        ? `Already at ${f}%: the open cleft has room for the ligand, it is just too open to count as a pocket. Binding needs the closing, not the space.`
        : `At ${f}%. Before that the ligand overlaps the protein: a structure from earlier on the path would reject the right answer.`) : "")
      + `<div class="task"><span class="tick" style="opacity:.4"></span><div><p>You only have the start structure and want to dock this ligand. What could you do?</p>
        <details class="small"><summary>One answer</summary>Let side chains or loops move during docking, dock into several structures (an ensemble, e.g. from MD or other crystal structures), or predict the complex with a co-folding model (Chapter 5). Slide 32.</details></div></div>`;
  }

  // endpoint summary, computed once per morph
  function endpointSummary() {
    return [0, 1].map((t) => {
      const res = Pockets.gridPockets({ ...S.baseAtoms, xyz: coordsAt(t) }, Pockets.GRID_DEFAULTS);
      const k = res.pockets.findIndex((p) => Pockets.evaluate(p, S.ligXyz).dcc <= Pockets.DCC_SUCCESS_A);
      return k >= 0 ? `#${k + 1} (${res.pockets[k].volume.toFixed(0)} Å³)` : "not found";
    });
  }

  function play() {
    if (S.playing) { S.playing = false; $("#mo-play").textContent = "▶ Play"; return; }
    S.playing = true;
    $("#mo-play").textContent = "❚❚ Pause";
    clearPockets();
    let dir = S.t >= 1 ? -1 : 1, last = performance.now();
    const step = (now) => {
      if (!S || !S.playing) return;
      const dt = (now - last) / 1000; last = now;
      let t = S.t + dir * dt / PLAY_SECONDS;
      if (t >= 1) { t = 1; dir = -1; } else if (t <= 0) { t = 0; dir = 1; }
      setFrame(t);
      requestAnimationFrame(step);
    };
    requestAnimationFrame(step);
  }

  function renderPanels() {
    const d = S.data;
    $("#left").innerHTML = `<h3>Apo to holo</h3>
      <p class="small">A structure is one snapshot. Many proteins change shape when a ligand binds (slide 31), and a pocket finder or a docking program only ever sees the snapshot you give it.</p>
      <div class="levels" id="mo-list"></div>
      <h3>Morph</h3>
      <div class="row"><button id="mo-play" class="primary">▶ Play</button><span class="chip" id="mo-t">0%</span></div>
      <div id="mo-slider"></div>
      <div class="row small" style="justify-content:space-between;align-items:flex-start;flex-wrap:nowrap"><span class="muted" style="max-width:48%">${esc(d.start.label)}</span><span class="muted" style="max-width:48%;text-align:right">${esc(d.end.label)}</span></div>
      <h3>Colour</h3><div id="mo-color"></div>
      <h3>Show</h3><div id="mo-show"></div>
      <p class="small muted">The frames in between are a straight line from one crystal structure to the other. Real proteins do not move like this; only the two end points are measured.</p>`;
    S.index.morphs.forEach((m) => {
      const b = el(`<button class="lvl ${m.id === d.id ? "on" : ""}"><span class="n">↔</span><span class="t">${esc(m.title)}<small>CA RMSD ${m.ca_rmsd_A} Å</small></span><span></span></button>`);
      b.addEventListener("click", () => openMorph(m.id));
      $("#mo-list").appendChild(b);
    });
    S.slider = slider({ label: "Position", min: 0, max: 100, step: 1, value: Math.round(S.t * 100), unit: "%",
      onInput: (v) => {
        S.playing = false; $("#mo-play").textContent = "▶ Play";
        if (S.dotModels && S.dotModels.length) clearPockets(); // the dots belong to the old frame
        setFrame(v / 100);
        if (S.fit !== null && v === S.fit) taskDone("q2");
      } });
    $("#mo-slider").appendChild(S.slider);
    $("#mo-play").addEventListener("click", play);
    $("#mo-color").appendChild(segmented([{ value: "motion", label: "How far it moves" }, { value: "plain", label: "Plain" }], S.colorBy, (v) => { S.colorBy = v; styleProtein(); const lg = $("#mo-legend"); if (lg) lg.classList.toggle("hidden", v !== "motion"); }));
    $("#mo-show").appendChild(checkbox("Ligand (from the end structure)", S.showLigand, (v) => { S.showLigand = v; styleProtein(); }));
    $("#mo-show").appendChild(checkbox("Surface (slow while playing)", S.showSurface, (v) => { S.showSurface = v; if (!v) S.viewer.removeAllSurfaces(); styleProtein(); }));

    const [startRank, endRank] = S.summary;
    $("#right").innerHTML = `<h3>${esc(d.title)}</h3>
      <div class="story">${d.story}</div>
      <table style="margin-top:10px"><thead><tr><th></th><th>Start</th><th>End</th></tr></thead><tbody>
        <tr><td>Structure</td><td><a href="https://www.rcsb.org/structure/${d.start.pdb_id}" target="_blank" rel="noopener">${d.start.pdb_id}</a></td><td><a href="https://www.rcsb.org/structure/${d.end.pdb_id}" target="_blank" rel="noopener">${d.end.pdb_id}</a></td></tr>
        <tr><td>State</td><td>${esc(d.start.label)}</td><td>${esc(d.end.label)}</td></tr>
        <tr><td>Grid pocket at the ligand site</td><td>${startRank}</td><td>${endRank}</td></tr></tbody></table>
      <p class="small muted">Superposed on the ${d.core_ca} of ${d.aligned_ca} CA atoms that stay within ${d.core_trim_A} Å (the rigid core). Overall CA RMSD ${d.ca_rmsd_A} Å.</p>
      <div class="box acc"><div class="readout"><span>Position</span><b id="mo-t2"></b>
        <span title="Ligand-protein heavy-atom pairs closer than ${S.clashA} Å, with the ligand held where it sits in the end structure">Ligand clashes at this frame</span><b id="mo-clash">–</b></div>
        <div class="row" style="margin-top:8px"><button id="mo-find">Find pockets at this frame</button><button id="mo-clear" class="ghost">Clear</button></div>
        <p class="small" id="mo-pocket">Grid method, Chapter 2 defaults.</p></div>
      <h3>Tasks</h3><div id="mo-tasks"></div>
      <p class="small muted">Start: ${citation(d.start.citation)}<br>End: ${citation(d.end.citation)}</p>`;
    $("#mo-find").addEventListener("click", findPocketsNow);
    renderTasks();
    $("#mo-clear").addEventListener("click", clearPockets);
  }

  async function openMorph(id) {
    if (S.playing) S.playing = false;
    PW.setArg(id);
    const d = await PW.get("morph", id);
    const v = PW.resetView();
    Object.assign(S, { data: d, viewer: v, t: 0, dotModels: [], pocketLabels: [] });
    PW.setHud("Apo to holo", d.title);
    S.model = v.addModel(d.pdb, "pdb");
    S.atoms3d = S.model.selectedAtoms({});
    S.baseAtoms = Pockets.atomsFromPdb(d.pdb);
    const lig = d.ligand_pdb ? Pockets.atomsFromPdb(d.ligand_pdb) : { n: 0, xyz: [] };
    S.ligXyz = [];
    for (let i = 0; i < lig.n; i++) S.ligXyz.push([lig.xyz[3 * i], lig.xyz[3 * i + 1], lig.xyz[3 * i + 2]]);
    S.ligModel = d.ligand_pdb ? v.addModel(d.ligand_pdb, "pdb") : null;
    S.dispColors = S.atoms3d.map((a, i) => dispColor(Math.hypot(d.end_xyz[3 * i] - d.start_xyz[3 * i], d.end_xyz[3 * i + 1] - d.start_xyz[3 * i + 1], d.end_xyz[3 * i + 2] - d.start_xyz[3 * i + 2])));
    S.summary = endpointSummary();
    S.fit = S.ligXyz.length ? firstFit() : null;
    renderPanels();
    setFrame(0);
    v.zoomTo();
    v.render();
    $("#stage-overlay").innerHTML = "";
    $("#stage-overlay").appendChild(el(`<div class="legend-card${S.colorBy === "motion" ? "" : " hidden"}" id="mo-legend">Colour: how far each atom moves between the two structures, after superposing the rigid core.
      <div style="height:8px;border-radius:6px;margin:5px 0 2px;background:linear-gradient(90deg,${dispColor(0)},${dispColor(DISP_MAX_A / 2)},${dispColor(DISP_MAX_A)})"></div>
      <div class="row small" style="margin:0"><span>0 Å</span><span style="flex:1"></span><span>≥ ${DISP_MAX_A} Å</span></div></div>`));
  }

  PW.chapter({
    id: "morph", title: "Apo to holo", color: "#ffd166", num: "Chapter 4",
    blurb: "Proteins move. Morph between the empty and the ligand-bound crystal structure and watch the pocket appear.",
    learn: "Conformational change · induced fit · why docking into the wrong structure fails",
    async enter(arg) {
      const index = await PW.get("index", "index");
      S = { index, colorBy: "motion", showLigand: true, showSurface: false, playing: false };
      S.clashA = index.docking.clash_cutoff_A;
      await openMorph(index.morphs.some((m) => m.id === arg) ? arg : index.morphs[0].id);
    },
    progress(index) {
      const done = index.morphs.filter((m) => PW.progress.tasks[`morph:${m.id}:q1`] && PW.progress.tasks[`morph:${m.id}:q2`]).length;
      return { text: `${done}/${index.morphs.length}`, done: done === index.morphs.length };
    },
    tour: () => [
      { sel: "#mo-list", title: "Three proteins that move", text: "Each morph goes from one crystal structure to another; the frames in between are a straight line, not a real pathway." },
      { sel: "#mo-slider", title: "Slide through the change", text: "Or press Play. The colour shows how far each part moves." },
      { sel: "#mo-find", title: "Run the pocket finder", text: "The grid method from Chapter 2, on the frame you are looking at." },
      { sel: "#mo-tasks", title: "Tasks", text: "They tick themselves; the answer appears when they do." },
    ],
    leave() { if (S) S.playing = false; S = null; },
    recentre() { if (S) { S.viewer.zoomTo(); S.viewer.render(); } },
  });
})();
