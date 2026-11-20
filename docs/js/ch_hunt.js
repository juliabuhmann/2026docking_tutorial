// Chapter 1: Pocket hunt. Click where you think the ligand binds; ten real proteins, easy to hard.
(function () {
  "use strict";
  const { $, esc, el, stars, citation, toast, segmented, slider, vec } = PW;

  // Feedback distances (A from the clicked atom to the nearest ligand atom). Feel only.
  const HOT_A = 8, WARM_A = 15;
  // Surface view: how far beyond the van der Waals radius a click still counts as on an atom. About
  // the solvent probe radius (1.4 A), which is how far the molecular surface can sit from the atoms.
  const SURF_PICK_PAD_A = 1.2;
  const CHAIN_COLORS = ["#8fb8ff", "#ffb38a", "#b4f0a8", "#e3b6ff"];

  let S = null; // state of the open level

  // ------------------------------------------------------------------ picking that respects the slice
  // 3Dmol's own picking ignores the clipping plane, so atoms that are sliced away could be hit.
  // Instead: project every atom, drop the clipped ones, take the frontmost atom under the cursor.
  function quatRotate(q, v) {
    const [x, y, z, w] = q;
    const ix = w * v[0] + y * v[2] - z * v[1], iy = w * v[1] + z * v[0] - x * v[2];
    const iz = w * v[2] + x * v[1] - y * v[0], iw = -x * v[0] - y * v[1] - z * v[2];
    return [ix * w + iw * -x + iy * -z - iz * -y, iy * w + iw * -y + iz * -x - ix * -z, iz * w + iw * -z + ix * -y - iy * -x];
  }
  // Depth of each atom along the view axis, relative to the rotation centre (positive = towards you).
  function viewDepths(viewer, xyz) {
    const v = viewer.getView(), q = [v[4], v[5], v[6], v[7]], out = new Float64Array(xyz.length / 3);
    for (let i = 0; i < out.length; i++) out[i] = quatRotate(q, [xyz[3 * i] + v[0], xyz[3 * i + 1] + v[1], xyz[3 * i + 2] + v[2]])[2];
    return out;
  }
  function pixelsPerAngstrom(viewer, p) {
    const [a, b] = viewer.modelToScreen([vec(p), { x: p[0] + 1, y: p[1], z: p[2] }]);
    const [c] = viewer.modelToScreen([{ x: p[0], y: p[1] + 1, z: p[2] }]);
    return Math.max(Math.hypot(b.x - a.x, b.y - a.y), Math.hypot(c.x - a.x, c.y - a.y));
  }
  function pickAtom(viewer, atoms, radii, pageX, pageY, cutDepth) {
    const pts = [];
    for (let i = 0; i < atoms.n; i++) pts.push({ x: atoms.xyz[3 * i], y: atoms.xyz[3 * i + 1], z: atoms.xyz[3 * i + 2] });
    const scr = viewer.modelToScreen(pts);
    const depth = viewDepths(viewer, atoms.xyz);
    const ppa = pixelsPerAngstrom(viewer, [atoms.xyz[0], atoms.xyz[1], atoms.xyz[2]]);
    let best = -1;
    for (let i = 0; i < atoms.n; i++) {
      if (depth[i] > cutDepth) continue;
      const r = radii[i] * ppa;
      if ((scr[i].x - pageX) ** 2 + (scr[i].y - pageY) ** 2 > r * r) continue;
      if (best < 0 || depth[i] > depth[best]) best = i;
    }
    return best;
  }

  // ------------------------------------------------------------------ level list
  function unlocked(index, i) {
    if (PW.progress.unlockAll || i === 0) return true;
    return PW.progress.stars[`hunt:${index.levels[i - 1].id}`] !== undefined;
  }

  function renderLeft(index, cur) {
    const left = $("#left");
    left.innerHTML = `<h3>Pocket hunt</h3>
      <p class="small">Each protein binds a small molecule somewhere. Rotate it, look for a pocket, and <b>click on the protein where you think the ligand sits</b>. First try: ★★★.</p>
      <div class="levels"></div>
      <h3>Views</h3>
      <div id="hunt-views"></div>
      <div id="hunt-slice"></div>
      <p class="small muted">Some pockets are buried. Slice into the protein to look inside; you can click on what the slice reveals.</p>`;
    const list = $(".levels", left);
    index.levels.forEach((L, i) => {
      const s = PW.progress.stars[`hunt:${L.id}`];
      const open = unlocked(index, i);
      const b = el(`<button class="lvl ${i === cur ? "on" : ""} ${open ? "" : "locked"}" ${open ? "" : "disabled"}>
        <span class="n">${i + 1}</span>
        <span class="t">${esc(L.title)}<small><span class="diff ${L.difficulty}">${L.difficulty}</span> PDB ${L.pdb_id}</small></span>
        ${s !== undefined ? stars(s) : `<span class="stars"><span class="off">★★★</span></span>`}</button>`);
      if (open) b.addEventListener("click", () => openLevel(index, i));
      list.appendChild(b);
    });
  }

  // ------------------------------------------------------------------ drawing
  function styleProtein() {
    const v = S.viewer, m = S.model;
    const chains = [...new Set(S.atoms.chain)];
    const colorOf = (c) => CHAIN_COLORS[chains.indexOf(c) % CHAIN_COLORS.length];
    const shown = S.solved || S.revealed;
    m.setStyle({}, {});
    if (S.mode === "surface") {
      // once the answer is shown: pocket-lining residues gold, surface see-through so the ligand shows
      const cf = (a) => (shown && S.pocketKeys.has(`${a.chain}:${a.resi}`) ? "#ffd166" : colorOf(a.chain));
      const op = shown ? 0.62 : 1;
      if (!S.surf) S.surf = PW.surface(v, m, cf, op);
      else S.surf.restyle(cf, op);
    } else if (S.surf) { S.surf.remove(); S.surf = null; }
    // a sliced surface is a hollow shell; fill it with the atoms so the cut face looks solid
    const atomsToo = S.mode === "atoms" || (S.mode === "surface" && S.slice > 0);
    S.styledSliced = S.slice > 0;
    chains.forEach((c) => {
      if (S.mode === "cartoon") {
        m.setStyle({ chain: c }, { cartoon: { color: colorOf(c), opacity: 1 } });
      } else if (atomsToo) {
        m.setStyle({ chain: c }, { sphere: { colorscheme: { prop: "elem", map: { C: colorOf(c), N: "#6d8df0", O: "#ee7b7b", S: "#f0cd62" } }, scale: 1 } });
      }
    });
    if (S.mode !== "surface" || atomsToo) m.setStyle({ hetflag: true }, { sphere: { color: "#c0c7d4", scale: 0.9 } });
    if (shown && S.mode !== "atoms") {
      S.prot.pocket_residues.forEach((r) => {
        m.setStyle({ chain: r.chain, resi: r.resi }, { stick: { radius: 0.18, colorscheme: "whiteCarbon" } }, true);
      });
    }
    v.render();
  }

  function drawLigand() {
    if (S.ligModel) S.viewer.removeModel(S.ligModel);
    S.ligModel = null;
    if (!(S.solved || S.revealed)) return;
    S.ligModel = S.viewer.addModel(S.prot.ligand.molblock, "sdf");
    S.ligModel.setStyle({}, { stick: { radius: 0.28, colorscheme: "greenCarbon" }, sphere: { scale: 0.25, colorscheme: "greenCarbon" } });
  }

  function drawProbes() {
    S.viewer.removeAllShapes();
    S.probes.forEach((p) => {
      S.viewer.addSphere({ center: vec(p.xyz), radius: p.ok ? 1.1 : 0.9, color: p.color, opacity: 0.95 });
    });
    if (S.hintPockets) {
      S.hintPockets.forEach((pk, i) => {
        const col = ["#ffd166", "#5eead4", "#b69cff"][i];
        pk.points.forEach((q, j) => { if (j % 2 === 0) S.viewer.addSphere({ center: vec(q), radius: 0.35, color: col, opacity: 0.8 }); });
        S.viewer.addLabel(`#${i + 1}`, { position: vec(pk.centre), fontSize: 13, fontColor: "#06131c", backgroundColor: col, backgroundOpacity: 0.95, inFront: true });
      });
    }
    S.viewer.render();
  }

  function applySlice() {
    // slider value = how far behind the front of the protein the cut sits, 0 = no cut
    const v = S.viewer;
    if (S.slice <= 0) { v.setSlab(-S.radius * 3, S.radius * 3); S.cutDepth = Infinity; }
    else {
      const near = -S.radius + S.slice;
      v.setSlab(near, S.radius * 3);
      S.cutDepth = -near;
    }
    v.render();
  }

  // ------------------------------------------------------------------ the right panel
  function renderRight() {
    const L = S.level, right = $("#right");
    const tries = S.tries, maxStars = S.hintLevel >= 3 ? 1 : 3;
    let html = `<h3>Level ${S.i + 1} of ${S.index.levels.length}</h3>
      <div class="big" style="font-size:22px">${esc(L.title)}</div>
      <p class="small muted"><span class="diff ${L.difficulty}">${L.difficulty}</span> &nbsp;PDB <a href="https://www.rcsb.org/structure/${L.pdb_id}" target="_blank" rel="noopener">${L.pdb_id}</a> · ${L.chains.length > 1 ? `${L.chains.length} chains` : "1 chain"} · ${L.resolution ? `${L.resolution} Å` : ""}</p>
      ${L.intro && L.intro.length ? `<button id="hunt-about" class="ghost small"><svg width="13" height="13" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.6" style="vertical-align:-2px;margin-right:5px"><circle cx="8" cy="8" r="6.6"/><path d="M8 7v4.4" stroke-linecap="round"/><circle cx="8" cy="4.6" r=".9" fill="currentColor" stroke="none"/></svg>About this protein</button>` : ""}`;
    if (!S.solved && !S.revealed) {
      html += `<div class="box acc"><b>Task:</b> click on the protein where ${esc(L.ligand_label)} binds.
        <div class="small muted" style="margin-top:4px">Tries so far: <b>${tries}</b>. Solve on the first try for ${stars(maxStars)}${maxStars < 3 ? " (capped by the hint)" : ""}.</div></div>
        <div id="hunt-log"></div>
        <div class="row"><button id="hunt-hint">Hint ${S.hintLevel + 1 <= 3 ? `(${S.hintLevel + 1}/3)` : ""}</button><button id="hunt-reveal" class="ghost">Show me</button></div>
        <div id="hunt-hints"></div>`;
    } else {
      const s = PW.progress.stars[`hunt:${L.id}`];
      html += `<div class="box gold"><div class="row" style="justify-content:space-between"><b>${S.solved ? "Found it!" : "Here it is."}</b>${s !== undefined ? stars(s) : ""}</div>
        <div class="small muted">${S.solved ? `${tries} ${tries === 1 ? "try" : "tries"}.` : "Revealed: no stars, but you can come back."} The ligand is shown in green, the residues lining the pocket as sticks.</div></div>
        <h3>The story</h3>
        <div class="story">${L.story}${L.caveat ? `<p class="small muted">${esc(L.caveat)}</p>` : ""}
          <div class="ref">Structure: ${citation(L.citation)}${L.story_refs ? `<br>Story: ${esc(L.story_refs)}` : ""}</div></div>
        <div class="row" style="margin-top:12px">${S.i + 1 < S.index.levels.length ? `<button class="primary" id="hunt-next">Next protein →</button>` : `<button class="primary" id="hunt-algo">On to the algorithms →</button>`}
        <button id="hunt-retry" class="ghost">Try again</button></div>`;
    }
    right.innerHTML = html;
    renderLog();
    renderHints();
    const q = (id) => document.getElementById(id);
    if (q("hunt-about")) q("hunt-about").addEventListener("click", () => intro(true));
    if (q("hunt-hint")) q("hunt-hint").addEventListener("click", nextHint);
    if (q("hunt-reveal")) q("hunt-reveal").addEventListener("click", reveal);
    if (q("hunt-next")) q("hunt-next").addEventListener("click", () => openLevel(S.index, S.i + 1));
    if (q("hunt-algo")) q("hunt-algo").addEventListener("click", () => PW.open("algo"));
    if (q("hunt-retry")) q("hunt-retry").addEventListener("click", () => openLevel(S.index, S.i));
  }

  function renderLog() {
    const box = document.getElementById("hunt-log");
    if (!box) return;
    box.innerHTML = S.log.slice(-5).reverse().map((m) => `<div class="box small" style="border-color:${m.color}">${m.html}</div>`).join("");
  }

  const HINTS = [
    "Pockets are concave: a groove, a cleft or a hole that a small molecule could sit in. Turn the protein around before you click.",
    "Use the <b>Surface</b> view to see the real shape of the protein. If nothing looks deep enough, use the <b>slice</b> slider: some pockets are completely buried.",
    "Ask an algorithm: the three largest pockets found by the grid method from Chapter 2 are now shown as dots. Is the biggest one always the right one? (This hint caps the level at ★.)",
  ];
  function renderHints() {
    const box = document.getElementById("hunt-hints");
    if (!box) return;
    box.innerHTML = HINTS.slice(0, S.hintLevel).map((h, i) => `<p class="small"><b>Hint ${i + 1}.</b> ${h}</p>`).join("");
  }
  function nextHint() {
    if (S.hintLevel >= 3) return;
    S.hintLevel++;
    if (S.hintLevel === 3) {
      const res = Pockets.gridPockets(S.atoms, Pockets.GRID_DEFAULTS);
      S.hintPockets = res.pockets.slice(0, 3);
      drawProbes();
    }
    renderRight();
  }

  // ------------------------------------------------------------------ clicking
  function onPick(i) {
    if (S.solved || S.revealed || i < 0) return;
    const p = [S.atoms.xyz[3 * i], S.atoms.xyz[3 * i + 1], S.atoms.xyz[3 * i + 2]];
    const chain = S.atoms.chain[i], resi = S.atoms.resi[i], resn = S.atoms.resn[i];
    const d = Math.min(...S.ligXyz.map((l) => Pockets.dist(l, p)));
    S.tries++;
    const inPocket = S.pocketKeys.has(`${chain}:${resi}`);
    const decoy = S.prot.decoys.find((dc) => dc.residues.some((r) => r.chain === chain && r.resi === resi));
    const where = `${resn} ${resi}${S.level.chains.length > 1 ? ` (chain ${chain})` : ""}`;
    if (inPocket) {
      S.solved = true;
      const n = S.hintLevel >= 3 ? 1 : S.tries === 1 ? 3 : S.tries <= 3 ? 2 : 1;
      const key = `hunt:${S.level.id}`;
      PW.progress.stars[key] = Math.max(n, PW.progress.stars[key] || 0);
      PW.save();
      S.probes.push({ xyz: p, color: "#7ee787", ok: true });
      toast(`<b>${where}</b> lines the pocket. ${stars(n)}`, "good");
      drawProbes(); drawLigand(); styleProtein(); renderRight(); renderLeft(S.index, S.i);
      return;
    }
    let color, html;
    if (decoy) {
      color = "#ffb454";
      html = `<b>${where}</b>: you found ${esc(decoy.name)}. ${esc(decoy.feedback)}`;
    } else if (d <= HOT_A) {
      color = "#ff8f8f"; html = `<b>${where}</b>: hot! The ligand is ${d.toFixed(0)} Å away.`;
    } else if (d <= WARM_A) {
      color = "#ffb454"; html = `<b>${where}</b>: warm. ${d.toFixed(0)} Å from the ligand.`;
    } else {
      color = "#60a5fa"; html = `<b>${where}</b>: cold. ${d.toFixed(0)} Å from the ligand.`;
    }
    S.probes.push({ xyz: p, color, ok: false });
    S.log.push({ html, color });
    drawProbes(); renderRight();
  }

  function reveal() {
    S.revealed = true;
    drawLigand(); styleProtein(); renderRight();
    S.viewer.zoomTo({ model: S.model, or: S.prot.pocket_residues.map((r) => ({ chain: r.chain, resi: r.resi })) }, 600);
    applySlice();
  }

  function bindClicks() {
    const stage = $("#stage");
    let down = null;
    S.onDown = (e) => { down = { x: e.pageX, y: e.pageY, t: Date.now() }; };
    S.onUp = (e) => {
      if (!down || e.button > 0) return;
      const moved = Math.hypot(e.pageX - down.x, e.pageY - down.y);
      if (moved < 5 && Date.now() - down.t < 600) {
        const radii = S.mode === "cartoon" ? S.cartoonRadii : S.mode === "surface" ? S.surfRadii : S.vdwRadii;
        onPick(pickAtom(S.viewer, S.atoms, radii, e.pageX, e.pageY, S.cutDepth));
      }
      down = null;
    };
    stage.addEventListener("pointerdown", S.onDown);
    stage.addEventListener("pointerup", S.onUp);
  }
  function unbindClicks() {
    if (!S) return;
    $("#stage").removeEventListener("pointerdown", S.onDown);
    $("#stage").removeEventListener("pointerup", S.onUp);
  }

  // ------------------------------------------------------------------ open a level
  async function openLevel(index, i, first) {
    unbindClicks();
    PW.closeCards();
    const L = index.levels[i];
    PW.setArg(i + 1);
    PW.setHud("Pocket hunt", `Level ${i + 1}: ${L.title}`);
    const prot = await PW.get("prot", L.id);
    const viewer = PW.resetView();
    const atoms = Pockets.atomsFromPdb(prot.receptor_pdb);
    S = {
      index, i, level: L, prot, viewer, atoms,
      model: viewer.addModel(prot.receptor_pdb, "pdb"),
      ligXyz: prot.ligand.atoms.map((a) => a.xyz),
      pocketKeys: new Set(prot.pocket_residues.map((r) => `${r.chain}:${r.resi}`)),
      vdwRadii: atoms.elem.map(Pockets.vdw),
      // in cartoon view the backbone is what you see; accept clicks within 2.2 A of any atom
      cartoonRadii: atoms.elem.map(() => 2.2),
      // the molecular surface bridges the gaps between atoms; accept clicks a little outside them
      surfRadii: atoms.elem.map((e) => Pockets.vdw(e) + SURF_PICK_PAD_A),
      surf: null,
      mode: S ? S.mode : "surface", slice: 0, cutDepth: Infinity,
      tries: 0, probes: [], log: [], hintLevel: 0, hintPockets: null, solved: false, revealed: false, ligModel: null,
    };
    // radius of the protein around its centre, for the slice slider
    const pts = [];
    for (let a = 0; a < atoms.n; a++) pts.push([atoms.xyz[3 * a], atoms.xyz[3 * a + 1], atoms.xyz[3 * a + 2]]);
    const c = Pockets.centroid(pts);
    let r = 0;
    for (const p of pts) r = Math.max(r, Pockets.dist(c, p));
    S.radius = r;

    renderLeft(index, i);
    const views = $("#hunt-views");
    views.appendChild(segmented([
      { value: "surface", label: "Surface", title: "Molecular surface: the shape a ligand actually sees" },
      { value: "atoms", label: "Atoms", title: "Every atom as a sphere, coloured by element" },
      { value: "cartoon", label: "Cartoon", title: "The backbone fold: helices and strands" }], S.mode, (m) => { S.mode = m; styleProtein(); }));
    $("#hunt-slice").appendChild(slider({
      label: "Slice into the protein", min: 0, max: Math.round(r * 1.4), step: 1, value: 0, unit: " Å",
      onInput: (val) => {
        S.slice = val;
        if (S.mode === "surface" && (val > 0) !== S.styledSliced) styleProtein();
        applySlice();
      },
    }));
    styleProtein();
    viewer.zoomTo();
    applySlice();
    viewer.render();
    renderRight();
    bindClicks();
    const legend = el(`<div class="legend-card"><b>Click</b> the protein to probe. <span style="color:#ff8f8f">hot</span> ≤ ${HOT_A} Å,
      <span style="color:#ffb454">warm</span> ≤ ${WARM_A} Å, <span style="color:#60a5fa">cold</span> further from the ligand.</div>`);
    $("#stage-overlay").innerHTML = "";
    $("#stage-overlay").appendChild(legend);
    if (!first) intro(false);
  }

  // ------------------------------------------------------------------ short intro before each level
  // Automatic the first time a level is opened in this browser; the "About" button replays it.
  function intro(force) {
    const L = S.level, key = `intro:hunt:${L.id}`;
    if (!L.intro || !L.intro.length || (!force && PW.seen(key))) return Promise.resolve();
    PW.markSeen(key);
    const n = L.intro.length;
    const facts = `<div class="facts">${esc(L.ligand_name)} · ${L.n_ligand_atoms} heavy atoms · PDB ${esc(L.pdb_id)}</div>`;
    return PW.cards(L.intro.map((c, j) => ({
      kicker: `Level ${S.i + 1} · ${esc(L.title)} · ${j + 1}/${n}`,
      title: c.title, html: c.text + (c.title === "The ligand" ? facts : ""),
    })), { done: S.solved || S.revealed ? "Close" : "Start hunting" });
  }

  PW.chapter({
    id: "hunt", title: "Pocket hunt", color: "#5eead4", num: "Chapter 1",
    blurb: "Ten real proteins, from wide-open clefts to pockets nobody could see for 40 years. Find where the drug binds.",
    learn: "Pockets · surfaces · buried cavities · allosteric and cryptic sites",
    progress(index) {
      const done = index.levels.filter((L) => PW.progress.stars[`hunt:${L.id}`] !== undefined).length;
      return { text: `${done}/${index.levels.length}`, done: done === index.levels.length };
    },
    async enter(arg) {
      const index = await PW.get("index", "index");
      let i = Math.max(0, Math.min(index.levels.length - 1, (parseInt(arg, 10) || 1) - 1));
      if (!unlocked(index, i)) i = 0;
      await openLevel(index, i, true);
    },
    afterTour() { if (S) return intro(false); },
    tour: () => [
      { sel: "#left .levels", title: "Ten levels", text: "From easy to hard. The next one unlocks when you solve (or reveal) the current one. Stars: ★★★ on the first click." },
      { sel: "#hunt-views", title: "Three views", text: "<b>Surface</b>: the shape a ligand sees, best for spotting pockets. <b>Atoms</b>: every atom as a sphere. <b>Cartoon</b>: the fold." },
      { sel: "#hunt-slice", title: "Slice", text: "Cuts away the front of the protein so you can look (and click) inside. Some pockets are completely buried." },
      { sel: "#stage", title: "Click to probe", text: "A short click on the protein is a guess. The mark is coloured hot, warm or cold by its distance to the ligand." },
      { sel: "#hunt-hint", title: "Stuck?", text: "Three hints, the last one runs a pocket finder (and caps the level at ★). <b>Show me</b> reveals the answer." },
      { sel: "#hunt-about", title: "About this protein", text: "Every level starts with a few short cards on the protein and its ligand. Reopen them here." },
    ],
    leave() { unbindClicks(); S = null; },
    recentre() { if (S) { S.viewer.zoomTo(); applySlice(); } },
  });
})();
