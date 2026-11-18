// Chapter 5: co-folding. (1) An illustration of the diffusion idea, (2) a hand-off to LocalFold,
// which runs AlphaFold 3 and similar models in the browser, (3) scoring the student's prediction
// against the crystal structure (js/compare.js, tested in tests/test_compare.js).
(function () {
  "use strict";
  const { $, esc, el, citation, toast, vec } = PW;
  const FFLATE_URL = "https://cdn.jsdelivr.net/npm/fflate@0.8.2/umd/index.js";
  const LOCALFOLD_URL = "https://localfold.org";
  // Illustration only. The shape follows the EDM / AlphaFold 3 schedule (rho = 7, Karras et al.
  // 2022), but the largest noise is shrunk from thousands of Angstrom to something you can watch.
  const ILLU_SIGMA_MAX_A = 30, ILLU_SIGMA_MIN_A = 0.05, ILLU_RHO = 7, ILLU_STEPS = 60, ILLU_SECONDS = 6;
  const THREE = { A: "ALA", R: "ARG", N: "ASN", D: "ASP", C: "CYS", Q: "GLN", E: "GLU", G: "GLY", H: "HIS", I: "ILE", L: "LEU", K: "LYS", M: "MET", F: "PHE", P: "PRO", S: "SER", T: "THR", W: "TRP", Y: "TYR", V: "VAL" };

  let S = null;

  // ------------------------------------------------------------------ illustration
  function gauss(rng) { // Box-Muller with a seeded generator, so every play looks the same
    const u = Math.max(rng(), 1e-12), v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }
  function mulberry32(a) { return () => { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
  function sigmaAt(f) { // f: 0 (start, pure noise) -> 1 (end)
    const a = ILLU_SIGMA_MAX_A ** (1 / ILLU_RHO), b = ILLU_SIGMA_MIN_A ** (1 / ILLU_RHO);
    return (a + f * (b - a)) ** ILLU_RHO;
  }

  function buildIllustration() {
    const d = S.cof, p = S.prot, lines = [];
    const keys = Object.keys(d.ca_by_seq).map(Number).sort((a, b) => a - b);
    keys.forEach((k, i) => {
      const [x, y, z] = d.ca_by_seq[k];
      lines.push(`ATOM  ${String(i + 1).padStart(5)}  CA  ${THREE[d.sequence[k - 1]] || "UNK"} A${String(k).padStart(4)}    ${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}  1.00  0.00           C`);
    });
    p.ligand.atoms.forEach((a, i) => {
      const [x, y, z] = a.xyz;
      lines.push(`HETATM${String(keys.length + i + 1).padStart(5)} ${a.name.padEnd(4)} ${d.ligand_id.padStart(3)} B   1    ${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}  1.00  0.00          ${a.elem.padStart(2)}`);
    });
    const v = PW.resetView();
    S.viewer = v;
    S.illuModel = v.addModel(lines.join("\n") + "\nEND\n", "pdb");
    S.illuAtoms = S.illuModel.selectedAtoms({});
    S.illuRef = S.illuAtoms.map((a) => [a.x, a.y, a.z]);
    const rng = mulberry32(42);
    S.illuNoise = S.illuAtoms.map(() => [gauss(rng), gauss(rng), gauss(rng)]);
    S.nCa = keys.length;
    illuFrame(1);
    v.zoomTo();
    v.render();
  }

  function illuFrame(f) {
    const sigma = f >= 1 ? 0 : sigmaAt(f);
    S.illuAtoms.forEach((a, i) => {
      a.x = S.illuRef[i][0] + sigma * S.illuNoise[i][0];
      a.y = S.illuRef[i][1] + sigma * S.illuNoise[i][1];
      a.z = S.illuRef[i][2] + sigma * S.illuNoise[i][2];
    });
    const m = S.illuModel, n = S.cof.sequence.length;
    // N terminus blue to C terminus violet, so the chain can be followed
    const tone = (a) => {
      const f = Math.min(1, Math.max(0, (a.resi - 1) / Math.max(1, n - 1)));
      const c0 = [86, 132, 230], c1 = [196, 140, 255];
      return `rgb(${c0.map((v, k) => Math.round(v + (c1[k] - v) * f)).join(",")})`;
    };
    m.setStyle({ chain: "A" }, { sphere: { radius: 0.6, colorfunc: tone } });
    m.setStyle({ chain: "A" }, { cartoon: { style: "trace", thickness: 0.3, colorfunc: tone, opacity: 0.85 } }, true);
    // bonds were assigned from the crystal geometry; while the atoms are scattered they would
    // stretch across the picture, so the ligand is drawn as atoms only until it has nearly settled
    m.setStyle({ chain: "B" }, sigma > 1 ? { sphere: { scale: 0.45, colorscheme: "greenCarbon" } }
      : { sphere: { scale: 0.32, colorscheme: "greenCarbon" }, stick: { radius: 0.24, colorscheme: "greenCarbon" } });
    S.viewer.render();
    const step = Math.round(f * ILLU_STEPS);
    const o = document.getElementById("cf-step");
    if (o) o.innerHTML = `step <b>${step}</b> / ${ILLU_STEPS} · noise σ = <b>${sigma.toFixed(sigma < 1 ? 2 : 0)} Å</b>`;
  }

  function playIllustration() {
    if (S.view !== "illu") showIllustration();
    S.illuPlaying = true;
    const t0 = performance.now();
    const tick = (now) => {
      if (!S || !S.illuPlaying) return;
      const f = Math.min(1, (now - t0) / 1000 / ILLU_SECONDS);
      // advance in discrete denoising steps, like the sampler does
      illuFrame(Math.floor(f * ILLU_STEPS) / ILLU_STEPS);
      if (f < 1) requestAnimationFrame(tick); else S.illuPlaying = false;
    };
    requestAnimationFrame(tick);
  }

  function showIllustration() {
    S.view = "illu";
    S.illuPlaying = false;
    buildIllustration();
    $("#stage-overlay").innerHTML = "";
    $("#stage-overlay").appendChild(el(`<div class="legend-card"><b>Illustration, not a prediction.</b> The crystal structure with noise added and then removed step by step: the idea behind the diffusion module of AlphaFold 3, Boltz and Chai. A real model starts from noise alone and has to <i>predict</i> each step.<div id="cf-step" style="margin-top:4px"></div></div>`));
    illuFrame(1);
  }

  // ------------------------------------------------------------------ comparison
  async function readUpload(file) {
    const name = file.name.toLowerCase();
    if (name.endsWith(".zip")) {
      await PW.loadScript(FFLATE_URL);
      const files = window.fflate.unzipSync(new Uint8Array(await file.arrayBuffer()));
      const names = Object.keys(files).filter((n) => /\.(pdb|cif)$/i.test(n) && !/template/i.test(n));
      const pick = names.find((n) => /model_0\.cif$/i.test(n)) || names.find((n) => /\.pdb$/i.test(n)) || names[0];
      if (!pick) throw new Error("no .pdb or .cif structure inside the zip");
      return { text: new TextDecoder().decode(files[pick]), name: `${file.name} › ${pick}` };
    }
    return { text: await file.text(), name: file.name };
  }

  async function handleFile(file) {
    const out = $("#cf-result");
    out.innerHTML = `<p class="small muted">Reading ${esc(file.name)}…</p>`;
    try {
      const { text, name } = await readUpload(file);
      const pred = Compare.parseStructure(text);
      const ref = { ca_by_seq: S.cof.ca_by_seq, pocket_seq: S.cof.pocket_seq, ligand: S.prot.ligand };
      const res = Compare.compareToCrystal(pred, ref);
      showComparison(pred, res);
      const okLig = res.ligandRmsd !== null && res.ligandRmsd <= S.successA;
      if (okLig) { PW.progress.tasks[`cofold:${S.cof.id}`] = true; PW.save(); }
      out.innerHTML = `<div class="box ${okLig ? "gold" : ""}">
        <div class="small muted">${esc(name)}</div>
        <div class="row" style="justify-content:space-between;margin-top:4px"><span>Ligand RMSD (pocket-aligned)</span>
          <span class="big" style="font-size:24px">${res.ligandRmsd === null ? "–" : `${res.ligandRmsd.toFixed(2)} Å`}</span></div>
        ${res.ligandRmsd !== null ? `<div class="small ${okLig ? "yes" : "no"}">${okLig ? `Within ${S.successA} Å of the crystal pose: a success by the usual docking criterion.` : `More than ${S.successA} Å from the crystal pose.`}</div>` : ""}
        <div class="readout" style="margin-top:8px"><span>Protein CA RMSD (whole chain)</span><b>${res.caRmsd.toFixed(2)} Å</b>
          <span>Pocket CA RMSD</span><b>${res.pocketCaRmsd.toFixed(2)} Å</b>
          <span>Residues matched</span><b>${res.nResidues}</b><span>Pocket residues used</span><b>${res.nPocket}</b></div>
        ${res.warnings.map((w) => `<p class="small" style="color:var(--warn)">${esc(w)}</p>`).join("")}
        <p class="small muted">Superposed on the CA atoms of the ${res.nPocket} pocket residues, then the ligand's heavy atoms compared by name, taking the best of ${S.prot.ligand.automorphisms.length} symmetry-equivalent matchings.</p></div>`;
    } catch (e) {
      console.error(e);
      out.innerHTML = `<div class="box" style="border-color:var(--bad)"><b>Could not compare.</b> <span class="small">${esc(e.message)}</span></div>`;
    }
  }

  function showComparison(pred, res) {
    S.view = "compare";
    S.illuPlaying = false;
    const v = PW.resetView();
    S.viewer = v;
    const rec = v.addModel(S.prot.receptor_pdb, "pdb");
    rec.setStyle({}, { cartoon: { color: "#9aa3b2", opacity: 0.7 } });
    const lig = v.addModel(S.prot.ligand.molblock, "sdf");
    lig.setStyle({}, { stick: { radius: 0.25, colorscheme: "greenCarbon" } });
    // predicted CA trace and ligand, moved into the crystal frame
    const T = res.transform, lines = [];
    let n = 1;
    [...pred.protein.entries()].sort((a, b) => a[0] - b[0]).forEach(([k, xyz]) => {
      const [x, y, z] = Compare.apply(T.R, T.t, xyz);
      lines.push(`ATOM  ${String(n++).padStart(5)}  CA  ${THREE[S.cof.sequence[k - 1]] || "UNK"} P${String(k).padStart(4)}    ${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}  1.00  0.00           C`);
    });
    pred.ligand.forEach((a) => {
      const [x, y, z] = Compare.apply(T.R, T.t, a.xyz);
      lines.push(`HETATM${String(n++).padStart(5)} ${a.name.padEnd(4).slice(0, 4)} LIG Q   1    ${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}  1.00  0.00          ${a.elem.padStart(2)}`);
    });
    const pm = v.addModel(lines.join("\n") + "\nEND\n", "pdb");
    pm.setStyle({ chain: "P" }, { cartoon: { style: "trace", thickness: 0.3, color: "#60a5fa" } });
    pm.setStyle({ chain: "Q" }, { stick: { radius: 0.22, colorscheme: "magentaCarbon" } });
    v.zoomTo({ model: lig }, 0);
    v.zoom(0.5);
    v.render();
    $("#stage-overlay").innerHTML = "";
    $("#stage-overlay").appendChild(el(`<div class="legend-card"><span style="color:#c7cfdc">Grey</span> crystal protein, <span style="color:#7ee787">green</span> crystal ligand.
      <span style="color:#60a5fa">Blue</span> your predicted backbone, <span style="color:#ff8fd0">pink</span> your predicted ligand, superposed on the pocket.</div>`));
  }

  // ------------------------------------------------------------------ panels
  function renderPanels() {
    const c = S.cof, L = S.index.levels.find((x) => x.id === c.id);
    $("#left").innerHTML = `<h3>Co-folding</h3>
      <p class="small">Docking needs a protein structure and moves a ligand into it. A co-folding model starts from the <b>sequence</b> and the <b>ligand</b> and builds both together (slides 39 and 57).</p>
      <div class="levels" id="cf-list"></div>
      <h3>1 · The idea</h3>
      <p class="small">AlphaFold 3, Boltz and Chai generate the structure by <i>denoising</i>: start from random atom positions and remove noise in many small steps, each one predicted by the network.</p>
      <div class="row"><button class="primary" id="cf-play">▶ Watch the idea</button></div>
      <h3>Why this is hard</h3>
      <p class="small">The network has to get the fold, the pocket and the ligand pose right at once. It has seen many complexes in training; how well it does on pockets and ligands unlike those is an open question (slide 57).</p>`;
    S.index.cofold.forEach((t) => {
      const done = PW.progress.tasks[`cofold:${t.id}`];
      const b = el(`<button class="lvl ${t.id === c.id ? "on" : ""}"><span class="n">${done ? "✓" : "◇"}</span>
        <span class="t">${esc(t.ligand_label)} + ${esc(S.index.levels.find((x) => x.id === t.id).title.split(" + ")[0])}<small>${t.length} residues · PDB ${t.pdb_id}</small></span><span></span></button>`);
      b.addEventListener("click", () => openTarget(t.id));
      $("#cf-list").appendChild(b);
    });
    $("#cf-play").addEventListener("click", playIllustration);

    $("#right").innerHTML = `<h3>2 · Fold it for real</h3>
      <p class="small"><a href="${LOCALFOLD_URL}" target="_blank" rel="noopener">LocalFold</a> runs AlphaFold 3 and other co-folding models <b>inside your browser</b> (it needs WebGPU: a recent Chrome, Edge or Safari) and lets you watch the structure being built.</p>
      <ol class="steps small">
        <li><a href="${c.job_file}" download="${esc(c.id)}_job.json">Download the job file</a> for ${esc(L.ligand_label)} + ${esc(L.title.split(" + ")[0])}.</li>
        <li>Open <a href="${LOCALFOLD_URL}" target="_blank" rel="noopener">localfold.org</a>. Under alignments choose <b>Load job JSON…</b> and pick the file.</li>
        <li>Keep a model that can place ligands (the default AF3, or e.g. Boltz-2 or Protenix-v2; AlphaFold 2 and ESMFold cannot). Press <b>Fold</b>.</li>
        <li>Watch the trajectory, then <b>Download All</b>.</li>
        <li>Drop the zip (or the structure file inside it) below.</li></ol>
      <p class="small muted">Or type it in yourself:</p>
      <div class="small">Protein, ${c.sequence.length} residues (the deposited sequence of PDB ${c.pdb_id} chain ${c.chain}) <button class="ghost" data-copy="seq">Copy</button></div>
      <div class="seqbox">${esc(c.sequence)}</div>
      <div class="small" style="margin-top:6px">Ligand: CCD code <code>${esc(c.ligand_id)}</code> <button class="ghost" data-copy="ccd">Copy</button></div>
      <div class="small muted" style="word-break:break-all">SMILES <code>${esc(c.ligand_smiles)}</code> <button class="ghost" data-copy="smi">Copy</button></div>
      <p class="small muted">Use the CCD code rather than the SMILES: then the predicted ligand keeps the crystal's atom names and can be scored below.</p>
      <h3>3 · Check your prediction</h3>
      <label class="drop" id="cf-drop"><input type="file" id="cf-file" accept=".zip,.pdb,.cif" hidden>
        Drop LocalFold's zip, or a <code>.pdb</code> / <code>.cif</code> from AlphaFold Server, here<br><span class="small">or click to choose a file</span></label>
      <div id="cf-result"></div>
      <p class="small muted" style="margin-top:10px">Crystal structure: ${citation(L.citation)}. ${esc(c.note)}</p>`;
    const copies = { seq: c.sequence, ccd: c.ligand_id, smi: c.ligand_smiles };
    $("#right").querySelectorAll("[data-copy]").forEach((b) => b.addEventListener("click", async () => {
      try { await navigator.clipboard.writeText(copies[b.dataset.copy]); toast("Copied.", "good", 1500); }
      catch (_) { toast("Copy failed: select the text by hand.", "warn"); }
    }));
    const drop = $("#cf-drop"), input = $("#cf-file");
    input.addEventListener("change", () => { if (input.files[0]) handleFile(input.files[0]); });
    drop.addEventListener("dragover", (e) => { e.preventDefault(); drop.classList.add("over"); });
    drop.addEventListener("dragleave", () => drop.classList.remove("over"));
    drop.addEventListener("drop", (e) => { e.preventDefault(); drop.classList.remove("over"); if (e.dataTransfer.files[0]) handleFile(e.dataTransfer.files[0]); });
  }

  async function openTarget(id) {
    PW.setArg(id);
    S.cof = await PW.get("cofold", id);
    S.prot = await PW.get("prot", id);
    const L = S.index.levels.find((x) => x.id === id);
    PW.setHud("Co-folding", `${L.ligand_label} + ${L.title.split(" + ")[0]}`);
    renderPanels();
    showIllustration();
  }

  PW.chapter({
    id: "cofold", title: "Co-folding", color: "#60a5fa", num: "Chapter 5",
    blurb: "Give a model only the sequence and the ligand, and watch it build the complex. Then score its prediction against the crystal.",
    learn: "AlphaFold 3 · diffusion · in-browser folding with LocalFold · pocket-aligned RMSD",
    progress(index) {
      const done = index.cofold.filter((t) => PW.progress.tasks[`cofold:${t.id}`]).length;
      return { text: `${done}/${index.cofold.length}`, done: done === index.cofold.length };
    },
    async enter(arg) {
      const index = await PW.get("index", "index");
      S = { index, successA: index.docking.success_rmsd_A };
      await openTarget(index.cofold.some((t) => t.id === arg) ? arg : index.cofold[0].id);
    },
    leave() { if (S) S.illuPlaying = false; S = null; },
    recentre() {
      if (!S) return;
      if (S.view === "illu") S.viewer.zoomTo(); else S.viewer.zoomTo();
      S.viewer.render();
    },
  });
})();
