// Chapter 5: co-folding. (1) AlphaFold 3's diffusion sampler, step by step (js/diffusion.js), (2) a hand-off to LocalFold,
// which runs AlphaFold 3 and similar models in the browser, (3) scoring the student's prediction
// against the crystal structure (js/compare.js, tested in tests/test_compare.js).
(function () {
  "use strict";
  const { $, esc, el, citation, toast, vec } = PW;
  const FFLATE_URL = "https://cdn.jsdelivr.net/npm/fflate@0.8.2/umd/index.js";
  const LOCALFOLD_URL = "https://localfold.org";
  const THREE = { A: "ALA", R: "ARG", N: "ASN", D: "ASP", C: "CYS", Q: "GLN", E: "GLU", G: "GLY", H: "HIS", I: "ILE", L: "LEU", K: "LYS", M: "MET", F: "PHE", P: "PRO", S: "SER", T: "THR", W: "TRP", Y: "TYR", V: "VAL" };

  let S = null;

  // ------------------------------------------------------------------ AlphaFold 3's sampler, watched
  // The trajectory comes from js/diffusion.js (AF3's sample(), with an oracle denoiser). Here it is
  // only drawn: three views of the same atoms, a step slider, and the numbers of the current step.
  const PLAY_SECONDS = 12;
  const CHAIN_FROM = [86, 132, 230], CHAIN_TO = [196, 140, 255]; // N terminus blue -> C terminus violet

  // The folded chain's protein atoms only: streptavidin's data hold the whole tetramer, but one
  // chain is folded here (as in the job file).
  const chainPdb = () => S.prot.receptor_pdb.split("\n").filter((l) => l.startsWith("ATOM") && l[21] === S.cof.chain).join("\n") + "\nEND\n";

  function buildSampler(seed) {
    S.illuPlaying = false;
    const v = PW.resetView();
    S.viewer = v; S.view = "illu";
    // the crystal's protein heavy atoms (bonds assigned by 3Dmol from the crystal geometry) and the
    // ligand with its bonds from the CCD: the connectivity is input to the model, not something it finds
    const protPdb = chainPdb();
    S.allModel = v.addModel(protPdb, "pdb");
    S.ligModel = v.addModel(S.prot.ligand.molblock, "sdf");
    S.protAtoms = S.allModel.selectedAtoms({});
    S.ligAtoms = S.ligModel.selectedAtoms({});
    const resis = S.protAtoms.map((a) => a.resi), lo = Math.min(...resis), hi = Math.max(...resis);
    S.tone = (resi) => {
      const f = (resi - lo) / Math.max(1, hi - lo);
      return `rgb(${CHAIN_FROM.map((c, k) => Math.round(c + (CHAIN_TO[k] - c) * f)).join(",")})`;
    };
    // C-alpha trace: one pseudo-atom per residue, bonded to the next residue when both are modelled
    S.caIdx = [];
    S.protAtoms.forEach((a, i) => { if (a.atom === "CA") S.caIdx.push(i); });
    S.caModel = v.addModel();
    S.caModel.addAtoms(S.caIdx.map((i, j) => {
      const a = S.protAtoms[i], nxt = S.caIdx[j + 1] !== undefined ? S.protAtoms[S.caIdx[j + 1]] : null;
      return { elem: "C", atom: "CA", resn: a.resn, resi: a.resi, chain: "T", serial: j, x: a.x, y: a.y, z: a.z,
        bonds: nxt && nxt.resi === a.resi + 1 ? [j + 1] : [], bondOrder: nxt && nxt.resi === a.resi + 1 ? [1] : [] };
    }));
    S.caAtoms = S.caModel.selectedAtoms({});
    const ref = [];
    S.protAtoms.forEach((a) => ref.push(a.x, a.y, a.z));
    S.ligAtoms.forEach((a) => ref.push(a.x, a.y, a.z));
    S.ref = Float64Array.from(ref);
    S.seed = seed;
    S.run = Diffusion.sample(S.ref, Compare, { seed });
    S.step = Math.min(S.step === undefined ? 0 : S.step, Diffusion.AF3.STEPS);
    renderSamplerOverlay();
    drawStep(S.step, true);
  }

  function setCoords(X) {
    const np = S.protAtoms.length;
    S.protAtoms.forEach((a, i) => { a.x = X[3 * i]; a.y = X[3 * i + 1]; a.z = X[3 * i + 2]; });
    S.ligAtoms.forEach((a, i) => { const k = np + i; a.x = X[3 * k]; a.y = X[3 * k + 1]; a.z = X[3 * k + 2]; });
    S.caAtoms.forEach((a, j) => { const i = S.caIdx[j]; a.x = X[3 * i]; a.y = X[3 * i + 1]; a.z = X[3 * i + 2]; });
  }

  function styleSampler() {
    const tone = (a) => S.tone(a.resi);
    // drawn size grows with the cloud, or atoms in a cloud thousands of A wide would be invisible dots
    const k = Math.max(1, Diffusion.spread(S.run.frames[S.step]) / 25);
    S.allModel.setStyle({}, {}); S.caModel.setStyle({}, {});
    if (S.simView === "trace") S.caModel.setStyle({}, { sphere: { radius: 0.55 * k, colorfunc: tone }, stick: { radius: 0.22 * k, colorfunc: tone } });
    else if (S.simView === "atoms") S.allModel.setStyle({}, { stick: { radius: 0.14 * k, colorfunc: tone }, sphere: { radius: 0.3 * k, colorfunc: tone } });
    else S.allModel.setStyle({}, { cartoon: { colorfunc: tone } });
    S.ligModel.setStyle({}, { stick: { radius: 0.2 * k, colorscheme: "greenCarbon" }, sphere: { radius: 0.42 * k, colorscheme: "greenCarbon" } });
    if (S.ghostModel) { S.viewer.removeModel(S.ghostModel); S.ghostModel = null; }
    if (S.showGhost) { // the oracle's answer, in the display frame: the crystal structure
      S.ghostModel = S.viewer.addModel(chainPdb(), "pdb");
      S.ghostModel.setStyle({}, { cartoon: { style: "trace", color: "#ffffff", opacity: 0.35, thickness: 0.25 } });
    }
  }

  function drawStep(step, refit) {
    S.step = step;
    const X = S.run.frames[step], f = S.run.info[step];
    setCoords(X);
    styleSampler();
    // AF3 has no fixed frame; each step is superposed onto the crystal and the camera fits the cloud
    if (refit !== false) S.viewer.zoomTo({ model: S.simView === "trace" ? S.caModel : S.allModel });
    S.viewer.render();
    const q = (id) => document.getElementById(id);
    if (S.slider) S.slider.set(step);
    const fmt = (x) => (x >= 100 ? x.toFixed(0) : x >= 1 ? x.toFixed(1) : x >= 0.01 ? x.toFixed(3) : x.toExponential(1));
    if (q("cf-r-step")) {
      q("cf-r-step").textContent = `${step} / ${Diffusion.AF3.STEPS}`;
      q("cf-r-sigma").textContent = `${fmt(f.sigma)} Å`;
      q("cf-r-gamma").textContent = step === 0 ? "–" : f.gamma.toFixed(1);
      q("cf-r-that").textContent = step === 0 ? "–" : `${fmt(f.tHat)} Å`;
      q("cf-r-size").textContent = `${fmt(2 * Diffusion.spread(X))} Å`;
    }
    if (q("cf-plot-dot")) {
      const [px, py] = plotXY(step, f.sigma);
      q("cf-plot-dot").setAttribute("cx", px); q("cf-plot-dot").setAttribute("cy", py);
    }
    if (q("cf-ov-step")) q("cf-ov-step").innerHTML = `step <b>${step}</b> / ${Diffusion.AF3.STEPS} · σ = <b>${fmt(f.sigma)} Å</b>`;
    drawScaleBar();
  }

  // a scale bar in the overlay, so the shrinking cloud keeps its size in Angstrom
  function drawScaleBar() {
    const bar = document.getElementById("cf-scale");
    if (!bar) return;
    const v = S.viewer, view = v.getView(), [x, y, z, w] = [-view[4], -view[5], -view[6], view[7]];
    // screen x axis in model coordinates: the inverse view rotation applied to (1, 0, 0)
    const right = [1 - 2 * (y * y + z * z), 2 * (x * y + w * z), 2 * (x * z - w * y)];
    const X = S.run.frames[S.step], n = X.length / 3, c = [0, 0, 0];
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) c[k] += X[3 * i + k] / n;
    const [a, b] = v.modelToScreen([{ x: c[0], y: c[1], z: c[2] }, { x: c[0] + right[0], y: c[1] + right[1], z: c[2] + right[2] }]);
    const ppa = Math.hypot(b.x - a.x, b.y - a.y);
    if (!(ppa > 0)) return;
    const target = 120 / ppa, p10 = 10 ** Math.floor(Math.log10(target));
    const len = [1, 2, 5, 10].map((m) => m * p10).reduce((best, L) => (Math.abs(L * ppa - 120) < Math.abs(best * ppa - 120) ? L : best));
    bar.querySelector("i").style.width = `${Math.round(len * ppa)}px`;
    bar.querySelector("span").textContent = `${len >= 1 ? len.toLocaleString("en") : len} Å`;
  }

  // sigma against step, log scale, for the panel
  const PLOT = { w: 270, h: 104, l: 34, r: 8, t: 18, b: 18 };
  function plotXY(step, sigma) {
    const L = Diffusion.noiseLevels(), lmax = Math.log10(L[0]), lmin = Math.log10(L[L.length - 1]);
    const x = PLOT.l + (step / (L.length - 1)) * (PLOT.w - PLOT.l - PLOT.r);
    const y = PLOT.t + ((lmax - Math.log10(sigma)) / (lmax - lmin)) * (PLOT.h - PLOT.t - PLOT.b);
    return [x.toFixed(1), y.toFixed(1)];
  }
  function schedulePlot() {
    const L = Diffusion.noiseLevels(), A = Diffusion.AF3;
    const path = L.map((s, i) => `${i ? "L" : "M"}${plotXY(i, s).join(" ")}`).join("");
    const gStep = L.findIndex((s) => s <= A.GAMMA_MIN), [gx] = plotXY(gStep, L[gStep]);
    const ticks = [1000, 10, 0.1, 0.01].filter((t) => t <= L[0] && t >= L[L.length - 1])
      .map((t) => { const [, y] = plotXY(0, t); return `<text x="${PLOT.l - 4}" y="${+y + 3}" text-anchor="end">${t >= 1 ? t : t}</text><path d="M${PLOT.l} ${y}H${PLOT.w - PLOT.r}" stroke="#24304a"/>`; }).join("");
    return `<svg viewBox="0 0 ${PLOT.w} ${PLOT.h}" width="100%" font-size="9" fill="#9aa7bf" font-family="ui-monospace,monospace">
      ${ticks}<path d="M${gx} ${PLOT.t}V${PLOT.h - PLOT.b}" stroke="#ffb454" stroke-dasharray="3 3"/>
      <text x="${+gx - 3}" y="${PLOT.t - 6}" fill="#ffb454" text-anchor="end">γ = 0 once σ ≤ ${A.GAMMA_MIN} Å →</text>
      <path d="${path}" fill="none" stroke="#60a5fa" stroke-width="1.8"/>
      <circle id="cf-plot-dot" r="4" fill="#5eead4" stroke="#06131c"/>
      <text x="${PLOT.l}" y="${PLOT.h - 4}">step 0</text><text x="${PLOT.w - PLOT.r}" y="${PLOT.h - 4}" text-anchor="end">${A.STEPS}</text>
      <text x="4" y="10" fill="#cdd6e7">σ (Å)</text></svg>`;
  }

  function renderSamplerOverlay() {
    $("#stage-overlay").innerHTML = "";
    $("#stage-overlay").appendChild(el(`<div class="legend-card"><b>AlphaFold 3's sampler, with one stand-in.</b>
      The steps, noise levels and random rotations are AlphaFold 3's. The network that guesses the clean structure is
      replaced by the answer: the crystal structure. <a href="#" id="cf-truth">What else differs?</a>
      <div id="cf-ov-step" style="margin-top:4px"></div></div>`));
    $("#stage-overlay").appendChild(el(`<div class="legend-card" id="cf-scale" style="align-self:flex-end"><i style="display:block;height:3px;background:var(--acc);border-radius:2px;width:60px"></i><span></span></div>`));
    $("#cf-truth").addEventListener("click", (e) => { e.preventDefault(); const d = $("#cf-truth-box"); if (d) { d.open = true; d.scrollIntoView({ block: "nearest" }); } });
  }

  function playSampler() {
    if (S.view !== "illu") buildSampler(S.seed || 1);
    if (S.illuPlaying) { S.illuPlaying = false; updatePlayButton(); return; }
    if (S.step >= Diffusion.AF3.STEPS) S.step = 0;
    S.illuPlaying = true; updatePlayButton();
    const start = S.step, t0 = performance.now();
    const tick = (now) => {
      if (!S || !S.illuPlaying) return;
      const step = Math.min(Diffusion.AF3.STEPS, start + Math.floor((Math.max(0, now - t0) / 1000 / PLAY_SECONDS) * Diffusion.AF3.STEPS)); // rAF time can predate t0
      if (step !== S.step) drawStep(step);
      if (step < Diffusion.AF3.STEPS) requestAnimationFrame(tick);
      else { S.illuPlaying = false; updatePlayButton(); toast("Done: switch the view to <b>Cartoon</b> to see the fold.", "good", 4000); }
    };
    requestAnimationFrame(tick);
  }
  function updatePlayButton() { const b = document.getElementById("cf-play"); if (b) b.textContent = S.illuPlaying ? "❚❚ Pause" : "▶ Play"; }

  function samplerPanel() {
    const A = Diffusion.AF3, wrap = $("#cf-sampler");
    wrap.innerHTML = `<div class="row"><button class="primary" id="cf-play">▶ Play</button><button id="cf-seed" class="ghost" title="Start again from different random noise">↻ New noise</button></div>
      <div id="cf-step-slider"></div>
      <div class="readout"><span>Step</span><b id="cf-r-step"></b><span>Noise level σ</span><b id="cf-r-sigma"></b>
        <span>Extra noise γ</span><b id="cf-r-gamma"></b><span>Raised level t̂ = σ<sub>prev</sub>(1+γ)</span><b id="cf-r-that"></b>
        <span>Cloud size (2 × RMS radius)</span><b id="cf-r-size"></b></div>
      <div style="margin-top:6px">${schedulePlot()}</div>
      <h3>View</h3><div id="cf-view"></div><div id="cf-ghost"></div>
      <details class="small" style="margin-top:8px"><summary><b>The algorithm</b> (sample() in AlphaFold 3's code)</summary>
<pre class="algo">x ← σ₀ · ε                    σ₀ = ${(A.SIGMA_DATA * A.S_MAX).toFixed(0)} Å, ε ~ N(0, I)
for τ = 1 … ${A.STEPS}:
  x ← rotate(x − mean) + shift   random, every step
  γ ← ${A.GAMMA_0} if σ<sub>τ</sub> > ${A.GAMMA_MIN} Å else 0
  t̂ ← σ<sub>τ−1</sub> · (1 + γ)
  x' ← x + ${A.NOISE_SCALE} · √(t̂² − σ<sub>τ−1</sub>²) · ε
  x̂ ← <b>network</b>(x', t̂, trunk)   the guess
  x ← x' + ${A.STEP_SCALE} · (σ<sub>τ</sub> − t̂) · (x' − x̂) / t̂</pre>
        <p>σ falls from ${(A.SIGMA_DATA * A.S_MAX).toFixed(0)} Å to ${(A.SIGMA_DATA * A.S_MIN).toFixed(4)} Å: σ<sub>data</sub> = ${A.SIGMA_DATA} Å times (${A.S_MAX}<sup>1/${A.P}</sup> + t(${A.S_MIN}<sup>1/${A.P}</sup> − ${A.S_MAX}<sup>1/${A.P}</sup>))<sup>${A.P}</sup>, so most steps are spent at small noise, where the details are decided.
        The last line moves from the noisy atoms towards the guess, ${A.STEP_SCALE} times as far as a plain Euler step would.</p></details>
      <details class="small" id="cf-truth-box"><summary><b>What is AlphaFold 3 here, and what is not</b></summary>
        <p><b>As in AlphaFold 3</b> (google-deepmind/alphafold3, diffusion_head.py): the noise schedule, ${A.STEPS} steps, γ<sub>0</sub> = ${A.GAMMA_0} above ${A.GAMMA_MIN} Å, noise scale ${A.NOISE_SCALE}, step scale ${A.STEP_SCALE}, the random rotation and 1 Å shift at every step, and the start from pure noise centred on the origin.</p>
        <p><b>Different:</b></p><ul>
          <li>The <b>guess x̂</b> is the crystal structure superposed onto the noisy atoms. AlphaFold 3's network has to predict it from the sequence, MSA and templates (via the trunk), and early guesses are blurry. This is why the result here is always right.</li>
          <li><b>One sample</b>. AlphaFold 3 runs ${A.NUM_SAMPLES} per seed and ranks them by confidence.</li>
          <li><b>Atoms</b>: the protein heavy atoms modelled in the crystal plus the ligand. AlphaFold 3 also builds residues missing from the crystal (which shows here as gaps in the trace), and ions if you give them.</li>
          <li><b>Display</b>: every step is superposed onto the crystal and the camera fits the cloud (see the scale bar). AlphaFold 3 has no fixed frame; the random rotations would make the cloud tumble.</li>
          <li><b>Speed</b>: ${PLAY_SECONDS} s for all steps, evenly; the real model needs a network evaluation per step.</li></ul></details>`;
    S.slider = PW.slider({ label: "Step", min: 0, max: A.STEPS, step: 1, value: S.step || 0,
      onInput: (v) => { S.illuPlaying = false; updatePlayButton(); if (S.view !== "illu") buildSampler(S.seed || 1); drawStep(v); } });
    $("#cf-step-slider").appendChild(S.slider);
    $("#cf-view").appendChild(PW.segmented([
      { value: "trace", label: "Cα trace", title: "One sphere per residue, joined along the chain" },
      { value: "atoms", label: "All atoms", title: "Every heavy atom that is diffused, with its bonds" },
      { value: "cartoon", label: "Cartoon", title: "The fold (best at the end)" }], S.simView, (m) => {
      S.simView = m; if (S.view !== "illu") buildSampler(S.seed || 1); else drawStep(S.step); }));
    $("#cf-ghost").appendChild(PW.checkbox("Crystal structure (the stand-in's answer)", S.showGhost, (on) => { S.showGhost = on; if (S.view === "illu") drawStep(S.step, false); }));
    $("#cf-play").addEventListener("click", playSampler);
    $("#cf-seed").addEventListener("click", () => { S.step = 0; buildSampler((S.seed || 1) + 1); playSampler(); });
  }

  // ------------------------------------------------------------------ intro
  function intro(force) {
    if (!force && PW.seen("intro:cofold")) return Promise.resolve();
    PW.markSeen("intro:cofold");
    const A = Diffusion.AF3, n = 5;
    const k = (i) => `Co-folding · ${i}/${n}`;
    return PW.cards([
      { kicker: k(1), title: "Docking or co-folding?",
        html: `Docking (Chapter 3) takes a <b>finished</b> protein structure and searches for the ligand pose. A co-folding model (AlphaFold 3, Boltz, Chai, Protenix) gets only the <b>sequence</b> and the <b>ligand</b> (a CCD code or SMILES) and builds both at once, so the pocket can form around the ligand.` },
      { kicker: k(2), title: "The trunk: what the network knows",
        html: `Every residue and every ligand atom becomes a token. The trunk reads the MSA and templates and fills a <b>pair representation</b>: for every two tokens, what the network believes about how they sit relative to each other. Think of it as a soft distance map. It does not yet contain coordinates.` },
      { kicker: k(3), title: "The diffusion module: noise to atoms",
        html: `Coordinates are made by denoising. Every atom starts at a random position: noise of σ = ${(A.SIGMA_DATA * A.S_MAX).toFixed(0)} Å around the origin, a cloud thousands of Å wide. In ${A.STEPS} steps the network guesses the clean structure from the noisy one, guided by the trunk, and the sampler moves part of the way towards the guess.<div class="facts">AlphaFold 2 worked differently: its structure module starts every residue at the origin (a "black hole") and moves them out.</div>` },
      { kicker: k(4), title: "Several answers, ranked",
        html: `Different noise gives different structures. AlphaFold 3 draws ${A.NUM_SAMPLES} samples per seed and ranks them by its own confidence (pLDDT, PAE, ipTM). <b>Confident is not the same as correct</b>, especially for ligands and pockets unlike those it saw in training.` },
      { kicker: k(5), title: "In this chapter",
        html: `<b>1</b> Watch AlphaFold 3's sampler step by step (left; one part is replaced, and the page says which). <b>2</b> Fold the complex for real in your browser with LocalFold (right). <b>3</b> Upload the prediction and score its ligand pose against the crystal.` },
    ], { done: "Let's go" });
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
    out.scrollIntoView({ behavior: "smooth", block: "nearest" });
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
    $("#left").innerHTML = `<button id="cf-intro" class="ghost small" style="margin-bottom:6px">What is co-folding?</button>
      <h3>Co-folding</h3>
      <p class="small">Docking needs a protein structure and moves a ligand into it. A co-folding model starts from the <b>sequence</b> and the <b>ligand</b> and builds both together (slides 39 and 57).</p>
      <div class="levels" id="cf-list"></div>
      <h3>1 · Watch AlphaFold 3 sample</h3>
      <p class="small">AlphaFold 3, Boltz and Chai make coordinates by <i>denoising</i>: start from random atom positions and remove the noise in many small steps. Press play, or drag the step slider.</p>
      <div id="cf-sampler"></div>
      <h3>Why this is hard</h3>
      <p class="small">The network has to get the fold, the pocket and the ligand pose right at once. It has seen many complexes in training; how well it does on pockets and ligands unlike those is an open question (slide 57).</p>`;
    S.index.cofold.forEach((t) => {
      const done = PW.progress.tasks[`cofold:${t.id}`];
      const b = el(`<button class="lvl ${t.id === c.id ? "on" : ""}"><span class="n">${done ? "✓" : "◇"}</span>
        <span class="t">${esc(t.ligand_label)} + ${esc(S.index.levels.find((x) => x.id === t.id).title.split(" + ")[0])}<small>${t.length} residues · PDB ${t.pdb_id}</small></span><span></span></button>`);
      b.addEventListener("click", () => openTarget(t.id));
      $("#cf-list").appendChild(b);
    });
    $("#cf-intro").addEventListener("click", () => intro(true));
    samplerPanel();

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
    S.illuPlaying = false;
    S.step = 0;
    renderPanels();
    buildSampler(S.seed || 1);
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
      S = { index, successA: index.docking.success_rmsd_A, simView: "trace", showGhost: false, seed: 1 };
      await openTarget(index.cofold.some((t) => t.id === arg) ? arg : index.cofold[0].id);
    },
    afterTour() { return intro(false); },
    tour: () => [
      { sel: "#cf-intro", title: "What is co-folding?", text: "Five short cards on how AlphaFold 3 builds a complex. They open by themselves the first time." },
      { sel: "#cf-list", title: "Three complexes", text: "Each can be watched, folded for real, and scored." },
      { sel: "#cf-sampler", title: "AlphaFold 3's sampler", text: "Play the 200 denoising steps or drag the slider. The numbers and the plot show the noise level at each step; the views switch between Cα trace, all atoms and cartoon." },
      { sel: "#right", title: "Fold it for real", text: "LocalFold runs co-folding models in your browser. Upload the result below to score the ligand pose." },
    ],
    leave() { if (S) S.illuPlaying = false; S = null; },
    recentre() {
      if (!S) return;
      if (S.view === "illu") S.viewer.zoomTo(); else S.viewer.zoomTo();
      S.viewer.render();
    },
  });
})();
