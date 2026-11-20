// Shared plumbing: data loading, progress, the 3D viewer, panels, toasts, chapter routing.
// Chapters register with PW.chapter({...}) and own the game screen while they are open.
window.PW = (function () {
  "use strict";

  // ------------------------------------------------------------------ data
  // Data files are <script>s that call PW.data(kind, key, obj); see scripts/build_data.py.
  const store = {}, loading = {};
  function data(kind, key, obj) { store[`${kind}:${key}`] = obj; }
  function loadScript(src) {
    if (!loading[src]) {
      loading[src] = new Promise((resolve, reject) => {
        const s = document.createElement("script");
        s.src = src;
        s.onload = () => resolve();
        s.onerror = () => reject(new Error(`could not load ${src}`));
        document.head.appendChild(s);
      });
    }
    return loading[src];
  }
  async function get(kind, key) {
    const file = kind === "index" ? "data/index.js" : `data/${kind}_${key}.js`;
    if (!store[`${kind}:${key}`]) await loadScript(file);
    return store[`${kind}:${key}`];
  }

  // ------------------------------------------------------------------ progress (this browser only)
  const SAVE_KEY = "pocketworks-v1";
  const progress = { stars: {}, tasks: {}, unlockAll: false, docking: {}, seen: {} };
  try { Object.assign(progress, JSON.parse(localStorage.getItem(SAVE_KEY) || "{}")); } catch (_) { /* storage blocked */ }
  function save() { try { localStorage.setItem(SAVE_KEY, JSON.stringify(progress)); } catch (_) { /* ignore */ } }
  function resetProgress() {
    progress.stars = {}; progress.tasks = {}; progress.unlockAll = false; progress.docking = {}; progress.seen = {};
    save();
  }

  // ------------------------------------------------------------------ DOM helpers
  const $ = (sel, root = document) => root.querySelector(sel);
  const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  function el(html) {
    const t = document.createElement("template");
    t.innerHTML = html.trim();
    return t.content.firstElementChild;
  }
  function stars(n, max = 3) {
    return `<span class="stars">${"★".repeat(n)}<span class="off">${"★".repeat(Math.max(0, max - n))}</span></span>`;
  }
  function citation(c) {
    if (!c || !c.authors) return "";
    const doi = c.doi ? ` · <a href="https://doi.org/${esc(c.doi)}" target="_blank" rel="noopener">paper</a>` : "";
    return `${esc(c.authors)}, ${esc(c.journal)} ${esc(c.year || "")}${doi}`;
  }

  function toast(html, kind = "", ms = 4500) {
    const t = el(`<div class="toast ${kind}">${html}</div>`);
    $("#toasts").appendChild(t);
    setTimeout(() => {
      t.style.transition = "opacity .4s";
      t.style.opacity = "0";
      setTimeout(() => t.remove(), 450);
    }, ms);
  }

  function modal(html, buttons = [{ label: "OK", primary: true }]) {
    return new Promise((resolve) => {
      const m = $("#modal"), box = $(".box", m);
      box.innerHTML = html + `<div class="row" style="justify-content:flex-end;margin-top:14px">` +
        buttons.map((b, i) => `<button data-i="${i}" class="${b.primary ? "primary" : ""}">${esc(b.label)}</button>`).join("") + `</div>`;
      m.classList.remove("hidden");
      const close = (i) => {
        m.classList.add("hidden");
        document.removeEventListener("keydown", esc_);
        const btn = buttons[i];
        resolve(btn.value !== undefined ? btn.value : i);
      };
      // Esc picks the first button (Cancel, or OK when it is the only one)
      const esc_ = (e) => { if (e.key === "Escape") { e.stopPropagation(); close(0); } };
      document.addEventListener("keydown", esc_);
      $$("button[data-i]", box).forEach((b) => b.addEventListener("click", () => close(+b.dataset.i)));
      $(".primary", box) ? $(".primary", box).focus() : $("button", box).focus();
    });
  }

  // Segmented control. options: [{value, label, title?}]; onPick(value) on change.
  function segmented(options, value, onPick) {
    const s = el(`<div class="seg" role="radiogroup">${options.map((o) =>
      `<button role="radio" data-v="${esc(o.value)}" aria-checked="${o.value === value}" title="${esc(o.title || "")}">${o.label}</button>`).join("")}</div>`);
    $$("button", s).forEach((b) => b.addEventListener("click", () => {
      $$("button", s).forEach((x) => x.setAttribute("aria-checked", String(x === b)));
      onPick(b.dataset.v);
    }));
    return s;
  }

  // Labelled range slider. o: {label, min, max, step, value, unit, why, onInput, onChange}
  function slider(o) {
    const fmt = o.format || ((v) => `${v}${o.unit || ""}`);
    const s = el(`<div class="slider"><label><span>${o.label}</span><b></b></label>
      <input type="range" aria-label="${esc(String(o.label).replace(/<[^>]+>/g, ""))}" min="${o.min}" max="${o.max}" step="${o.step}" value="${o.value}">
      ${o.why ? `<div class="why">${o.why}</div>` : ""}</div>`);
    const input = $("input", s), out = $("b", s);
    const show = () => { out.textContent = fmt(+input.value); };
    show();
    input.addEventListener("input", () => { show(); if (o.onInput) o.onInput(+input.value); });
    input.addEventListener("change", () => { if (o.onChange) o.onChange(+input.value); });
    s.set = (v) => { input.value = v; show(); };
    s.get = () => +input.value;
    return s;
  }

  function checkbox(label, checked, onChange) {
    const c = el(`<label class="check"><input type="checkbox" ${checked ? "checked" : ""}><span>${label}</span></label>`);
    $("input", c).addEventListener("change", (e) => onChange(e.target.checked));
    return c;
  }

  // ------------------------------------------------------------------ viewer
  let viewer = null;
  function view() {
    if (!viewer) {
      const bg = getComputedStyle(document.documentElement).getPropertyValue("--viewer-bg").trim() || "#0d1426";
      viewer = $3Dmol.createViewer($("#stage"), { backgroundColor: bg, antialias: true });
      window.addEventListener("resize", () => viewer.resize());
    }
    return viewer;
  }
  function resetView() {
    const v = view();
    v.clear();
    v.removeAllShapes();
    v.removeAllLabels();
    v.removeAllSurfaces();
    v.resize();
    v.render();
    return v;
  }
  const vec = (p) => ({ x: p[0], y: p[1], z: p[2] });
  // Many small spheres drawn as a model of pseudo-atoms: much faster than shape spheres, which
  // recompute their bounding sphere on every addition. Returns the model (remove with removeModel).
  function dots(v, points, color, radius = 0.3, opacity = 1) {
    const m = v.addModel();
    m.addAtoms(points.map((p, i) => ({ elem: "X", atom: "X", resn: "DOT", serial: i, x: p[0], y: p[1], z: p[2], hetflag: true, bonds: [] })));
    m.setStyle({}, { sphere: { radius, color, opacity } });
    return m;
  }

  // Molecular surface of a model's protein atoms. colorfunc(atom) -> colour. Returns a handle whose
  // restyle() and remove() wait for the (asynchronous) surface computation to finish.
  function surface(v, model, colorfunc, opacity = 1) {
    const p = v.addSurface($3Dmol.SurfaceType.MS, { colorfunc, opacity }, { model, hetflag: false }, { model, hetflag: false });
    let gone = false;
    return {
      restyle(cf, op) { p.then((id) => { if (!gone && v.getSurface(id)) { v.setSurfaceMaterialStyle(id, { colorfunc: cf, opacity: op }); v.render(); } }); },
      remove() { gone = true; p.then((id) => { if (v.getSurface(id)) { v.removeSurface(id); v.render(); } }); },
    };
  }

  // ------------------------------------------------------------------ "seen" flags (tours, intros)
  const seen = (key) => !!progress.seen[key];
  function markSeen(key) { progress.seen[key] = true; save(); }

  // ------------------------------------------------------------------ guided tour (coach marks)
  // steps: [{sel, title, text}]; sel may be omitted for a centred step. Steps whose element is
  // missing or hidden are skipped. Resolves when the tour ends (finished or skipped).
  function tour(steps) {
    steps = steps.filter((st) => !st.sel || visible($(st.sel)));
    if (!steps.length) return Promise.resolve();
    return new Promise((resolve) => {
      const root = el(`<div id="tour" role="dialog" aria-modal="true" aria-label="Guided tour">
        <div class="spot"></div>
        <div class="tip"><div class="count"></div><h4></h4><div class="txt"></div>
          <div class="row"><button class="ghost skip">Skip tour</button><span style="flex:1"></span>
          <button class="back">Back</button><button class="primary next">Next</button></div></div></div>`);
      document.body.appendChild(root);
      const spot = $(".spot", root), tip = $(".tip", root);
      const scrolled = ["#left", "#right", "#menu"].map((q) => [$(q), $(q) && $(q).scrollTop]);
      let i = 0;
      function place() {
        const st = steps[i], t = st.sel && $(st.sel);
        $(".count", tip).textContent = `${i + 1} / ${steps.length}`;
        $("h4", tip).innerHTML = st.title;
        $(".txt", tip).innerHTML = st.text;
        $(".back", tip).disabled = i === 0;
        $(".next", tip).textContent = i === steps.length - 1 ? "Got it" : "Next";
        const W = window.innerWidth, H = window.innerHeight, pad = 6, gap = 12;
        if (!t) {
          spot.style.cssText = `left:${W / 2}px;top:${H / 2}px;width:0;height:0`;
          tip.style.left = `${Math.max(12, (W - tip.offsetWidth) / 2)}px`;
          tip.style.top = `${Math.max(12, (H - tip.offsetHeight) / 2)}px`;
          return;
        }
        t.scrollIntoView({ block: "nearest", inline: "nearest" });
        const r = t.getBoundingClientRect();
        const x0 = Math.max(4, r.left - pad), y0 = Math.max(4, r.top - pad);
        const x1 = Math.min(W - 4, r.right + pad), y1 = Math.min(H - 4, r.bottom + pad);
        spot.style.cssText = `left:${x0}px;top:${y0}px;width:${x1 - x0}px;height:${y1 - y0}px`;
        const tw = tip.offsetWidth, th = tip.offsetHeight;
        // prefer below, then above, then right, then left of the target; else over its middle
        const opts = [
          [x0, y1 + gap, y1 + gap + th < H], [x0, y0 - gap - th, y0 - gap - th > 0],
          [x1 + gap, y0, x1 + gap + tw < W], [x0 - gap - tw, y0, x0 - gap - tw > 0],
          [(x0 + x1 - tw) / 2, (y0 + y1 - th) / 2, true]];
        const [lx, ly] = opts.find((o) => o[2]);
        tip.style.left = `${Math.min(Math.max(12, lx), W - tw - 12)}px`;
        tip.style.top = `${Math.min(Math.max(12, ly), H - th - 12)}px`;
      }
      function end() {
        window.removeEventListener("resize", place); document.removeEventListener("keydown", key, true);
        root.remove();
        scrolled.forEach(([e, top]) => { if (e) e.scrollTop = top; });
        resolve();
      }
      const go = (d) => { i += d; if (i >= steps.length) end(); else place(); };
      function key(e) {
        if (e.key === "Escape") { e.stopPropagation(); end(); }
        else if (e.key === "ArrowRight" || e.key === "Enter") { e.preventDefault(); e.stopPropagation(); go(1); }
        else if (e.key === "ArrowLeft" && i > 0) { e.preventDefault(); e.stopPropagation(); go(-1); }
      }
      $(".next", tip).addEventListener("click", () => go(1));
      $(".back", tip).addEventListener("click", () => go(-1));
      $(".skip", tip).addEventListener("click", end);
      window.addEventListener("resize", place);
      document.addEventListener("keydown", key, true);
      place();
      $(".next", tip).focus();
    });
  }
  function visible(e) { return !!e && e.getClientRects().length > 0 && getComputedStyle(e).visibility !== "hidden"; }

  // ------------------------------------------------------------------ intro cards over the 3D view
  // cards: [{kicker, title, html, art?}]. Short, skippable, slowly spinning protein behind them.
  // Resolves when closed. o.done: label of the last button.
  function cards(list, o = {}) {
    closeCards();
    return new Promise((resolve) => {
      const root = el(`<div class="intro-scrim"><div class="intro" role="dialog" aria-modal="false" aria-label="Introduction">
        <button class="x ghost" title="Skip (Esc)" aria-label="Skip">✕</button>
        <div class="art"></div><div class="kicker"></div><h3></h3><div class="txt"></div>
        <div class="row foot"><div class="dots">${list.map(() => "<i></i>").join("")}</div><span style="flex:1"></span>
          <button class="ghost skip">Skip</button><button class="back">Back</button><button class="primary next"></button></div></div></div>`);
      $("#stage-wrap").appendChild(root);
      // Slow turn of the protein behind the cards. requestAnimationFrame rather than viewer.spin(),
      // whose fixed 25 ms timer piles up on machines that need longer than that per frame.
      const v = viewer;
      let spinning = o.spin !== false && !!v, last = performance.now();
      const turn = (t) => {
        if (!spinning) return;
        v.rotate(Math.min(2, (t - last) * 0.015), "y");
        last = t;
        requestAnimationFrame(turn);
      };
      if (spinning) requestAnimationFrame(turn);
      let i = 0;
      function show() {
        const c = list[i];
        $(".art", root).innerHTML = c.art || ""; $(".art", root).classList.toggle("hidden", !c.art);
        $(".kicker", root).innerHTML = c.kicker || "";
        $("h3", root).innerHTML = c.title;
        $(".txt", root).innerHTML = c.html;
        $$(".dots i", root).forEach((d, j) => d.classList.toggle("on", j === i));
        $(".back", root).classList.toggle("hidden", i === 0);
        $(".skip", root).classList.toggle("hidden", i === list.length - 1);
        $(".next", root).textContent = i === list.length - 1 ? (o.done || "Start") : "Next →";
      }
      function end() {
        document.removeEventListener("keydown", key);
        spinning = false;
        root.remove(); cardsClose = null; resolve();
      }
      cardsClose = end;
      function key(e) {
        if (document.getElementById("tour") || !$("#modal").classList.contains("hidden")) return;
        if (e.key === "Escape") end();
        else if (e.key === "ArrowRight" || e.key === "Enter") { e.preventDefault(); if (++i >= list.length) end(); else show(); }
        else if (e.key === "ArrowLeft" && i > 0) { i--; show(); }
      }
      $(".next", root).addEventListener("click", () => { if (++i >= list.length) end(); else show(); });
      $(".back", root).addEventListener("click", () => { i--; show(); });
      $(".skip", root).addEventListener("click", end);
      $(".x", root).addEventListener("click", end);
      document.addEventListener("keydown", key);
      show();
    });
  }
  let cardsClose = null;
  function closeCards() { if (cardsClose) cardsClose(); }

  // ------------------------------------------------------------------ chapters and routing
  const chapters = [];
  let current = null;
  function chapter(c) { chapters.push(c); }

  function setHud(title, sub) {
    $("#hud-chapter").textContent = title;
    $("#hud-sub").textContent = sub || "";
  }

  async function open(id, arg) {
    const ch = chapters.find((c) => c.id === id);
    if (!ch) return home();
    closeCards();
    if (current && current.leave) current.leave();
    current = ch;
    $("#menu").classList.add("hidden");
    $("#game").classList.remove("hidden");
    $("#left").innerHTML = ""; $("#right").innerHTML = ""; $("#hud-mid").innerHTML = ""; $("#stage-overlay").innerHTML = "";
    $("#stage-hint").innerHTML = ""; $("#stage-hint").classList.add("hidden");
    setHud(ch.title, "");
    resetView();
    const hash = `#/${id}${arg !== undefined && arg !== null ? `/${arg}` : ""}`;
    if (location.hash !== hash) history.replaceState(null, "", hash);
    try {
      await ch.enter(arg);
      if (current === ch) await chapterTour(false);
      if (current === ch && ch.afterTour) await ch.afterTour();
    } catch (e) {
      console.error(e);
      toast(`Something went wrong: ${esc(e.message)}`, "bad", 8000);
    }
  }
  function setArg(arg) {
    if (!current) return;
    history.replaceState(null, "", `#/${current.id}/${arg}`);
  }
  function home() {
    closeCards();
    if (current && current.leave) current.leave();
    current = null;
    $("#game").classList.add("hidden");
    $("#menu").classList.remove("hidden");
    history.replaceState(null, "", "#/");
    if (window.PWMenu) window.PWMenu.render();
  }
  function route() {
    const parts = location.hash.replace(/^#\/?/, "").split("/").filter(Boolean);
    if (parts.length && chapters.some((c) => c.id === parts[0])) open(parts[0], parts[1]);
    else home();
  }

  // The game-screen tour: the general controls once (first chapter opened), then each chapter's own
  // steps the first time that chapter is opened. force: replay everything (the "Tour" button).
  const GAME_TOUR = [
    { title: "A quick tour", text: "Thirty seconds on the controls. <kbd>→</kbd> or <b>Next</b> to continue, <kbd>Esc</kbd> to skip. You can replay it any time with the <b>Tour</b> button at the top." },
    { sel: "#btn-home", title: "Back to the chapters", text: "Returns to the chapter menu. Your progress is kept." },
    { sel: "#hud-title", title: "Where you are", text: "The chapter and the current protein or level." },
    { sel: "#stage", title: "The 3D view", text: "<b>Drag</b> to rotate, <b>scroll</b> or pinch to zoom, <b>right-drag</b> (or <kbd>Ctrl</kbd>+drag) to move sideways." },
    { sel: "#btn-recentre", title: "Lost the protein?", text: "Re-centres the view." },
    { sel: "#btn-help", title: "Help", text: "How to move the view, if you forget." },
    { sel: "#left", title: "Controls", text: "The left panel holds the levels and settings of each chapter." },
    { sel: "#right", title: "Task and feedback", text: "The right panel says what to do, gives feedback, hints and, once you are done, the story behind the structure." },
  ];
  async function chapterTour(force) {
    if (!current) return;
    const ch = current, own = ch.tour ? ch.tour() : [];
    const general = force || !seen("tour:game");
    const mine = force || !seen(`tour:${ch.id}`);
    const steps = [...(general ? GAME_TOUR : []), ...(mine ? own : [])];
    if (!steps.length) return;
    markSeen("tour:game"); markSeen(`tour:${ch.id}`);
    await tour(steps);
  }
  const MENU_TOUR = [
    { title: "Welcome to Pocketworks", text: "Where does a drug bind, and how would a computer know? A 20-second tour of this page. <kbd>→</kbd> or <b>Next</b> to continue, <kbd>Esc</kbd> to skip." },
    { sel: "#chapter-cards .card", title: "Five chapters", text: "Each card is a chapter. Start with <b>Chapter 1</b>: the later ones build on it. Click a card to play." },
    { sel: "#chapter-cards .card .prog", title: "Your progress", text: "Levels or challenges done in each chapter. Progress is stored in this browser only." },
    { sel: "#btn-unlock", title: "Locked levels", text: "Levels unlock one after another. This button opens them all, e.g. if your teacher wants to jump ahead." },
    { sel: "#btn-reset", title: "Start over", text: "Deletes your stars and completed challenges." },
    { sel: "#btn-menu-tour", title: "Replay the tour", text: "Shows this tour again. Inside a chapter, the <b>Tour</b> button at the top explains the controls there." },
  ];
  function menuTour(force) {
    if (!force && seen("tour:menu")) return Promise.resolve();
    markSeen("tour:menu");
    return tour(MENU_TOUR);
  }

  function help() {
    modal(`<h2>Moving around</h2>
      <p><b>Rotate</b>: drag with the left mouse button (one finger on a touchscreen).</p>
      <p><b>Zoom</b>: scroll wheel, or pinch.</p>
      <p><b>Move the view sideways</b>: drag with the right mouse button, or with <kbd>Ctrl</kbd> held.</p>
      <p><b>Re-centre</b>: the target button at the top right.</p>
      <p class="muted small">Chapter-specific controls are explained in the left panel of each chapter.</p>`);
  }

  function boot() {
    $("#btn-home").addEventListener("click", home);
    $("#btn-help").addEventListener("click", help);
    $("#btn-tour").addEventListener("click", () => { closeCards(); chapterTour(true); });
    $("#btn-recentre").addEventListener("click", () => { if (current && current.recentre) current.recentre(); });
    window.addEventListener("hashchange", route);
    route();
  }

  return {
    data, get, loadScript, progress, save, resetProgress,
    $, $$, esc, el, stars, citation, toast, modal, segmented, slider, checkbox,
    view, resetView, vec, dots, surface, seen, markSeen, tour, menuTour, cards, closeCards,
    chapter, chapters, open, setArg, home, setHud, boot,
    get current() { return current; },
  };
})();
