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
  const progress = { stars: {}, tasks: {}, unlockAll: false, docking: {} };
  try { Object.assign(progress, JSON.parse(localStorage.getItem(SAVE_KEY) || "{}")); } catch (_) { /* storage blocked */ }
  function save() { try { localStorage.setItem(SAVE_KEY, JSON.stringify(progress)); } catch (_) { /* ignore */ } }
  function resetProgress() {
    progress.stars = {}; progress.tasks = {}; progress.unlockAll = false; progress.docking = {};
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
      $$("button[data-i]", box).forEach((b) => b.addEventListener("click", () => {
        m.classList.add("hidden");
        const btn = buttons[+b.dataset.i];
        resolve(btn.value !== undefined ? btn.value : +b.dataset.i);
      }));
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
      <input type="range" min="${o.min}" max="${o.max}" step="${o.step}" value="${o.value}">
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
    if (current && current.leave) current.leave();
    current = ch;
    $("#menu").classList.add("hidden");
    $("#game").classList.remove("hidden");
    $("#left").innerHTML = ""; $("#right").innerHTML = ""; $("#hud-mid").innerHTML = ""; $("#stage-overlay").innerHTML = "";
    setHud(ch.title, "");
    resetView();
    const hash = `#/${id}${arg !== undefined && arg !== null ? `/${arg}` : ""}`;
    if (location.hash !== hash) history.replaceState(null, "", hash);
    try {
      await ch.enter(arg);
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
    $("#btn-recentre").addEventListener("click", () => { if (current && current.recentre) current.recentre(); });
    window.addEventListener("hashchange", route);
    route();
  }

  return {
    data, get, loadScript, progress, save, resetProgress,
    $, $$, esc, el, stars, citation, toast, modal, segmented, slider, checkbox,
    view, resetView, vec, dots, chapter, chapters, open, setArg, home, setHud, boot,
    get current() { return current; },
  };
})();
