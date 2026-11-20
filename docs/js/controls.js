// Mouse / touch / wheel / keyboard control of a ligand pose in a 3Dmol viewer, for Chapter 3.
// Copied from the ligand_docking_game repository (docs/js/controls.js, the controls agreed for that
// game); the only change is detach(), so the chapter can remove its listeners when it is left.
//
// The pose is x' = R (x - center) + center + t (see Metrics.applyPose). LigandControls.attach()
// mutates getPose().R / .t and calls refresh(); `pivot` is the model-space point the ligand is
// being rotated around (for drawing a marker), or null.
(function (root) {
  "use strict";

  // UI tuning (feel only, no scientific meaning).
  const KEY_STEP_A = 0.25, KEY_STEP_FINE_A = 0.05;
  const KEY_ROT_DEG = 5, KEY_ROT_FINE_DEG = 1;
  const DRAG_ROT_RAD_PER_PX = 0.01;
  const WHEEL_MOVE_A_PER_UNIT = 0.01, WHEEL_ROT_RAD_PER_UNIT = 0.003;
  const GRAB_RADIUS_A = 1.5, GRAB_MIN_PX = 12; // how close to an atom counts as "on the ligand"

  // `enabled()` (optional) can switch ligand control off, e.g. while the ligand is hidden.
  function attach({ viewer, el, M, center, getPose, getXyz, refresh, enabled = () => true }) {
    const ctl = { pivot: null };

    // Model-space directions of screen right / down / toward-viewer at the ligand centroid,
    // from the pseudo-inverse of the local model->screen Jacobian (units: Å per pixel).
    function screenBasis() {
      const c = M.add(center, getPose().t);
      const toP = (p) => ({ x: p[0], y: p[1], z: p[2] });
      const s0 = viewer.modelToScreen(toP(c));
      const cols = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map((e) => viewer.modelToScreen(toP(M.add(c, e))));
      const Jx = cols.map((s) => s.x - s0.x), Jy = cols.map((s) => s.y - s0.y);
      const a = M.dot(Jx, Jx), b = M.dot(Jx, Jy), d = M.dot(Jy, Jy), det = a * d - b * b;
      const pinv = (sx, sy) => {
        const u = (d * sx - b * sy) / det, v = (-b * sx + a * sy) / det;
        return M.add(M.scale(Jx, u), M.scale(Jy, v));
      };
      const right = pinv(1, 0), down = pinv(0, 1);
      const up = M.scale(down, -1);
      return { right, down, up, toward: M.normalize(M.cross(right, up)) };
    }

    function translateScreen(dx, dy) {
      const B = screenBasis(), pose = getPose();
      pose.t = M.add(pose.t, M.add(M.scale(B.right, dx), M.scale(B.down, dy)));
    }
    function translateDepth(angstrom) {
      const pose = getPose();
      pose.t = M.add(pose.t, M.scale(screenBasis().toward, angstrom));
    }
    // Rotate by Q about a model-space pivot (default: the ligand centroid). With x' = R(x-c)+c+t,
    // Q(x'-p)+p = QR(x-c) + c + [Q(c+t-p) + p - c], so R <- QR and t <- Q(c+t-p) + p - c.
    function rotate(axisModel, angle, pivot) {
      if (M.norm(axisModel) === 0 || angle === 0) return;
      const Q = M.axisAngle(axisModel, angle), pose = getPose();
      if (pivot) {
        const arm = M.sub(M.add(center, pose.t), pivot);
        pose.t = M.sub(M.add(M.matVec(Q, arm), pivot), center);
      }
      pose.R = M.orthonormalize(M.matMul(Q, pose.R));
    }
    function rotateScreen(dx, dy, pivot) {
      const B = screenBasis();
      const axis = M.add(M.scale(M.normalize(B.up), dx), M.scale(M.normalize(B.right), dy));
      rotate(axis, Math.hypot(dx, dy) * DRAG_ROT_RAD_PER_PX, pivot);
    }

    // PyMOL-like: drag ON the ligand moves it, Shift+drag (anywhere) or right-drag on it rotates it
    // (around the atom under the cursor, if any, else around the centroid),
    // Ctrl/Cmd+drag (anywhere) moves it. Everything else goes to the 3Dmol camera
    // (left-drag rotate, middle-drag pan, right-drag / scroll zoom).
    // Index of the ligand atom under the cursor (closest on screen within the grab radius), or -1.
    function atomUnder(pageX, pageY) {
      const xyz = getXyz();
      if (!xyz) return -1;
      const pxPerA = 1 / Math.max(M.norm(screenBasis().right), 1e-9);
      const r = Math.max(GRAB_MIN_PX, GRAB_RADIUS_A * pxPerA);
      let best = -1, bestD = r;
      viewer.modelToScreen(xyz.map((p) => ({ x: p[0], y: p[1], z: p[2] }))).forEach((s, i) => {
        const d = Math.hypot(s.x - pageX, s.y - pageY);
        if (d <= bestD) { best = i; bestD = d; }
      });
      return best;
    }
    const overLigand = (pageX, pageY) => atomUnder(pageX, pageY) >= 0;
    const pt = (ev) => (ev.touches ? ev.touches[0] : ev);
    function actionFor(ev) {
      const p = pt(ev);
      if (ev.shiftKey) return "rotate";
      if (ev.ctrlKey || ev.metaKey) return "move";
      if (!overLigand(p.pageX, p.pageY)) return null;
      return ev.button === 2 ? "rotate" : "move";
    }
    let drag = null;
    function onDown(ev) {
      if (!enabled()) return;
      if (ev.touches && ev.touches.length > 1) return;
      if (ev.button === 1) return; // middle button: always camera pan
      const a = actionFor(ev);
      if (!a) return;
      ev.stopPropagation(); ev.preventDefault();
      drag = { a, x: pt(ev).clientX, y: pt(ev).clientY };
      // Rotating with the cursor on an atom turns the ligand around that atom (it stays in place),
      // so an end that is already placed can be kept fixed. Elsewhere: around the centroid.
      if (a === "rotate") {
        const i = atomUnder(pt(ev).pageX, pt(ev).pageY);
        if (i >= 0) { drag.pivot = getXyz()[i].slice(); ctl.pivot = drag.pivot; refresh(); }
      }
      el.style.cursor = "grabbing";
    }
    function onMove(ev) {
      if (!drag) {
        if (!ev.touches && ev.target && el.contains(ev.target)) {
          el.style.cursor = ev.shiftKey ? "crosshair" : overLigand(ev.pageX, ev.pageY) ? "grab" : "";
        }
        return;
      }
      ev.stopPropagation(); ev.preventDefault();
      const p = pt(ev), dx = p.clientX - drag.x, dy = p.clientY - drag.y;
      drag.x = p.clientX; drag.y = p.clientY;
      if (drag.a === "move") translateScreen(dx, dy); else rotateScreen(dx, dy, drag.pivot);
      refresh();
    }
    const onUp = () => {
      if (!drag) return;
      const hadPivot = !!ctl.pivot;
      drag = null; ctl.pivot = null; el.style.cursor = "";
      if (hadPivot) refresh();
    };
    const onContext = (ev) => ev.preventDefault();
    el.addEventListener("mousedown", onDown, true);
    el.addEventListener("touchstart", onDown, { capture: true, passive: false });
    root.addEventListener("mousemove", onMove, true);
    root.addEventListener("touchmove", onMove, { capture: true, passive: false });
    root.addEventListener("mouseup", onUp, true);
    root.addEventListener("touchend", onUp, true);
    // Right-drag on the ligand rotates it, so suppress the context menu over the viewer.
    el.addEventListener("contextmenu", onContext, true);
    // Shift+scroll pushes/pulls the ligand along the view axis, Alt+scroll spins it; plain scroll zooms.
    const onWheel = (ev) => {
      if ((!ev.shiftKey && !ev.altKey) || !enabled()) return;
      ev.stopPropagation(); ev.preventDefault();
      const delta = ev.deltaY || ev.deltaX; // Shift+wheel arrives as deltaX in some browsers
      if (ev.shiftKey) translateDepth(-delta * WHEEL_MOVE_A_PER_UNIT);
      else rotate(screenBasis().toward, delta * WHEEL_ROT_RAD_PER_UNIT);
      refresh();
    };
    el.addEventListener("wheel", onWheel, { capture: true, passive: false });

    const onKey = (ev) => {
      if (ev.target.closest && ev.target.closest("input, textarea, select")) return;
      if (ev.ctrlKey || ev.metaKey || ev.altKey) return; // leave browser shortcuts (Ctrl+F, Cmd+R, …) alone
      if (!enabled()) return;
      const step = ev.shiftKey ? KEY_STEP_FINE_A : KEY_STEP_A;
      const rot = ((ev.shiftKey ? KEY_ROT_FINE_DEG : KEY_ROT_DEG) * Math.PI) / 180;
      const B = screenBasis();
      const px = step / Math.max(M.norm(B.right), 1e-9); // pixels that correspond to `step` Å
      const k = ev.key.toLowerCase();
      const moves = {
        w: () => translateScreen(0, -px), s: () => translateScreen(0, px),
        a: () => translateScreen(-px, 0), d: () => translateScreen(px, 0),
        r: () => translateDepth(step), f: () => translateDepth(-step),
        arrowleft: () => rotate(B.up, -rot), arrowright: () => rotate(B.up, rot),
        arrowup: () => rotate(B.right, -rot), arrowdown: () => rotate(B.right, rot),
        q: () => rotate(B.toward, rot), e: () => rotate(B.toward, -rot),
      };
      if (!moves[k]) return;
      ev.preventDefault();
      moves[k]();
      refresh();
    };
    root.addEventListener("keydown", onKey);

    function detach() {
      el.removeEventListener("mousedown", onDown, true);
      el.removeEventListener("touchstart", onDown, true);
      root.removeEventListener("mousemove", onMove, true);
      root.removeEventListener("touchmove", onMove, true);
      root.removeEventListener("mouseup", onUp, true);
      root.removeEventListener("touchend", onUp, true);
      el.removeEventListener("contextmenu", onContext, true);
      el.removeEventListener("wheel", onWheel, true);
      root.removeEventListener("keydown", onKey);
      el.style.cursor = "";
    }
    return Object.assign(ctl, { screenBasis, rotate, detach });
  }

  root.LigandControls = { attach };
})(typeof window !== "undefined" ? window : this);
