// Menu: one card per chapter, with a small hand-drawn illustration each. Then boot.
(function () {
  "use strict";
  const { $, esc, el } = PW;
  const svg = (b) => `<svg viewBox="0 0 260 132" preserveAspectRatio="xMidYMid slice" xmlns="http://www.w3.org/2000/svg">${b}</svg>`;
  const blob = (cx, cy, r, n, seed) => {
    let d = "";
    for (let i = 0; i <= n; i++) {
      const a = (i / n) * Math.PI * 2, rr = r * (1 + 0.18 * Math.sin(a * 3 + seed) + 0.08 * Math.cos(a * 5 + seed * 2));
      d += `${i ? "L" : "M"}${(cx + Math.cos(a) * rr).toFixed(1)} ${(cy + Math.sin(a) * rr).toFixed(1)}`;
    }
    return d + "Z";
  };

  const ART = {
    hunt: svg(`<path d="${blob(130, 70, 52, 60, 1)}" fill="#1f3a5c" stroke="#8fb8ff" stroke-width="2"/>
      <path d="M150 34 C134 52 136 74 152 86 C160 72 166 54 150 34Z" fill="#0c1426" stroke="#8fb8ff" stroke-width="1.5"/>
      <g fill="#7ee787"><circle cx="150" cy="58" r="4"/><circle cx="146" cy="68" r="4"/><circle cx="153" cy="74" r="3.5"/></g>
      <circle cx="96" cy="88" r="6" fill="#60a5fa" opacity=".9"/><circle cx="118" cy="44" r="6" fill="#ffb454" opacity=".9"/>
      <circle cx="150" cy="64" r="15" fill="none" stroke="#5eead4" stroke-width="2" stroke-dasharray="4 3"/>`),
    algo: svg(`${(() => { let s = ""; for (let x = 20; x < 260; x += 14) for (let y = 12; y < 132; y += 14) s += `<circle cx="${x}" cy="${y}" r="1.3" fill="#33436a"/>`; return s; })()}
      <path d="${blob(130, 68, 50, 60, 2.3)}" fill="#1f3a5c" stroke="#b69cff" stroke-width="2" opacity=".95"/>
      <circle cx="122" cy="70" r="15" fill="rgba(255,209,102,.18)" stroke="#ffd166" stroke-width="1.6"/>
      <circle cx="140" cy="60" r="10" fill="rgba(255,209,102,.18)" stroke="#ffd166" stroke-width="1.6"/>
      <g fill="#ffd166"><circle cx="122" cy="70" r="2"/><circle cx="140" cy="60" r="2"/></g>
      <g stroke="#5eead4" stroke-width="1.3" opacity=".9"><path d="M60 70 H200"/><path d="M130 14 V122"/><path d="M80 20 L182 118"/></g>`),
    dock: svg(`<path d="M40 110 C40 40 110 20 130 60 C150 20 220 40 220 110 Z" fill="#1f3a5c" stroke="#8fb8ff" stroke-width="2"/>
      <path d="M112 110 C112 82 148 82 148 110Z" fill="#0c1426"/>
      <g transform="translate(150 30) rotate(25)"><g stroke="#ff8fd0" stroke-width="4" stroke-linecap="round" fill="none">
        <path d="M0 0 L14 8 L14 24 L0 32 L-14 24 L-14 8 Z"/><path d="M0 32 L0 46"/></g>
        <g fill="#60a5fa"><circle cx="-6" cy="52" r="4"/><circle cx="6" cy="52" r="4"/></g></g>
      <path d="M128 94 C120 80 150 64 160 56" stroke="#5eead4" stroke-width="1.5" stroke-dasharray="4 3" fill="none"/>
      <text x="196" y="122" fill="#7ee787" font-size="13" font-family="monospace">1.4 Å</text>`),
    morph: svg(`<g opacity=".45"><path d="M60 100 C60 50 90 30 110 40 L120 96 Z" fill="#1f3a5c" stroke="#8fb8ff" stroke-width="1.5"/>
      <path d="M200 100 C200 50 170 30 150 40 L140 96 Z" fill="#1f3a5c" stroke="#8fb8ff" stroke-width="1.5"/></g>
      <path d="M92 104 C92 60 112 46 126 52 L128 100 Z" fill="#244a73" stroke="#5eead4" stroke-width="2"/>
      <path d="M168 104 C168 60 148 46 134 52 L132 100 Z" fill="#244a73" stroke="#5eead4" stroke-width="2"/>
      <circle cx="130" cy="76" r="7" fill="#7ee787"/>
      <path d="M70 30 Q90 18 104 30" stroke="#ffd166" stroke-width="2" fill="none"/>
      <path d="M190 30 Q170 18 156 30" stroke="#ffd166" stroke-width="2" fill="none"/>`),
    cofold: svg(`${(() => { let s = ""; const N = 22; for (let i = 0; i < N; i++) { const t = i / (N - 1);
        const x = 40 + t * 180, y = 66 + Math.sin(t * 9) * 28 * (1 - t * .6), jx = (Math.sin(i * 12.9) * 30) * (1 - t), jy = (Math.cos(i * 7.3) * 26) * (1 - t);
        s += `<circle cx="${(x + jx).toFixed(1)}" cy="${(y + jy).toFixed(1)}" r="${(2.4 + t * 1.8).toFixed(1)}" fill="${t < .5 ? "#3a4c74" : "#60a5fa"}" opacity="${(.35 + t * .65).toFixed(2)}"/>`; } return s; })()}
      <path d="M150 66 C168 40 186 40 200 58 C212 74 196 96 178 92" stroke="#60a5fa" stroke-width="3" fill="none" stroke-linecap="round"/>
      <g fill="#7ee787"><circle cx="186" cy="70" r="3.4"/><circle cx="178" cy="62" r="3.4"/><circle cx="192" cy="60" r="3"/></g>
      <text x="34" y="122" fill="#66748f" font-size="11" font-family="monospace">noise → structure</text>`),
  };

  function render() {
    const wrap = $("#chapter-cards");
    PW.get("index", "index").then((index) => {
      wrap.innerHTML = "";
      PW.chapters.forEach((c) => {
        const p = c.progress ? c.progress(index) : null;
        const card = el(`<button class="card" style="--c:${c.color}">
          <div class="art">${ART[c.id] || ""}</div>
          ${p ? `<span class="prog ${p.done ? "done" : ""}">${esc(p.text)}</span>` : ""}
          <span class="num">${esc(c.num)}</span>
          <h2>${esc(c.title)}</h2>
          <p>${c.blurb}</p>
          <span class="learn">${esc(c.learn)}</span>
          <span class="go">Play</span></button>`);
        card.addEventListener("click", () => PW.open(c.id));
        wrap.appendChild(card);
      });
    });
  }

  window.PWMenu = { render };
  window.addEventListener("DOMContentLoaded", () => {
    $("#btn-unlock").addEventListener("click", () => { PW.progress.unlockAll = true; PW.save(); PW.toast("All levels unlocked.", "good"); render(); });
    $("#btn-reset").addEventListener("click", async () => {
      const ok = await PW.modal("<h2>Reset progress?</h2><p>Stars and completed tasks in this browser are deleted.</p>",
        [{ label: "Cancel", value: false }, { label: "Reset", primary: true, value: true }]);
      if (ok) { PW.resetProgress(); render(); }
    });
    PW.boot();
    render();
  });
})();
