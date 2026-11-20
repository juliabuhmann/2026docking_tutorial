// The diffusion sampler of AlphaFold 3, written to be read by students. No DOM access, so
// tests/test_diffusion.js can run it with node.
//
// Follows sample(), random_augmentation(), random_rotation() and noise_schedule() in
// src/alphafold3/model/network/diffusion_head.py of github.com/google-deepmind/alphafold3
// (commit 5870e09dd26a1541bd622fbdfa28dd518cc56192), with the constants of its SampleConfig.
//
// The one deliberate difference: AlphaFold 3's denoiser is a trained network that sees the noisy
// atoms, the noise level and the trunk's single and pair representations (from sequence, MSA,
// templates). It has to *predict* the clean structure. Here the denoiser is an oracle that returns
// the crystal structure, optimally superposed onto the noisy atoms. Everything around it, the
// schedule, the extra noise, the step size and the random rotations, is as in AlphaFold 3.
(function (root) {
  "use strict";

  const AF3 = {
    SIGMA_DATA: 16.0, S_MIN: 0.0004, S_MAX: 160.0, P: 7, // noise_schedule()
    STEPS: 200, GAMMA_0: 0.8, GAMMA_MIN: 1.0, NOISE_SCALE: 1.003, STEP_SCALE: 1.5, // SampleConfig, eval
    NUM_SAMPLES: 5, // AF3 draws 5 samples per seed and ranks them by confidence; one is shown here
  };

  // sigma(t) for t in [0, 1]: from SIGMA_DATA * S_MAX (2560 A) down to SIGMA_DATA * S_MIN (0.0064 A)
  function noiseSchedule(t) {
    const { SIGMA_DATA, S_MIN, S_MAX, P } = AF3;
    return SIGMA_DATA * (S_MAX ** (1 / P) + t * (S_MIN ** (1 / P) - S_MAX ** (1 / P))) ** P;
  }
  const noiseLevels = (steps = AF3.STEPS) => Array.from({ length: steps + 1 }, (_, i) => noiseSchedule(i / steps));

  // seeded random numbers, so a seed always gives the same trajectory
  function mulberry32(a) {
    return () => { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  }
  function normal(rng) { // Box-Muller
    const u = Math.max(rng(), 1e-12), v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  // random_rotation(): Gram-Schmidt on two random normal vectors; rows e0, e1, e2
  function randomRotation(rng) {
    const v0 = [normal(rng), normal(rng), normal(rng)], v1 = [normal(rng), normal(rng), normal(rng)];
    const n0 = Math.max(1e-10, Math.hypot(...v0)), e0 = v0.map((x) => x / n0);
    const d = v1[0] * e0[0] + v1[1] * e0[1] + v1[2] * e0[2];
    const w = v1.map((x, k) => x - e0[k] * d), n1 = Math.max(1e-10, Math.hypot(...w)), e1 = w.map((x) => x / n1);
    const e2 = [e0[1] * e1[2] - e0[2] * e1[1], e0[2] * e1[0] - e0[0] * e1[2], e0[0] * e1[1] - e0[1] * e1[0]];
    return [e0, e1, e2];
  }

  // random_augmentation(): centre on the mean, rotate (positions @ rot), add a N(0, 1 A) translation
  function augment(X, rng) {
    const n = X.length / 3, c = [0, 0, 0];
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) c[k] += X[3 * i + k] / n;
    const R = randomRotation(rng), t = [normal(rng), normal(rng), normal(rng)];
    for (let i = 0; i < n; i++) {
      const p = [X[3 * i] - c[0], X[3 * i + 1] - c[1], X[3 * i + 2] - c[2]];
      for (let j = 0; j < 3; j++) X[3 * i + j] = p[0] * R[0][j] + p[1] * R[1][j] + p[2] * R[2][j] + t[j];
    }
  }

  // Superpose flat coordinate array A onto B (least squares, Compare.superpose); returns A moved.
  function superposeFlat(A, B, superpose, apply) {
    const n = A.length / 3, P = [], Q = [];
    for (let i = 0; i < n; i++) { P.push([A[3 * i], A[3 * i + 1], A[3 * i + 2]]); Q.push([B[3 * i], B[3 * i + 1], B[3 * i + 2]]); }
    const { R, t } = superpose(P, Q), out = new Float64Array(A.length);
    P.forEach((p, i) => { const q = apply(R, t, p); out[3 * i] = q[0]; out[3 * i + 1] = q[1]; out[3 * i + 2] = q[2]; });
    return out;
  }

  // Run the sampler. ref: flat crystal coordinates (the oracle's answer). geom: {superpose, apply}
  // from js/compare.js. Returns, per step i = 0..steps: the noise level, gamma, t_hat, and the
  // positions superposed onto the crystal for display (AF3 itself keeps no fixed frame).
  function sample(ref, geom, { seed = 1, steps = AF3.STEPS } = {}) {
    const rng = mulberry32(seed), levels = noiseLevels(steps), n = ref.length / 3;
    let X = new Float64Array(ref.length);
    for (let i = 0; i < X.length; i++) X[i] = levels[0] * normal(rng); // positions *= noise_levels[0]
    const frames = [superposeFlat(X, ref, geom.superpose, geom.apply)], info = [{ sigma: levels[0], gamma: 0, tHat: levels[0] }];
    for (let s = 1; s <= steps; s++) {
      const prev = levels[s - 1], level = levels[s];
      augment(X, rng);
      const gamma = level > AF3.GAMMA_MIN ? AF3.GAMMA_0 : 0;
      const tHat = prev * (1 + gamma);
      const scale = AF3.NOISE_SCALE * Math.sqrt(Math.max(tHat * tHat - prev * prev, 0));
      const noisy = new Float64Array(X.length);
      for (let i = 0; i < X.length; i++) noisy[i] = X[i] + scale * normal(rng);
      const denoised = superposeFlat(ref, noisy, geom.superpose, geom.apply); // the oracle (see top)
      const dT = level - tHat;
      for (let i = 0; i < X.length; i++) X[i] = noisy[i] + AF3.STEP_SCALE * dT * (noisy[i] - denoised[i]) / tHat;
      frames.push(superposeFlat(X, ref, geom.superpose, geom.apply));
      info.push({ sigma: level, gamma, tHat });
    }
    return { frames, info, n };
  }

  // Root-mean-square distance of the atoms from their centre: how big the cloud is.
  function spread(X) {
    const n = X.length / 3, c = [0, 0, 0];
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) c[k] += X[3 * i + k] / n;
    let s = 0;
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) s += (X[3 * i + k] - c[k]) ** 2;
    return Math.sqrt(s / n);
  }

  const api = { AF3, noiseSchedule, noiseLevels, randomRotation, sample, spread, mulberry32 };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Diffusion = api;
})(typeof window !== "undefined" ? window : this);
