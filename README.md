# Docking tutorial

A browser game for the ETH Zürich Bioinformatics Concept Course (Protein Bioinformatics II):
where does a small molecule bind, and how would a computer know? Five chapters, all on real
structures from the PDB. No installation, no server: it is a static site.

**Status: first draft.** Every chapter has been exercised in headless Chromium (see
[What has been checked](#what-has-been-checked)), and the numerical parts have tests. It has not yet
been tried by students, on touch devices, or in Safari or Firefox.

| Chapter | What students do | Lecture slides |
|---|---|---|
| 1 · Pocket hunt | Click where the ligand binds on 10 proteins, easy to hard: open clefts, a dimer interface, buried cavities, an allosteric site, the "undruggable" KRAS pocket. Three views (surface, atoms, cartoon) and a slice. Short skippable intro cards on the protein and ligand before each level; hot/cold feedback, stars, a story per protein after it. | 21 |
| 2 · Pocket algorithms | Run a grid method (LIGSITE-like) and alpha spheres (fpocket-like) live on the same proteins, change the parameters and the score, benchmark on all ten. Eight challenges. Opens with a four-card illustrated intro to both methods (2D sketches computed on a toy protein). | 22–29 |
| 3 · Dock it yourself | Move and rotate the crystal ligand back into its pocket, watching symmetry-corrected RMSD, clashes and contacts. Five ligands, three difficulty modes. Ported from the `docking_game` prototype. | 32–35 |
| 4 · Apo to holo | Morph between two crystal structures (adenylate kinase, maltose-binding protein, the Abl DFG flip), watch the pocket appear, count ligand clashes per frame. | 31 |
| 5 · Co-folding | Short intro cards on co-folding, then step through AlphaFold 3's own diffusion sampler (its 200-step schedule, extra noise and random rotations; only the network is replaced by the crystal structure, and the page lists what differs), fold the complex for real in the browser with [LocalFold](https://localfold.org), then upload the prediction to get its pocket-aligned ligand RMSD. | 39, 57 |

The look and the menu-of-chapters structure were inspired by
[Cellworks](https://dennisgankin.github.io/cellworks/); no code was taken from it.

A guided tour (coach marks) introduces the buttons on the first visit to the menu and to the game
screen, plus each chapter's own controls the first time that chapter opens. The **Tour** button
(top bar) and **Take the tour** (menu) replay it. What has been seen is stored with the progress.

## Run it

Open `docs/index.html` in a browser. The data files are `<script>`s rather than JSON, so this works
straight from disk. Or serve the folder: `python -m http.server -d docs 8000`.

Two things are loaded from CDNs: [3Dmol.js](https://3dmol.csb.pitt.edu) 2.5.3 (always) and
[fflate](https://github.com/101arrowz/fflate) 0.8.2 (only when a zip is uploaded in Chapter 5).

## Publish on GitHub Pages

Push the repository, then **Settings → Pages → Deploy from a branch → `main`, folder `/docs`**.
`docs/.nojekyll` makes GitHub serve the files unchanged.

## Data

Nothing structural is typed in by hand. `scripts/build_data.py` reads the four files in `configs/`,
downloads everything from RCSB (entries, ligand definitions, sequences, citations; cached in
`data_cache/`, not committed) and writes `docs/data/`:

| Config | Output | Contents |
|---|---|---|
| `proteins.json` | `prot_<id>.js` | receptor, ligand with bond orders and symmetry, pocket residues (≤ 4.5 Å of the ligand), decoy sites, polar contacts, and the alpha spheres (Delaunay circumspheres, radius 2.8 to 8 Å) |
| `morphs.json` | `morph_<id>.js` | both structures, residues matched by sequence alignment, superposed on the rigid core |
| `cofold.json` | `cofold_<id>.js`, `jobs/<id>_job.json` | deposited sequence, CA coordinates by sequence position, AlphaFold-Server-style job file |
| `docking.json` | copied into `index.js` | docking cases and feedback thresholds |

```bash
conda activate ../conda_envs/docking      # rdkit, gemmi, scipy, numpy
python scripts/build_data.py
```

The texts in `proteins.json` and `morphs.json` (stories, intro cards, hints, decoy feedback) are written for
students. Each structure's own paper is fetched from RCSB and shown with it.

## Tests

```bash
node tests/test_pockets.js                  # pocket ranks, Chapter 2 challenges solvable, morph endpoints
node tests/test_metrics.js                  # pose RMSD maths (from the prototype), writes random poses
python tests/check_rmsd_rdkit.py            # ... and checks them against RDKit CalcRMS
python tests/check_automorphisms.py         # every symmetry matching is a real graph automorphism
python tests/make_compare_fixtures.py       # fake predictions from the crystal structures
node tests/test_compare.js                  # Chapter 5 scoring on them: 0.00 A, 3.00 A, and a warning
node tests/test_diffusion.js                # AF3 schedule endpoints, rotations, sampler ends on the crystal
```

All pass as of this commit. `tests/out/` is generated and not committed.

### How the two pocket finders do on the ten proteins

Rank of the true pocket (pocket centre within 4 Å of the ligand centre, the usual DCC
criterion; for the streptavidin tetramer, of any of its four biotins) with the default settings; generated by `tests/test_pockets.js`:

| Protein | Grid, volume | Grid, burial | Alpha, size | Alpha, hydrophobicity | Alpha, burial |
|---|---|---|---|---|---|
| Streptavidin + biotin (a tetramer) | #5 | #1 | #11 | #3 | #7 |
| Carbonic anhydrase II + acetazolamide | #1 | #1 | #23 | #14 | #11 |
| HIV-1 protease + indinavir (a dimer) | #1 | #1 | #1 | #1 | #1 |
| Trypsin + benzamidine | #2 | #2 | #1 | #3 | #1 |
| Abl kinase + imatinib | #1 | #3 | – | – | – |
| Acetylcholinesterase + donepezil | #1 | #1 | #1 | #1 | #1 |
| Estrogen receptor + estradiol | #1 | #1 | #2 | #1 | #2 |
| T4 lysozyme L99A + benzene | #2 | #1 | #3 | #1 | #1 |
| PTP1B + an allosteric inhibitor | – | – | #1 | #1 | #1 |
| KRAS G12C + sotorasib | – | – | #1 | #1 | #1 |
| **Rank 1** | 5/10 | 6/10 | 5/10 | 6/10 | 6/10 |

The defaults (grid: 1 Å spacing, 1.2 Å probe, 6 of 7 directions within 8 Å; alpha spheres:
fpocket's 3.4 to 6.2 Å radii, with a simpler clustering than fpocket's: single linkage at 2 Å, at least 10 spheres) were picked with `tests/explore_pockets.js` so that
neither method is perfect and the misses are instructive. The best setting found by
`tests/search_benchmark.js` reaches 7/10. Both are simplified teaching versions of the published
methods, not reimplementations.

## What has been checked

In headless Chromium 1254 (SwiftShader WebGL), driving the page with real mouse, keyboard and
file-upload events:

- the menu and all five chapters load without console errors (desktop size; the menu and Chapter 2 also at 390 px wide)
- Chapter 1: clicks pick the atom under the cursor (also through the slice), give hot/cold
  feedback, award and save stars, unlock the next level
- Chapter 2: both methods, scoring switches, the benchmark, and challenge completion
- Chapter 3: keyboard moves, RMSD readouts, submit and history
- Chapter 4: slider frames (about 30 ms each), pocket search per frame, clash counts
- Chapter 5: the illustration, and upload of the synthetic PDB and mmCIF predictions

Not checked yet:

- a **real** LocalFold result. The parser follows LocalFold's PDB writer as of its commit
  `7f66b52` (protein numbered 1…N along the input sequence, ligand as HETATM with CCD atom names),
  and the job files load in LocalFold's own `jobFromJson`, but no fold has been run and uploaded.
  Folding needs a WebGPU browser; there is none on the cluster.
- a real AlphaFold Server mmCIF (tested only with gemmi-written mmCIF)
- touch interaction, Safari, Firefox, and speed on student laptops
- whether the stories are pitched right for the course

Added later (guided tour, intro cards, surface view), checked in the same headless Chromium:
the menu tour, the Chapter 1 and 2 tours and intro cards (keyboard and click), the three views in
Chapter 1 including slicing and reveal in surface view, the surface view in Chapter 2, and the
tour and intro cards at 390 px width with touch taps, without page errors. Not checked: how the
surface view performs on student laptops (it is slow in software WebGL), and whether the intro
texts are pitched right.

A second round of review, one independent reviewer per chapter (facts and usability, driving the
page in headless Chromium at 1400, 1280 and 390 px), found and led to fixes of, among others: a
crash in docking case 4, stars that ignored aids switched off before Submit, buried pockets that
could be solved from outside, a grid score that never changed a ranking, and AlphaFold Server
offered as a route although it cannot take these ligands. The review did not cover real touch
hardware, Safari/Firefox, real GPUs, or a real LocalFold fold.

## Known simplifications

- The pocket finders are simplified (see `docs/js/pockets.js`); "found" means DCC ≤ 4 Å, which
  is strict for long ligands such as imatinib.
- Morph frames are a straight-line interpolation between two crystal structures, not a pathway.
- Chapter 5 runs AlphaFold 3's sampler (`docs/js/diffusion.js`, following `diffusion_head.py` of
  google-deepmind/alphafold3 at commit `5870e09`) with an oracle in place of the network: the crystal
  structure superposed onto the noisy atoms. So it always ends on the crystal. One sample instead of
  five, crystal heavy atoms only, and each frame superposed onto the crystal for display.
- Streptavidin is shown as its tetramer (biological assembly 1, built with gemmi) in Chapters 1 to 3,
  and its four biotin sites all count. Chapter 5 folds one chain.

## Credits

Structures: RCSB PDB. 3D viewer: 3Dmol.js. In-browser co-folding: LocalFold (Steinegger,
Ovchinnikov and contributors). Chapter 3 started as the `docking_game` prototype in this workspace.
Algorithms: LIGSITE (Hendlich et al. 1997), fpocket (Le Guilloux et al. 2009).
