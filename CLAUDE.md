# Guide for coding sessions

Read README.md first (roadmap, how to run, findings). This file is the
context a new session needs to continue.

## Goal

Artificial life from the bottom up: physics and chemistry rules, then
protocells, then evolving cells (phototaxis, engulfing / endosymbiosis),
then neurons and spiking brains with local learning (STDP, dopamine), not
backprop. Nothing is coded to "want" anything: behavior must emerge from
energy rules + copying with errors + selection.

## How the user works

- Learning the science while building. Show the mechanism first (picture,
  small interactive demo, or a short run), then turn it into rules and code.
- Ideas come as intuitions ("asteroids shake things up", "the world is
  connected"). Take them seriously, correct the science gently and briefly
  (e.g. motivation emerges from selection; the sun doesn't orbit Earth),
  then build the closest honest version.
- Wants plain explanations, short. Ask before starting long reruns; check
  saved outputs first.
- Runs heavy jobs on their own server (headless). Repo is private:
  github.com/mindaniel/is_this_life.

## Code map

- `alife/chem/elements.py`: the periodic table of bead types (charge,
  hydrated size, oiliness, color). Add an element = add a row.
- `alife/chem/world.py`: engine. Rules built from the table; WCA repulsion,
  cos^2 oily attraction, screened Coulomb, bonds/bends, Langevin thermostat
  with water flow drag (`flow`), per-position temperature (`kT_field`),
  frozen beads (rocks). Verlet neighbour lists (short + charged long list).
  `add_bond` / `break_bond` for reactions (molecule ids merge/split).
  `repel_inside` = bead types that also repel beads of their own molecule
  (floppy sugar chains folded through themselves and crashed without it).
- `alife/chem/reactions.py`: `Formose` (carb beads join into sugar chains
  with Ca2+ catalyst, chains split), kinetic Monte Carlo every 50 steps.
- `alife/chem/membrane.py`: `Leak`, small molecules hop across lipid walls
  by size/charge (C1 fast, ions slow, C5+ never). `Leak.wet` (pos -> bool)
  rejects hops that would land on dry land; e05 sets it to the pond.
- `alife/chem/molecules.py`: lipids (neutral or fatty acid), ions, rock walls,
  mineral surface.
- `alife/chem/fast.py`: numba versions of the force loops and a cell-list
  pair search. World uses them automatically when numba is installed.
- `alife/chem/analysis.py`: tail clusters; `water_pockets` = vesicle detector
  (sealed water inside one lipid cluster, area >= 20).
- `alife/env.py`: Environment (seasons, day/night light, drifting currents,
  asteroid shocks global or local, calm ending), Terrain (zone map of
  warmth and sun), Pond (one or many ponds with own timing; shore that
  shrinks/refills; rain swirls; full neighbours overlap), Crust (molecules
  on dry land freeze until rewetted; reactions and leaks skip frozen beads).
  Weather (4 seasons per `year`, steady wind = gyre in each pond, random
  storms; `carry_at` adds the water-pressure push so heavy beads ride
  the current instead of being flung to the shore). e05 `--weather`.
- `alife/viz.py`: matplotlib view; temperature map; red circles on vesicles;
  sugar carbons coloured by chain length; `legend()`; `Movie` (mp4 via
  ffmpeg, or .gif). e05 `--video out.mp4` films a run (the user likes
  watching runs: offer a video for new experiments).
- `experiments/e01..e04`: membrane, environment cycles, ions, big world.
  Each has live mode and `--headless`; e02/e03 have `--plot`.
  e04 logs vesicle x,y,area (`pockets` column) and saves `frame_*.npz`
  whenever one is detected; `e04_check.py` tracks lifetimes and draws a
  zoomed gallery of them.
- `experiments/e05_wetdry.py`: wet-dry pond. `--sugar`, `--chem` (formose
  with rain formaldehyde, Ca2+, meteorites), `--leak`, `--name`, `--ponds N`
  (`--r-wet` > spacing/2 so they connect), `--crust`, `--video`. One pond
  with default flags reproduces old runs bit for bit (checked on e08 seed 7). Pocket log
  per vesicle: x,y,area,lining,cover,trapped big sugar,small sugar C,ions.
- `experiments/e07_reactor.py`: one ring vesicle, seeds + Ca2+ inside, food
  outside; `--sealed`, `--ca none|outside` controls; `--video`.
- `experiments/e08_compete.py`: two ring vesicles (R 12), sugar in A and/or B,
  `--drip N` adds one loose lipid every N steps; logs lipids per ring, areas.
  Ring layers are counted at their middle radius (head-radius counting
  crammed the outer tails and some starts blew up).
- `alife/basin.py` + `experiments/e11_basin.py`: landscape level in real
  units (m, hours, mol): ponds on hills, rain/evaporation/springs, streams
  and floods, lipid/formaldehyde/sugar per pond, wet-dry film. Says which
  ponds reach the CVC. Bead world is ~100 nm, so landscape questions go here.
- 3D: `World(dim=3)` (2D stays bit-identical; numba kernels take 2 or 3
  columns, `cell_pairs3`). `molecules.add_lipid3 / bilayer_disc3 / vesicle3 /
  gas3`, `analysis.water_pockets3` (validated on spheres/disc/sheet).
  `World.rescale(f)` = water leaves/returns (molecules move as a whole).
  `e12_3d.py` (disc closing, self-assembly, vesicle lifetime; edge-on slice
  video), `e13_drywet3d.py` (wet-dry by box shrink/grow, lipid in mM).
- `alife/chem/nitrogen.py` + `experiments/e14_nitrogen.py`: HCN/NH3 beads (end of
  element table), Strecker -> glycine (nh3-carb-hcn molecule), glycolonitrile,
  HCN polymer (pure-HCN chains), hydrolysis retypes HCN -> NH3. e14 `--video`
  films it. In the pond: e05 `--nitrogen` (rain brings `--rain-hcn`/`--rain-nh3`
  per refill; an always-wet pond gets the same rain spread evenly). Formose and
  `chain_lengths` skip carbs bonded to other atoms (`pure_carb`); lone HCN/NH3 leak
  through membranes like formaldehyde.
- `experiments/e06_formose.py`: formose alone in a small box (checks lag,
  autocatalysis, Ca needed, heat).
- `experiments/e15_split3d.py`: 3D grow-and-split (Zhu & Szostak 2009). One
  vesicle, `--mix` of its lipids short-tail (`--short-oily`), sugar osmolyte
  inside and out (`--osm`), feeds loose short-tail lipids (`--feed`/`--batch`,
  `--near`), optional shear. Logs biggest cluster, pockets + volumes, length
  (eigenvalue ratio, 1 = ball); saves `*_split_<t>.npz` when 2+ pockets.
  Note: stail at oiliness 1.0 = normal tail.

## Rules of thumb (learned the hard way)

- Lipids: kT 0.6, 3 tails, bend 50. kT 1.0 too hot, 0.4 glassy.
- Random starts overlap: always `world.step(200, max_force=50)` before real
  physics (builders in experiments do this).
- Validate every detector on synthetic shapes and on a real snapshot before
  trusting numbers. The vesicle detector gave false positives twice (gaps in
  jammed clumps; loose ions as walls).
- New force terms: check force == -dE/dx by finite difference.
- Same seed = same history. Compare worlds over several seeds (3+).
- Parallel runs: `OMP_NUM_THREADS=1`, one process per core. Background
  shell kills may leave python processes alive; check before claiming stopped.
- `runs/` is git-ignored (outputs). Log what you will need later (e.g.
  vesicle positions), the pictures alone can't be re-analysed.
- Speed: numba kernels in `alife/chem/fast.py` (forces + cell-list neighbour
  search), ~4x faster: ~4.5 ms/step for the 1800-lipid world (8 with ions).
  `world.fast = False` uses the numpy reference (`forces_numpy`); after any
  force change, check both agree (they match to 1e-13 now).
  Flow and kT_field are re-read every `world.field_every` (10) steps.
- Numba is only in `C:\Apps\Python314\python.exe`; plain `python` (3.10) has none
  and silently uses the numpy forces (same numbers, slower). The 2026-10 e05
  nitrogen runs used 3.10.
- Machine has 64 cores: run many seeds at once, one process each.
- Always send run output to a file (`> runs/stdout/x.txt 2>&1`) and grep it
  for Traceback: runs sent to /dev/null crashed silently.
- New beads/molecules must go on free spots (cKDTree check, gap >= 1):
  overlaps blow up bonds.
- Vesicle lining scores drift to ~0.78 on real rings: use 0.7/0.7 cutoffs
  when counting, not the strict 0.8/0.9 defaults.
- If nothing can cross a membrane, any trapped bead (ion, formaldehyde)
  acts as osmotic pressure. Use `--leak` for realistic results.
- Any jump or placement (leak hops, new beads) must land in water: the
  uncapped shore spring flings a bead on dry land back into crowds and the
  run blows up (e05 seed 2 crash, 2026-10-30).

## Where things stand (2026-10-22)

- Membranes self-assemble (bilayer strips, discs). Vesicles are rare and
  short-lived in a single droplet; only full weather made any.
- Fatty acids (charged) make smaller round clusters; Ca2+ bridges heads,
  Mg2+ mostly doesn't (thick water shell). No vesicles with ions yet.
- Big connected world (e04), 8 seeds x (neutral, fatty), checked by eye:
  detections are real head-lined water cores but tiny (~3 beads across),
  mostly embedded in clumps of bilayer patches, median life <= 1000 steps.
  Clean single-ring vesicles: a few per 8 seeds. Fatty-acid pockets trap Na+.
  `lining()` in analysis.py scores membrane quality (validated on synthetic
  rings: ~0.95 score, cover ~1.2-1.5; tangle gap: 0.68 / 0.61).
- (2026-10-22) Wet-dry pond (e05) makes ~25x more vesicles than always-wet;
  trapped sugar keeps them open (osmotic). Formose works (e06: lag, then
  autocatalysis; seeding skips the lag). Pond chemistry + leaky membranes:
  meteorite sugar seeds formose in the pond; sugar-filled vesicles last
  longest (77600 steps). Details and tables in README.

- (2026-10-23) Lipid stealing (e08_compete, 8 seeds x 3 setups, 200k):
  no sign. Sugar-filled rings swell (area ~700 vs ~550) but gain no more
  lipids than empty ones; gains are 0-11, same as noise. Dripped lipids
  mostly clump with each other, not into rings. Rings (even empty) open
  after ~10-150k steps. 3-tail lipids almost never leave a membrane;
  2-tail ones leave all the time but can't hold a ring (hop test).

- (2026-10-26) e09 river basin (`--layout basin --ponds 9 --land-seed 3`,
  3 seeds, 160k steps): ponds on rolling ground, big lake at the lowest
  spot, streams to the nearest lower pond. Expected lipids to run down into
  the lake; the opposite happened: the lake fell from ~628 to 90-290 lipids,
  and the lipids piled into one or two small ponds (a different pond each seed,
  even the "isolated" one, 317 -> 917 via floods). Wind drift + floods beat
  the stream current. Vesicles about the same as the 4-pond loop (21% vs 20% of logs),
  but bigger clumps (largest cluster 682 vs 293). Speedup: `Pond.gaps` does
  shoreline trig only for near (bead, pond) pairs (bit-identical).

- (2026-10-27) e11 basin (3 years): lipid passes the CVC only a few hours a
  year, only in hot ponds as they nearly dry. 3D (e12, kT 1.1, 3 seeds):
  ready-made vesicles R8/R10 stay sealed 100k steps; flat discs R10-22 stay
  flat (none closed in 100k); 1700 gas lipids only merge into patches
  (500-860) in 200k, no vesicle. e13 15 mM: wet-dry makes one big branched
  membrane (largest ~1400 vs ~140 always-wet) but no vesicle on rewetting
  (box only 54 nm, network spans it). Film now = close-up cutaway of the
  biggest cluster.

- (2026-10-27) e13 now logs as it goes and saves `*_pocket_<step>.npz`. Flood reruns
  (box 60/100/140 x 3 seeds, last-wet 200k) started detached; one 140 run showed a
  410-volume pocket before (unchecked). e14 nitrogen alone works (README table):
  dilute HCN makes no polymer, NH3 picks glycine, HCN alone self-starts via NH3.

- (2026-10-30) Nitrogen in the e05 pond (README table): formaldehyde eats the HCN
  (glycine + glycolonitrile); with `--rain-carb 5` HCN polymerises, wet-dry ~2x
  the polymer of always-wet, chains up to 83. e13 floods (rewet box 100/140)
  close the dried lipid stack into one big 3D vesicle in 5/6 runs (~140-180k
  steps after rewet, stays sealed); rewet 60 never. No-dry control (e13
  `--no-dry`, box 100/140 x 3 seeds, 290k) done: 0/6 vesicles, only small
  clusters (largest 69-94 lipids). Drying matters, not just dilution.

- (2026-10-30) e08 with half fatty-acid walls (`--mix 0.5 --hop-oily 1.8`,
  5 seeds): still no stealing (with drip both rings gain the same; without
  drip nobody leaves a wall). e10: fed ring grows to ~210-290 lipids, never
  splits in 2D. e05 seed 2 crash fixed (`Leak.wet`); rerun finished 200k cleanly.
  e15 (3D grow + split, 18 runs, README table): feeding grows the wall
  (2225 -> ~3050) but the vesicle stays a ball (volume follows; 35 sugar beads
  too weak to hold water back). Shear stretches it into a tube that tears and
  wraps across the periodic box; "2 pocket" logs are narrow spots in one tube,
  checked by eye. No split. Next: stronger osmolyte, box long along shear.

## Next steps (user's choice; ask before starting a big one)

1. (done) e04 vesicles confirmed: real but tiny and short-lived.
2. (done) Wet-dry cycles, sugar, formose, leaky membranes.
2b. (done 2026-10-23: yes, 6 -> 27.6 sugar C vs ~6-10 in controls) Test the "vesicle as reactor" guess: does formaldehyde leak into
   sugar-filled compartments and get built into trapped sugar? Then:
   swelling -> taking up lipids -> growth -> splitting (stage 2).
3. Pigment chemistry: a light-absorbing bead that gets excited by `env.light`
   and either releases heat or drives a reaction (bond forming / pumping).
4. Stage 2 protocell: vesicle that grows by taking in lipids and splits
   (shear from currents can split it), something trapped inside.
5. Later stages per README ladder: genome, selection, phototaxis, engulfing,
   neurons.
