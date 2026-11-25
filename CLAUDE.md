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
- Plans a YouTube channel about the project: after every experiment, add to
  `JOURNEY.md` (episode log: question, what happened, caveats, footage paths).
  Film new experiments. `runs/` (videos) is not in git: remind to back it up.
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
  `--charged f`: fraction of short-tail heads that are COO- (`fa_head`, counter-ions
  implicit via Debye screening); 0 draws no extra random numbers, so old seeds replay.
  `e15_check.py` tells real daughters from rooms in one body (wall cluster, start
  sugar per pocket, head/tail layering between pockets, `--png` cut).
  Note: stail at oiliness 1.0 = normal tail. `--box-x` = long box (x only).
- Long box: `World(box=(Lx, Ly, Lz))` makes `w.L` an array (a number = cube, bit-identical
  to before, checked by replaying e15d). Kernels take `w.box()` (length per axis).
- `reactions.LipidJoin` + `experiments/e16_catalyst.py`: "cat" bead (end of element
  table) joins two single-chain lipids head to head into a two-chain lipid
  (phospholipid-like); no fuel yet (flagged). e16 = vesicle A with catalysts vs B with
  as many sugar beads; logs chains / two-chain lipids per wall, A_in/B_in; film = slice + counts.
  Use `--length 1.2 --fa-only` (default bond squeezes joined heads: see status).
- `reactions.Protonation` (constant-pH Monte Carlo, kernel `fast.titrate`): fatty-acid heads
  switch `fa_head` (COO-) <-> `fa_acid` (COOH) every 50 steps by Metropolis (pH, pKa 4.8,
  Coulomb, pair energy). `elements.PAIR_EPS` = per-pair attraction (fa_acid-fa_head hydrogen
  bond), `Rules.set_pair_eps`. e16 `--ph --pka --hbond`, logs A_q/B_q (charged share).
  `experiments/check_protonation.py` = energy + titration-curve checks.
  `Protonation.closed` (heads facing a sealed inside: no swap with outside, only `hand_over` hops
  between neighbouring heads, `fast.hand_over`); `Protonation.facing(w)` finds them
  (`water_pockets3(points=, near_water=True)`). `reactions.Burn` (cat burns fuel -> waste, H+ onto
  a COO-), `reactions.PHShocks` (outside pH jumps + relax). Bead `fa_oh` = fatty alcohol (H-bond
  to fa_head, never joins). e16 flags `--sealed-h --fuel --p-burn --shocks --relax --alcohol
  --acid-oily`; logs A_qin/B_qin, A_oh/B_oh, ph, burns. e15 `fa_head(..., alcohol)`.

- Save / carry on / join: `World.save(path, **extra)` / `World.load(path)` (pickle; flow, kT_field,
  external not saved; carrying on replays bit for bit, `check_save.py`). `World.merge(a, b, axis)`
  cuts each open at its emptiest slab and puts b after a (a's rules, time, markers first), then
  step(200, max_force=50). `experiments/merge_worlds.py a b out`, `experiments/hcn_pond.py`
  (concentrated HCN + NH3 pond, `--like` copies cross-section/kT/rules). e16 `--save-every`,
  `--start state.pkl` (+ `--steps` more, `--reseed`), `--hcn/--nh3` (adds beads, runs Nitrogen at
  kT0 = world kT), `--glow` (dark HCN chains `nitrogen.dark_beads` absorb light: `--glow-heat`,
  photo-acid `--glow-gain`, via `LightPump(pigment=, any_head=True, heat=)`). e16 `--rain` uses
  `World.rescale(f, group, unwrap=True)` (BFS unwrap; old e13 floods used no unwrap: big networks
  wider than half the box may have torn a little).
- Bead `tar` (end of table): `Nitrogen(tar_len=n)` retypes HCN chains of n+ to tar (oily, value a guess;
  e16 `--tar-len --tar-oily`, logs dark_tail). `Rules.grow()` (called by `World.load`) adds rows for element
  types added after a save. `experiments/e21_tar.py` = tar branches from e20's joined worlds.
- `alife/chem/sun.py` (`Sun(d, period, kT_earth)`): smooth day brightness (pi/2 max(0, sin)/d^2) and kT(t)
  (kT_earth d^-0.5, small lagged swing). `LightPump(sun=)` scales p_abs by brightness, gain/heat in Earth kT.
  e16 `--sun-d --kT-earth` sets w.kT every 50 steps; `hcn_pond.py --kT0`, `merge_worlds.py --kT`. `experiments/e22_sun.py`.

- `alife/chem/copying.py` + `experiments/e23_copy.py`: RNA-like copying. Monomer = backbone `nbb` + base `nx`/`ny`
  (end of element table). `Copying`: Monte Carlo base pairing (X-Y eps, X-X/Y-Y eps_mis, one partner per base,
  flagged), pair springs (base-base + backbone to partner base), stacking bonds between neighbour bases, joining of
  touching backbone ends whose bases are within stack_reach (counts templated / correct), hydrolysis. e23 flags
  `--eps --k-pair --stack --p-join --p-break --cycle --hot --dry --dry-f`; logs joins/templated/ok, `_chains.csv`
  (copy score vs shuffled).
- `experiments/e24_stretch.py` (+ `e24_plot.py`): osmotic stretch test, bag with `--osm-in` sugar vs `--osm` outside,
  probe carbs inside, no Leak; logs area per lipid, holes, probes/sugar inside.
- Beads `pep` (water-loving) / `pho` (oily) = peptide (end of element table). `molecules.peptide3` = stiff two-stripe
  ladder (rungs, diagonals sqrt2, bends), `peptide_points` = its bead positions. `experiments/e25_peptide.py`
  (`--start wall|water`, `--peptides 0` = control; wall start leaves room for each peptide, ~5 lipids) + `e25_plot.py`.
- Bead `wat` (end of table) = explicit water (~4 molecules). `experiments/e26_water.py`: flat periodic wall + water beads,
  `--ew` water-water/water-polar eps (0.5; 1 too sticky), counts crossings vs the wall's LOCAL middle (`local_offset`),
  `--bundle k` ready-made peptide rings, `--start pkl` carry on. ~64 ms/step for 6k beads. `e26_plot.py`, `e26_look.py`.

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
  checked by eye. No split.
  e15b (`runs/e15b/`, sugar 0.15/0.3, food anywhere `--near 0 --clear 1.0`,
  box 80 for shear): fed vesicles gain wall while water stays/shrinks, go
  long (shape ~1.9 vs ~1.5 unfed), stay sealed; no thread, no split. Shear 0.1
  tears the tube flat. ~1/3 of food makes its own blobs. Look script:
  slices along each pocket's long axis (scratchpad `e15_look.py`).

- (2026-11-01) e15c (`runs/e15c/`, food x5/x10 per 200 steps, current 0.02/0.05 in
  box 80): current 0.05 divides the vesicle in 3/3 seeds (pinch at ~72-88k, two
  sealed sacs to the end, each with start sugar, ~8% leaked by 200k). The sacs still
  touch (two bilayers, one cluster): not yet drifted apart. No current: food piles
  into thick lipid masses that make 2-4 rooms in one body (not a split). Check
  pockets by wall cluster + head/tail layering along the line between them, not by
  pocket count (touching vesicles = one cluster).
  Filmed 400k rerun (seed 1, `runs/e15c/long/`): split at ~88k, the two sacs drift
  back together and fuse at ~225-235k, then the vesicle stretches and splits again
  (~300k). Division happens but isn't final; ~15% of start sugar leaked by 400k.
  `experiments/e15_check.py` = the pocket/layering/sugar check.

- (2026-11-10) e15d (`runs/e15d/`, 400k, `--charged 0.5` x 3 seeds vs neutral x 3):
  the "fusion" is a periodic-box artefact. Daughters ride opposite bands of the sin
  current and move apart to half the box (~40), then meet again across the wrap
  (`runs/e15d/drift.png`). In an open pond they'd stay apart. Charge: water layer
  between charged daughters' walls (neutral: wall to wall), splits later (93-118k vs
  82-88k), merges later; weak effect. Every run: split, merge, split again. Sugar loss
  at 400k 39-98/173 (opening/resealing). e15_check counts fa_head as head now.

- (2026-11-11) e16 first try (`runs/e16/`): catalyst vesicle collapsed and lost all
  contents in 6/6 runs. Artefact: head-head bond (r0 1.0) squeezed pairs to 0.89 vs 1.15
  between neighbours, so joins shrank the inner leaflet. Fix: per-bond rest length
  (`add_bond(i, j, r0)`, `World.rest`), `LipidJoin(length=1.2, chains=(STAIL,))`, e16
  `--length 1.2 --fa-only`, logs A_in/B_in. Rerun (6 runs, `e16fix_*` logs) started.

- (2026-11-11) e16 fixed rerun: A stays sealed 6/6, but no win: fed A ends 41-69 chains
  behind B (3/3; controls -3..+119). Two-chain lipids form wall-spanning patches, flip to the
  outer leaflet. e16 logs A_fa/B_fa and A_left/B_left (fatty acids leaving each wall). Leaving
  rises with lower stail oiliness (0.85: 1.4%/1000 steps, 0.7: 8%, 0.65 leaks); A's extra
  leaving comes from its inner leaflet (into its own lumen), outer about equal. Beads lack
  the real effect (two-chain lipids order the wall and hold fatty acids tighter). Scratch
  scripts: leaflet_hops.py (per-leaflet leaving), pack.py (chain order per leaflet).

- (2026-11-13) `--fa-lengths "1:w,2:w,3:w"` (e15 `fa_tails`/`parse_lengths`, e16 flag, logs
  A_len/B_len): random fatty-acid tail lengths. 30k tests: walls keep the long ones (mean 2.0 in ->
  2.2 in wall), leaving 18% (25:50:25) / 41% (50:30:20, B opens once) per 1000 steps vs 1.5% before;
  A and B leave alike. 200k catalyst rerun with the 25:50:25 mix started 2026-11-13
  (12 runs, `e16mix_*` logs, CSVs `*_mix1-2-3_*`).
  Result (2026-11-14): all sealed, still no win. A-B at 200k: catalyst unfed +35/-2/+24, fed
  -33/+25/+2; controls unfed +6/-50/+16, fed -2/+7/-98. ~55% of A's fatty acids joined, but
  leaving per fatty acid same in A and B (13-16%/1000). Beads lack two-chain ordering of the wall.

- (2026-11-14) e15e long box (x 160, 300k, 3 charged + 3 neutral): split at ~75-100k, sacs
  move apart (to 60-80) without merging, but stay joined by a thin tether in 6/6 (thins to
  ~14 tail beads per unit length, then holds; never snaps). 2/6 end with one sealed sac.
  `runs/e15e/neck.png`. Neck measured along the min-image line: past half the box it reads
  the empty side (check the wall profile along x instead).

- (2026-11-15) e17 protonation (`runs/e17/`, 6 x 30k, pH 4/5/6.5 x pairing 0/1, no catalyst):
  all sealed. Pairing cuts loose fatty acids (225-316 vs 326-382) but leaving per 1000 steps is
  higher (25% vs 20%) with trips half as long; no half-charge sweet spot. Walls half charged at
  pH ~5 (real 7-9: no low-dielectric interface), so compare by charged share. 3 seeds agree.

- (2026-11-16) e18 (`runs/e18/`, 30k): pH shocks to 2 or 8 never break a bag; inner heads lag the
  outer ones ~2-3k steps (brief pH difference, then leaks by flip-flop). Fuel burnt by A's catalyst
  gives only a small short acid dip inside (0.04-0.06 vs B, 2/3 seeds): H+ leaks out as fast as made.
  Fatty alcohol 30%: slightly steadier wall. `--acid-oily 0.5`: walls fully neutral, bag squeezes to
  a thick oily shell (volume 900 -> 260); strength likely too high. Heat from fuel not done yet.

- (2026-11-16) e19 (`runs/e19/`, 24 runs): heat, held pH 2/8 (30k), random shocks, rain: no bag
  breaks. Hand-placed pump (30 pigments in A, light 10k on/off): A's inner charged share 0.2-0.3
  below B while lit, 6/6, back at night; A ends 30-80 fatty acids below B. e20 (`runs/e20/`,
  save + merge, 12 runs): pond HCN polymer doubles after joining (~1300 dark beads), touches walls,
  but heat/photo-acid change nothing in the bags; dilute HCN ("wet") ~50x less polymer.
  Untested: sticky (oily) HCN polymer that would sit in walls.

- (2026-11-17) e21 (`runs/e21/`, tar_len 5, oily 0.5/0.8, from e20 dry joins, 9 runs): tar sinks into walls
  (84-172 beads touching tails vs ~7 non-oily), unevenly; in light the tar-rich bag's inside goes acid by day
  (-0.05..-0.09 charged share vs other bag, dips to -0.17), 6/6; dark control ~0. Emergent crude pump.
- (2026-11-18) e22 habitable zone (`runs/e22/`, `alife/chem/sun.py` smooth sun, d 0.7/0.85/1/1.2/1.5 x 3 seeds,
  e16 `--sun-d`): bags sealed everywhere; tar and tar-in-walls fall with distance (~4x hot vs cold); pump
  (tar-rich minus other qin, day) -0.04..-0.08 at 0.85/1, ~0..-0.03 at 1.2/1.5. d 0.7 (3/3) and 0.85 (1/3)
  crash at peak light: photon heat (10 kT x brightness into 8 beads) overheats tar clumps, integrator
  blows up (checked: heat 0 survives). Fix before rerunning: cap/spread photon heat or smaller dt.
- (2026-11-20) e22w (`e22_sun.py --warm`, e16 `--glow-warm`, `LightPump(warm=True)` -> `World.warm` = extra kT per
  bead for 50 steps, same energy as the kick; reactions ignore it): the same 4 runs crash (0.7 x3, 0.85 s1, at
  sun 2-3x). The heat amount (10 kT/photon) is the problem, not delivery. Tar/loose trends repeat; pump vs
  distance does not (seed spread ~0.05). `e22_sum.py/e22_plot.py --warm` select `_h10w_` runs.
  Options: heat 0 at all d (sun kT already includes the light's warming), smaller heat, or dt 0.005.
  User decision (2026-11-20): not needed. Goal is life on Earth; Venus/other planets were a side question.
  Earth runs (d 1) all finish fine. Don't chase the Venus crash.
- (2026-11-20) Earth 9 seeds (d 1, warm heat, seeds 3-8 with own e16 bags in runs/e20/pre): pump day -0.036 +- 0.033,
  night -0.003 +- 0.006, 6/9 pump; size follows tar gap between the two walls (r -0.85; even tar = no pump).
  `experiments/e22_earth.py` -> `runs/e22/e22w_earth.png`. Next offered: gated channels.

- (2026-11-23) e23 copying (5 seeds x 100k): ~13% of joins template-guided, ~99% correct, but random joins win ~7:1;
  warm/cool cycles ~9%, drying doubles both (11.5% vs 13.3%). No full copies yet. Idea: short pieces (oligos) hold on
  longer than monomers.
- (2026-11-23) e24 stretch (4 sugar levels x 3 seeds x 60k): wall stretches up to ~7% area but never opens (0 holes,
  1/480 probes out). Implicit water: no water inflow, so tension may be too weak.
- (2026-11-24) e24 harder push (inside 0.8 / 1.6 x 3 seeds): wall tears (area +20% / +60%, holes 15-33%),
  ~90% of sugar and nearly all probes leave, then the bag shrinks and reseals (3/6 fully by 60k). Stretch =
  safety valve, threshold between 0.4 and 0.8; not selective. Next options: peptide pores, oligos for copying.
- (2026-11-24) e25 peptides (40 per bag, wall / water / control x 3 seeds x 60k): no pores. In the wall they stay upright
  and alone (groups <= 3), some slide out to lie on the surface; from water none get in (land on the surface, clump in
  water). Probes 39-40/40 inside like the control. Likely cause: implicit water makes a buried polar bead cheap and a
  channel brings no water reward; also no voltage. Options: polar-in-oil penalty, explicit water beads, voltage.

- (2026-11-25) e26 explicit water (`runs/e26full/`, ew 0.5, 4 setups x 3 seeds x 100k + 30k recount): no water crosses
  without peptides or at 1:40 peptide:lipid; at 1:17 all 3 seeds leak (~5-6 water beads / 1000 steps after a 40-80k wait)
  around tipped-over peptides. Ready-made 6-rings fall apart in ~10k (leak, then slow/stop). Threshold "carpet" leak, no
  barrel pore. Next options: voltage, bag version, peptides entering from water.

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
6. (agreed 2026-11-10, after e16) Fuel for the e16 catalyst (used-up food molecule;
   by-products heat or gas that feed back), then light-driven charge pump pigment,
   then gated channels (stretch first; light, heat, chemical, voltage later). The
   user's idea: "ports" that light/heat open so charge can pass; reactions making
   heat or gas that change what happens next. README "Next protocell steps".
