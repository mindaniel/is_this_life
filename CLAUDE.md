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
  by size/charge (C1 fast, ions slow, C5+ never).
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
- `experiments/e06_formose.py`: formose alone in a small box (checks lag,
  autocatalysis, Ca needed, heat).

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
- Machine has 64 cores: run many seeds at once, one process each.
- Always send run output to a file (`> runs/stdout/x.txt 2>&1`) and grep it
  for Traceback: runs sent to /dev/null crashed silently.
- New beads/molecules must go on free spots (cKDTree check, gap >= 1):
  overlaps blow up bonds.
- Vesicle lining scores drift to ~0.78 on real rings: use 0.7/0.7 cutoffs
  when counting, not the strict 0.8/0.9 defaults.
- If nothing can cross a membrane, any trapped bead (ion, formaldehyde)
  acts as osmotic pressure. Use `--leak` for realistic results.

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
