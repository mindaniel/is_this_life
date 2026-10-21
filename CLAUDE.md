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
- `alife/chem/molecules.py`: lipids (neutral or fatty acid), ions, rock walls,
  mineral surface.
- `alife/chem/analysis.py`: tail clusters; `water_pockets` = vesicle detector
  (sealed water inside one lipid cluster, area >= 20).
- `alife/env.py`: Environment (seasons, day/night light, drifting currents,
  asteroid shocks global or local, calm ending) and Terrain (zone map of
  warmth and sun).
- `alife/viz.py`: matplotlib view; temperature map; red circles on vesicles.
- `experiments/e01..e04`: membrane, environment cycles, ions, big world.
  Each has live mode and `--headless`; e02/e03 have `--plot`.

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
- Speed: ~2 ms/step for 200 lipids, ~15 ms/step for the 1800-lipid world.

## Where things stand (2026-10-21)

- Membranes self-assemble (bilayer strips, discs). Vesicles are rare and
  short-lived in a single droplet; only full weather made any.
- Fatty acids (charged) make smaller round clusters; Ca2+ bridges heads,
  Mg2+ mostly doesn't (thick water shell). No vesicles with ions yet.
- Big connected world (e04): far more vesicle detections (not yet confirmed
  by eye; new snapshots circle them), bigger clusters, lipids drift out of
  the hot vent zone.

## Next steps (user's choice; ask before starting a big one)

1. Confirm e04 vesicles with circled snapshots, more seeds on the server.
2. Wet-dry cycles: evaporation shrinks the pond (crowding), rain refills.
3. Pigment chemistry: a light-absorbing bead that gets excited by `env.light`
   and either releases heat or drives a reaction (bond forming / pumping).
4. Stage 2 protocell: vesicle that grows by taking in lipids and splits
   (shear from currents can split it), something trapped inside.
5. Later stages per README ladder: genome, selection, phototaxis, engulfing,
   neurons.
