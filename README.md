# ai_brain: artificial life from the bottom up

Goal: grow life-like behavior from simple physics and chemistry rules, then
cells, then nervous systems. Nothing "wants" anything; behavior must emerge
from energy rules + copying with errors + selection.

## Ladder (roadmap)

| Stage | Unit | Question | Status |
|---|---|---|---|
| 1. Chemistry | particles with pair rules | Do membranes self-assemble? | in progress |
| 2. Protocell | membrane bubble + autocatalysis | Can a bubble grow and split? | planned |
| 3. Cell world | cells with genome, energy, light | Does phototaxis evolve? Engulfing? | planned |
| 4. Endosymbiosis | cell inside cell | Does "keep instead of digest" win? | planned |
| 5. Multicell + neurons | cells with ion channels, spikes | Do signaling cells specialize? | planned |
| 6. Brain | spiking networks, STDP, dopamine | Does a creature learn in its lifetime? | planned |

## Stage 1: membrane self-assembly

Model: 2D coarse-grained lipids (Cooke-Kremer-Deserno style), implicit water.

- Lipid = 1 head bead (likes water) + 2 tail beads (oily), bonded by springs.
- Water is not simulated as particles. Its effect lives in the rules: oily
  tails attract each other (the hydrophobic effect), heads do not.
- Thermal noise (Langevin thermostat) jiggles everything.

Nobody codes "make a membrane". Tails hide from water by clumping, heads stay
outside, and micelles, bilayers, and vesicles (rings in 2D) appear.

The chemistry lives in `alife/chem/world.py` as the `Rules` table: data, not code.

## Run

```bash
pip install -r requirements.txt
python experiments/e01_membrane.py              # live window
python experiments/e01_membrane.py --headless   # saves PNG snapshots to runs/
```

Useful flags: `--lipids 300`, `--kT 0.6` (heat), `--stick 1.0` (tail attraction),
`--bend 50` (stiffness), `--tails 3`, `--seed 1`.

### Findings so far

- kT 1.0: too hot. Lipids stay loose, only tiny clumps.
- kT 0.4: too cold. Glassy blobs with heads trapped inside.
- kT 0.6 + stiff 3-tail lipids: micelles and bilayer fragments (tails parallel
  inside, heads on both edges). Default setting.
- Floppy lipids (`--bend 10`) pack badly: shape matters, not just stickiness.

## Environment: cycles and shocks

A planet is never steady. `alife/env.py` wraps any World in an `Environment`:

- seasons: Earth orbits the sun. Temperature (`kT`) swings slowly (`amp`, `period`).
- day/night: Earth spins. `light` goes 0..1 and back (`day_period`). Sunlight
  warms matter (`day_heat`). Lipids feel light only as heat; a pigment
  particle (later) will absorb it and use the energy.
- currents: uneven sun heating makes water turn in swirls (`current`); the
  pattern drifts like weather, stronger by day. Beads are dragged toward the
  local water speed. Diffusion alone spreads matter only as sqrt(time).
- shocks: rare random asteroid. Heat spike that fades + outward blast near
  the impact point, which mixes matter that was stuck.

Idea under test: steady worlds get stuck. Heating breaks bad structures,
cooling lets them reform better (annealing). Shocks shake loose jams.

Experiment 2 compares three worlds with the same lipids. Each run ends with
a calm stretch at base temperature, so all are judged equally:

```bash
python experiments/e02_environment.py --env shocks              # live window
python experiments/e02_environment.py --env steady --headless   # also: cycle
python experiments/e02_environment.py --env shocks --day-heat 0.1   # + day/night sun
python experiments/e02_environment.py --env shocks --day-heat 0.1 --current 0.4   # + currents
python experiments/e02_environment.py --plot                    # compare all logs
```

Vesicle = sealed water pocket inside one membrane cluster
(`water_pockets` in `alife/chem/analysis.py`).

### Findings (3 seeds each, partial runs for steady/cycle/shocks/sun)

- Steady water packs tightest (lowest energy).
- Only full weather (seasons + sun + currents + asteroids) ever made vesicles:
  4 times, each broke again within ~1000 steps. Disturbance makes and breaks.
- Currents keep clusters smaller (10-12 vs 7).
- Limit: one 60 nm droplet. A real asteroid is local; Earth has billions of
  ponds where survivors reseed. Next idea: many ponds linked by rare transport.

## Elements and charge

`alife/chem/elements.py` is the sim's periodic table: each bead type has a
charge, a size (with its water shell for ions), and how oily it is. Adding an
element is adding a row. Types: head, tail, fa_head (fatty acid, -1), Na, K,
Cl, Mg, Ca, pigment, mineral (frozen surface).

Charges use screened Coulomb: `coulomb * qi * qj * exp(-r/debye) / r`.
`debye` ~2 is a pond, ~0.5 is seawater (salt hides charges).

```bash
python experiments/e03_ions.py --world ca     # also: neutral, na, mg
python experiments/e03_ions.py --plot
```

### Findings (100k steps, 3 seeds, averages)

| world | clusters | largest | Ca/Mg bridges | energy |
|---|---|---|---|---|
| neutral heads | 9.5 | 45 | - | -6.27 |
| fatty acid + Na | 11.9 | 37 | 0 | -5.98 |
| + Ca2+ | 10.5 | 34 | 5.4 | -6.13 |
| + Mg2+ | 12.7 | 38 | 1.9 | -5.95 |

- Charged heads repel: fatty acids make smaller, rounder clusters than neutral lipids.
- Ca2+ grabs heads and bridges lipids (glue). Mg2+ mostly doesn't: its thick
  water shell keeps it away. Matches real chemistry, from size + charge alone.
- No vesicles in any element world yet.

## Big connected world

Real Earth is one connected place, but not the same everywhere, and not
evenly connected. Open ocean dilutes everything; fully closed pools trap
discoveries. Life likely needed something in between.

`experiments/e04_world.py`: one 192 x 192 box (3 x 3 zones, 1800 lipids).

- `Terrain` (in `alife/env.py`): each zone has its own warmth and sunlight,
  blended smoothly. Hot vent left-middle, cold top-right, shaded spots.
- Each bead feels its local temperature (`World.kT_field`).
- Rock walls (`add_rock_wall`, frozen `rock` beads): an almost closed pool,
  an open pool, a ridge with a strait.
- Asteroids are local (`shock_local=True`): heat and blast near the impact,
  the rest of the world survives.
- Log column `zones`: lipids per zone, to see matter move and refill.

```bash
python experiments/e04_world.py              # live (slow, ~15 ms/step)
python experiments/e04_world.py --headless   # runs/e04_world/seed0
python experiments/e04_world.py --fatty      # fatty acids + Na+
```
