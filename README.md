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

Check vesicles by eye: `python experiments/e04_check.py runs/e04_world`
tracks each pocket over time and draws `vesicle_gallery.png` (zoomed crops,
with a lining score: do lipid heads face into the pocket, and how densely).

### Findings (100k steps, 8 seeds each, 1000-step logs)

| world | pockets detected | lived >= 5000 steps | membrane-lined |
|---|---|---|---|
| neutral lipids | 44 | 1 | 4 |
| fatty acid + Na+ | 33 | 2 | 1 |

- Detections are real water cores lined by lipid heads, but tiny
  (area 20-40, ~3 beads across) and short-lived (median <= 1000 steps).
- Most sit inside a bigger clump of bilayer patches, not a free single ring.
  Clean single-ring vesicles are rare (a few per 8 seeds).
- Fatty-acid pockets often trap Na+ ions inside: first "something inside".
- So: vesicles form in the big world, but they don't grow or last.
  Bigger, longer-lived ones probably need a push (wet-dry cycles, more lipid).

## Wet-dry pond

Ponds on early Earth dried out and refilled. Drying crowds everything
together; rain brings it back apart.

`experiments/e05_wetdry.py`: one pond of 300 lipids. `Pond` (in
`alife/env.py`) pushes beads back inside a shore whose radius shrinks
(46 -> 22) and grows back every `--period` steps. `--rain` adds swirls while
it refills.

```bash
python experiments/e05_wetdry.py --pond wetdry --headless   # runs/e05_wetdry/wetdry_seed0
python experiments/e05_wetdry.py --plot                     # compare runs
python experiments/e05_wetdry.py --chem --leak --meteor-sugar 10 --meteor-every 40000     --seed 7 --video runs/video/rain_meteor_seed7.mp4         # watch a run (one frame per 200 steps)
```

Many ponds and shore crust: `--ponds 4 --r-wet 50` puts 4 ponds on a 2x2
grid of land, each drying on its own timing (start near full, speeds differ
by up to 10%). Full neighbours overlap and share water, so molecules can
move between ponds. `--crust`: molecules stranded on dry land stick there
(frozen, no reactions, no leaking) and dissolve back when the water returns.

```bash
python experiments/e05_wetdry.py --ponds 4 --r-wet 50 --crust --chem --leak --meteor-sugar 10     --meteor-every 40000 --video runs/video/ponds.mp4
```

`--video` saves an mp4 (needs ffmpeg; `.gif` works without it). Colours:
orange = lipid tails, red = heads, pale / mid / dark purple = formaldehyde /
small sugar / big sugar, grey = Ca2+, green = Cl-, red circles = vesicles,
faded = dried into the shore crust. Default 10 frames/s: one dry-refill
cycle = 10 s.

### Findings (200k steps, 8 seeds, vesicle = lining >= 0.7 and cover >= 0.7)

| pond | vesicles | present | longest (steps) |
|---|---|---|---|
| always wet | 43 | 2% | 6000 |
| wet-dry | 1022 | 34% | 10200 |
| wet-dry + gentle rain | 1034 | 34% | 11200 |
| wet-dry + strong rain | 749 | 19% | 7800 |

- Drying merges the lipids into one sheet; refilling leaves it folded with
  sealed water inside. Wet-dry makes ~25x more vesicles than an always-wet pond.
- Small rings (radius 3-5) break; rings of radius 8+ survive. Nothing holds
  the inside open, so most still collapse.

Osmotic pressure: dissolved stuff that can't cross the membrane pushes
water in and the vesicle open. `--sugar N` adds N sugar beads.

| wet-dry pond | vesicles | present | lived >= 2000 | longest |
|---|---|---|---|---|
| 150 sugar | 2815 | 90% | 461 | 65000 |
| 300 sugar | 2072 | 84% | 368 | 43400 |

- Drying traps sugar inside folds; refilling dilutes the outside. The
  inside > outside imbalance comes by itself.
- Long-lived vesicles are bigger and hold more sugar. Too small + too much
  sugar = it bursts (osmotic lysis).
- The structure is one lump with many sugar-filled compartments, not free
  single vesicles.

## Sugar from formaldehyde (formose)

Instead of adding sugar by hand, make it. Formaldehyde (CH2O) forms in the
air (lightning, UV) and rains into ponds. With Ca2+ as catalyst it joins
into sugars: the formose reaction.

`alife/chem/reactions.py`: each carbon is one `carb` bead; bonds form and
break while the world runs (`World.add_bond` / `break_bond`).

- join: two chain ends that touch bond. Needs Ca2+ nearby, faster when hot.
  C1 + C1 is very slow (the hard first step); C1 + chain is fast.
- split: a chain of 4+ can break into two pieces of 2+ (C4 -> 2 C2).

The cycle is not coded. Grow + split together should multiply sugar.

```bash
python experiments/e06_formose.py --headless --name base   # runs/e06_formose/base_seed0
python experiments/e06_formose.py --plot                   # compare.png
```

### Findings (100k steps, 8 seeds, 400 formaldehyde, 20 Ca2+)

| condition | C1 left | C2-C4 | C5+ (ribose-sized) |
|---|---|---|---|
| base (kT 0.6) | 14% | 69% | 18% |
| cold (kT 0.5) | 96% | 1% | 3% |
| hot (kT 0.7) | 0% | 93% | 7% |
| no Ca2+ | 97% | 3% | 0% |
| 5 C2 added at start | 0% | 80% | 20% |

- Slow start, then sugar explodes (autocatalysis), like real formose.
  Half the formaldehyde used by step 26k-88k (one run never started).
- A few C2 at the start skip the wait: all runs take off by 14k-21k.
- No catalyst or cold: almost nothing. Heat is faster but breaks big sugars.

## Pond chemistry and leaky membranes

Everything in one pond: `e05_wetdry.py --chem`. Each refill brings
formaldehyde (`--rain-carb`), Ca2+ sits in the water (`--ca`), and
meteorites sometimes drop ready-made C5 sugars (`--meteor-sugar`,
`--meteor-every`).

`--leak` (`alife/chem/membrane.py`): real membranes leak, by size and
charge. formaldehyde fast >> C2 > C3-C4 >> ions; C5+ never. Done as Monte
Carlo hops that only count if they cross lipid tails and land in clear water.
Checked on one ring vesicle: half the formaldehyde leaks out in 20k steps;
C5 sugars stay, ions almost all stay; sealed rings keep everything.

```bash
python experiments/e05_wetdry.py --headless --chem --leak --meteor-sugar 10 --meteor-every 40000 \
    --name rain_meteor --out runs/e08_leak
```

### Findings (200k steps, 8 seeds, wet-dry)

Sealed membranes (nothing crosses):

| pond | present | longest | carbon at end C1 / C2-4 / C5+ |
|---|---|---|---|
| only Ca2+ / Cl- ions | 82% | 35800 | - |
| rain | 70% | 36400 | 280 / 96 / 9 |
| rain + meteorites | 72% | 42200 | 122 / 321 / 46 |
| meteorites only | 75% | 52200 | 0 / 100 / 11 |

Here ions alone kept vesicles open best: when nothing leaks, anything
trapped (salt, formaldehyde) pushes like sugar. Not realistic.

Leaky membranes:

| pond | present | longest | carbon at end | inside long-lived vesicles |
|---|---|---|---|---|
| only ions | 80% | 29400 | - | ~5 ions |
| rain | 72% | 24200 | 346 / 28 / 10 | no big sugar |
| meteorites only | 87% | 66600 | 0 / 186 / 31 | 3 big-sugar C |
| rain + meteorites | 86% | 77600 | 51 / 370 / 65 | 6 big-sugar C, 10 small; bigger |

- Rain alone barely reacts in a pond (4x more dilute than the e06 box).
  Meteorite sugar seeds the cycle and the pond turns its rain into sugar.
- With leaks, formaldehyde no longer props vesicles; sugar-filled ones last
  longest (record 77600 steps) and are bigger.
- With leaks the pond made more big sugar (C5+ 65 vs 46). Guess, not yet
  tested: food leaks into compartments holding sugar + Ca2+, gets built into
  sugar, and the sugar can't leave.
- Ions leak very slowly, so trapped salt still holds vesicles open too.
- Leak rates are hand-set; only their order is meant to be real.

## Vesicle as a reactor

`experiments/e07_reactor.py`: one ring vesicle with 3 C2 sugar seeds and
4 Ca2+ inside, 150 formaldehyde (food) outside, leaky membrane. Food leaks
in, gets joined onto the seeds, and the bigger sugar can't leak out.

```bash
python experiments/e07_reactor.py --seed 0            # reactor
python experiments/e07_reactor.py --sealed --seed 0   # control: nothing crosses
python experiments/e07_reactor.py --ca none           # control: no catalyst
python experiments/e07_reactor.py --ca outside        # control: catalyst outside the bag
```

### Findings (100k steps, 8 seeds; sugar carbons inside, start 6)

| condition | sugar C inside at end |
|---|---|
| reactor (leaky, Ca2+ inside) | 27.6 (every seed grew: 18-39) |
| sealed | 6.5 |
| no catalyst | 9.0 |
| catalyst outside | 10.5 |

- Nothing "wants" food. Hops in and out are equally likely; joining food
  onto trapped sugar removes it from the free food inside, so more random
  hops go in than out (like cells trapping glucose as glucose-6-phosphate).
- Needs both the leaky wall and the catalyst inside.
- Mostly C2-C4 inside (C2 leaks slowly, so some escapes: ~6 C outside).
- Here the outside has no catalyst. In the pond, Ca2+ is everywhere, so
  the vesicle advantage there is smaller.

## Do full vesicles steal lipid? (e08) Can a ring grow? (e10)

Lab result (Chen, Roberts & Szostak 2004): fatty-acid vesicles swollen by
trapped solute take lipid from empty neighbours. `experiments/e08_compete.py`:
two rings (R 12), sugar in A and/or B, optional lipid drip.

```bash
python experiments/e08_compete.py --sugar-a 30 --sugar-b 0 --drip 500 --mix 0.5 --hop-oily 1.8 --seed 0
```

### Findings (200k steps)

- 3-tail walls (8 seeds x 3 setups): no stealing. Sugar rings swell (area
  ~700 vs ~550) but gain no more lipid; gains 0-11, same as noise.
- Half 2-tail "fatty-acid" walls (`--mix 0.5 --hop-oily 1.8`, 5 seeds x 5 setups):
  still no stealing. Without drip every ring keeps exactly 191 lipids: at
  oiliness 1.8 short-tail lipids never leave a wall in 2D, so nothing is
  traded. With drip (one every 500 steps) both rings gain the same
  (at 100k: sugar ring 203.6, empty ring 203.2; ~+20-40 by 200k). Sugar rings
  are bigger (area at 50k ~590-750 vs ~480-650) but don't last longer;
  rings open at random times (24k to never).
- Why no stealing here: in the lab a stretched membrane has lower tension
  energy per added lipid; our 2D rings are floppy and barely stretched, and
  lipids hardly ever leave a wall to be shared.

`experiments/e10_grow.py`: one mixed ring fed by drip, with shear.
The fed ring grows (191 -> ~210-240 lipids at oiliness 1.4-2.0, 256-294 with
`--near 6` drip next to it) but never splits; rings open after ~50-200k
steps (o1.6-1.8 last longest). In 2D a long "tube" is two parallel lines
with nothing pulling them into a neck, so splitting needs 3D (e15).

## Nitrogen: HCN and ammonia (first step towards pigments)

`alife/chem/nitrogen.py`, `experiments/e14_nitrogen.py`: HCN, ammonia and
formaldehyde alone in a small box (like e06). Rules: Strecker (formaldehyde
+ HCN + NH3 -> glycine), glycolonitrile side product (no NH3), HCN chains
(faster with NH3 nearby as a base, slow first dimer), hydrolysis (lone HCN
-> formic acid + NH3). Long HCN chains = dark HCN polymer, a crude pigment.

```bash
python experiments/e14_nitrogen.py --name base --seed 0
python experiments/e14_nitrogen.py --name dilute --box 80 --seed 0
python experiments/e14_nitrogen.py --plot
```

### Findings (100k steps, 3 seeds; start 300 HCN, 100 NH3, 100 formaldehyde)

| condition | dark polymer (HCN in chains 5+) | glycine | glycolonitrile | HCN lost to water |
|---|---|---|---|---|
| base | 114 | 53 | 47 | 48 |
| no NH3 | 32 | 11 | 89 | 69 |
| HCN only | 90 | 0 | 0 | 81 |
| HCN + NH3 | 202 | 0 | 0 | 57 |
| dilute (4x water) | 0 | 19 | 66 | 70 |
| hot (kT 0.7) | 103 | 52 | 48 | 87 |

- Dilute HCN makes no polymer at all: it is lost to water first (one molecule
  is enough to hydrolyse, two must meet to join). Real chemists say HCN had
  to be concentrated (drying ponds, freezing): same result.
- Ammonia picks glycine over glycolonitrile, and speeds the polymer.
- HCN alone starts slowly (lag ~30-40k), then polymerises: hydrolysis
  releases ammonia, which catalyses the polymer. Nobody coded that loop.
- Formaldehyde competes for HCN (polymer 114 with it vs 202 without).
- Hot: faster, but more HCN lost to water.
- Rates are hand-set; only their order is meant to be real.

### Nitrogen in the wet-dry pond (e05 `--nitrogen --leak`, 400k steps, 3 seeds)

Rain brings 40 HCN + 10 NH3 per refill (~800 HCN in total), plus
formaldehyde (`--rain-carb`). The always-wet pond gets the same rain spread
evenly over time.

```bash
python experiments/e05_wetdry.py --headless --pond wetdry --nitrogen --leak --rain-carb 5 --steps 400000 --seed 0
```

| at the end (per seed) | 40 formaldehyde, wet | 40, wet-dry | 5 formaldehyde, wet | 5, wet-dry |
|---|---|---|---|---|
| HCN in dark polymer (chains 5+) | 5-7 | 12-47 | ~204 | 406-466 |
| longest HCN chain | 5-7 | 6-19 | 22-24 | 79-83 |
| glycine | 161-178 | 152-177 | 34-46 | 51-53 |
| glycolonitrile | 423-456 | 344-359 | 49-57 | 24-29 |
| HCN lost to water | 112-125 | 49-71 | ~285 | ~105 |

- Formaldehyde and pigment compete for the same HCN: lots of formaldehyde
  gives amino acid (glycine) and a dead end (glycolonitrile), little gives polymer.
- Drying concentrates HCN, so it joins before water breaks it: wet-dry gives
  ~2x the polymer, 3-4x longer chains and half the loss to water.
- Low formaldehyde + wet-dry: about half of all rained HCN ends as polymer.
- 5-formaldehyde wet-dry seed 2 blew up at step 189k (two beads on one spot);
  up to then it matched the other seeds (dark 175 vs 184/165), so its
  numbers are from 189k. Cause: a membrane hop (`Leak`) put an NH3 bead
  7 units onto dry land; the shore spring dragged it back fast into a packed
  clump and two lipid tails were pushed together. Fix: `Leak.wet` (hops must
  land in water; e05 sets it to the pond). Same random numbers, so runs are
  unchanged until the first hop that would have landed on land; older
  `--leak` pond runs may differ slightly after that point. Rerun with the
  fix (`--name nlow_fix_wetdry`) finished 200k: dark polymer 200, longest 54,
  in line with the other seeds.

## 3D wet-dry with a flood (e13)

`experiments/e13_drywet3d.py`: lipid gas in 3D (15 mM, ~1420 lipids, box 60).
The box shrinks (pond dries, lipids stack into sheets), then grows back.
A flood = rewet to a bigger box than the start (more water).

```bash
python experiments/e13_drywet3d.py --mM 15 --cycles 1 --rewet-box 100 --last-wet 200000 --seed 0
```

### Findings (3 seeds per flood size, 200k steps after rewetting)

| rewet box | runs with a vesicle | formed (steps after rewet) | sealed until end |
|---|---|---|---|
| 60 (no flood) | 0 / 3 | - | - |
| 100 | 2 / 3 | 136k, 150k | yes (97-100% of logs) |
| 140 | 3 / 3 | 139k, 166k, 178k | yes (95-100% of logs) |

- First 3D vesicles that formed on their own. Drying makes one big lipid
  stack (~1410 of 1420 lipids); after a flood it peels off and closes into a
  single bilayer sphere (outer radius ~11, water core ~490), checked by eye
  in three slices (heads on both faces, water inside).
- Without a flood the sheet spans the whole box (periodic) and can't curl up.
- Control (`--no-dry`, same lipids straight into box 100/140, 290k steps,
  3 seeds each): 0 / 6 vesicles. Lipids stay as ~80-160 small clusters,
  largest only 69-94 lipids, from 100k to the end. So the drying matters,
  not just the dilution: only the dried stack is big enough to close.

## 3D grow and split (e15)

Lab result (Zhu & Szostak 2009): fatty-acid vesicles fed more fatty acid grow
wall faster than volume (solutes inside hold the water), stretch into threads,
and gentle shear breaks the threads into daughter vesicles.
`experiments/e15_split3d.py`: one vesicle (R 10, 2225 lipids, half short-tail),
sugar inside and out (0.03 per unit volume, 35 inside), one short-tail lipid
added near it every 200 steps (1000 in 200k), optional shear bands.

```bash
python experiments/e15_split3d.py --short-oily 0.85 --feed 200 --shear 0.1 --seed 0 --video runs/video/e15.mp4
```

### Findings (200k steps, 3 seeds each, short-tail oiliness 1.0 and 0.85 alike)

| setup | wall lipids (start -> end) | water inside (start -> end) | shape (1 = ball) | split? |
|---|---|---|---|---|
| no food | 2225 -> 2200-2225 | ~1650 -> 1420-1500 | 1.1-1.2 | no |
| food | 2225 -> 3000-3075 | ~1650 -> 2570-2900 | 1.2-1.4 | no |
| food + shear 0.1 | 2225 -> 2620-2790 | ~1650 -> 0 at 30k, 590-770 at end (open tube) | 3.2-3.4 (max 4.7) | no |

- Food works: ~830 of 1000 added lipids join the wall (rest loose or small clumps).
- Without shear the vesicle just gets bigger and stays round, keeping its 35 sugar.
  Volume grows as fast as a ball's would (area x1.37 -> volume x1.6), so the
  sugar does not hold the water back. In the lab it does; here 35 sugar beads
  push far too weakly (implicit water: no real osmotic pressure).
- With shear the vesicle is pulled into a long tube within ~20k steps, tears
  (sugar escapes; no sealed water at 30k) and the tube joins itself across the
  periodic box (box 50 = tube length). The "2 vesicles" logs (~19-27 per run)
  are checked by eye: narrow spots along one tube, not separate daughters;
  the wall stays one cluster.
- Next tries: much stronger osmolyte (or explicit water) so feeding makes
  excess wall; a box long in the shear direction so a tube can't wrap.
