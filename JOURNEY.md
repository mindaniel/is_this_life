# Journey log (for the YouTube channel)

The story of the project in order, written to be turned into videos later.
Each entry: the question, what we built, what happened (surprises and
failures included: they make the best stories), the honest caveats, and the
footage. Numbers and details live in README.md; this file is the narrative.

Footage is in `runs/video/` and `runs/**/` pictures. **`runs/` is not in git**
(~370 MB of video as of 2026-11-08): back it up separately before cleaning
the server. Videos are filmed at fps 6 unless noted.

Running themes to come back to on camera:
- Nothing is told what to do. Only physics rules (oil hides from water,
  charges, heat jiggling) and chemistry rules (which bonds form), then we watch.
- Shortcuts are allowed as steps, but we say them out loud (implicit water,
  ready-made start vesicle, ...).
- Same seed = same history: we can replay any run exactly and film it.
- Failures are results: many "vesicles" turned out to be fake when checked by eye.

---

## Ep 1. Oil and water make a skin (2026-10-21)

- **Question:** can membranes build themselves from nothing but "oily tails
  hide from water"?
- **Built:** 2D lipids (a head that likes water + oily tail beads on springs),
  heat noise. No rule says "make a membrane".
- **Happened:** too hot: nothing sticks. Too cold: glassy blobs. In between
  (kT 0.6, stiff 3-tail lipids): micelles and bilayer strips appear by themselves.
  Shape matters: floppy lipids pack badly.
- **Footage:** e01 live window / PNG snapshots (no video yet; easy to film).

## Ep 2. A planet is never calm (2026-10-21)

- **Question:** do seasons, day/night, currents and asteroid hits help?
- **Happened:** a calm pond packs tightest, but only the full weather (seasons +
  sun + currents + asteroids) ever made vesicles: 4 times, each broke within
  ~1000 steps. "Disturbance makes and breaks."
- **Then charge:** a periodic table of beads (Na, K, Cl, Mg, Ca, fatty acid).
  Calcium grabs fatty-acid heads and glues lipids; magnesium doesn't (its thick
  water shell keeps it away). Real chemistry fell out of size + charge alone.

## Ep 3. A big world with hot vents and pools (2026-10-21)

- **Built:** 3x3 zones, hot vent, cold corner, rock pools, local asteroids.
- **Happened:** "vesicles" detected, but checked by eye they are tiny water
  holes inside clumps, living ~1000 steps. Lesson on camera: always check
  your detector (it fooled us twice).

## Ep 4. Ponds that dry and refill (2026-10-22)

- **Question:** do wet-dry cycles (a real early-Earth idea) help?
- **Happened:** drying squeezes lipids into one sheet; refilling leaves it
  folded with water sealed inside: **25x more vesicles** than an always-wet pond.
  Adding sugar that can't cross the wall holds them open (osmosis): present
  90% of the time, record life 65000 steps. Nobody put sugar "inside": drying
  trapped it.
- **Footage:** `rain_meteor_seed7.mp4` (+ `_slow`), `ponds_test.mp4`.

## Ep 5. Making sugar from formaldehyde (2026-10-22)

- **Built:** formose chemistry: formaldehyde (rains from the sky) joins into
  sugars with calcium as catalyst; long chains split.
- **Happened:** a long quiet wait, then sugar explodes (autocatalysis, like the
  real reaction). A few small sugars at the start skip the wait. No calcium or
  cold: nothing. Meteorite sugar can kick-start a whole pond.

## Ep 6. Leaky bags and the first "metabolism" (2026-10-22/23)

- **Built:** membranes that leak by size and charge (formaldehyde fast, ions
  slow, big sugar never).
- **Happened:** sugar-filled vesicles live longest (77600 steps). The vesicle as
  a reactor: food leaks in, gets built onto trapped sugar, can't leave. Sugar
  inside grew 6 -> 27.6 carbons (every seed), vs 6-10 in all controls. Nothing
  "wants" food; trapping it just makes more hops go in than out (like cells
  trapping glucose).
- **Footage:** `reactor_seed4.mp4`.

## Ep 7. Do full bubbles steal lipid? (2026-10-23, 2026-10-30)

- **Question:** (Szostak's lab result) do sugar-filled vesicles grow by taking
  lipid from empty ones?
- **Happened:** no sign, with normal or half fatty-acid walls. Fed rings grow
  (to ~250 lipids) but never split in 2D. A 2D "tube" can't pinch: we need 3D.
- **Footage:** `compete_seed0.mp4`, `e08_fatty_a30_b0_s0.mp4`, `grow_mix50.mp4`,
  `grow_sticky18.mp4`.

## Ep 8. Weather, many ponds, a river basin (2026-10-23 to 10-27)

- **Built:** 4 seasons, wind gyres, storms; many ponds that overflow into each
  other; shore crust; a landscape with streams.
- **Surprise:** we expected lipids to wash down into the big lake. The opposite:
  the lake lost lipid and it piled into one or two small ponds (wind + floods
  beat the stream). A landscape model in real units says lipid reaches vesicle-
  forming concentration only a few hours a year, in hot ponds nearly dry.
- **Footage:** `weather_year_seed7.mp4`, `landscape_year_seed0.mp4`,
  `landscape_land1_seed0.mp4`, `basin_year_seed0.mp4`, `e11_basin.mp4`,
  `e11_basin_zoom.mp4`, `basin_maps.png`.

## Ep 9. Going 3D (2026-10-27/30)

- **Happened:** a ready-made 3D vesicle stays sealed. Flat discs don't close on
  their own, and loose lipid only makes patches. But **dry it and flood it**: the
  dried lipid stack closes into one big 3D vesicle in 5/6 runs. Control with no
  drying: 0/6. Drying is what does it.
- **Footage:** `e12_vesicle_close.mp4`, `e12_gas.mp4`, `e12_disc18.mp4`,
  `e13_drywet.mp4`, `e13_flood140_s0.mp4`, `e13_nodry140_s0.mp4` (control).

## Ep 10. Nitrogen: towards pigments (2026-10-27/30)

- **Built:** HCN and ammonia; glycine (an amino acid) forms; HCN polymer (a dark,
  light-absorbing stuff) forms.
- **Happened:** formaldehyde eats the HCN; with more formaldehyde rain, HCN
  polymerises, and wet-dry makes ~2x the polymer of always-wet. A crash taught a
  rule: anything that jumps must land in water.
- **Footage:** `e14_base.mp4`, `e05_nitrogen_wetdry_s0.mp4`,
  `e05_nitrogen_lowcarb_wetdry_s0.mp4`.

## Ep 11. Can a cell divide without being told to? (2026-10-31 to 11-08)

- **Question:** (Zhu & Szostak 2009) fed fatty-acid vesicles grow wall faster
  than volume, go long and thread-like, and a gentle shake splits them.
- **Try 1:** feeding grows the wall but it stays a ball. Strong current tears it.
- **Try 2:** more sugar inside: it goes long and stays sealed, but no split.
- **Try 3:** much more food + a *gentle* current: **the vesicle pinches into two
  sealed sacs, each keeping sugar from the parent, in 3/3 seeds.** Nothing cuts
  it. Without current, food piles into lumps that make rooms in one body: a fake
  split we caught by slicing it open.
- **Longer film (400k steps):** the two daughters stay pressed together and
  **fuse back** at ~225k, then the vesicle stretches and splits again. Division
  happens but isn't final.
- **Next (running 2026-11-08):** charged fatty-acid heads (real fatty acids are
  half charged at pond pH): do the daughters repel and stay apart?
- **Footage:** `e15_o0.85_f200_sh0.1_s0.mp4` (tears), `e15b_*.mp4`,
  `e15c_box80_osm0.15_x5_sh0.02_s0.mp4`, `e15c_box50_osm0.3_x5_s0.mp4`,
  **`e15c_sh0.05_s1_400k.mp4` (split + fusion, the key film)**, pictures
  `runs/e15c/long/check.png` (cut views).
