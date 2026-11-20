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
- **Charged heads (2026-11-08 to 11-10):** real fatty acids are half charged at
  pond pH, so we charged half the heads (3 seeds vs 3 neutral, 400k each). Charged
  daughters keep a thin water layer between their walls; neutral ones press wall
  to wall. They still came back together.
- **Plot twist: the fusion was our box.** Tracking the two daughters showed
  they never drifted back together. The current runs in bands: one daughter
  rides the band going left, the other the band going right, and they keep
  separating. But our world is a periodic box (leave one side, come in the
  other, like Pac-Man), so after half a box they meet again from behind. In a
  real pond they would just drift apart. Good on-camera lesson: check if a
  result comes from nature or from the simulation's walls. (`drift.png` shows it:
  every curve climbs to 40 = half the box, then falls.)
- **So:** a gentle current divides a fed vesicle and carries the daughters
  apart, each keeping some of the parent's sugar. Opening and resealing again
  and again loses contents (up to half the sugar by 400k).
- **Long box (2026-11-14):** world twice as long so the daughters can't come
  round the back. They split and pull apart for 200k steps, but stay tied by a
  thin oily thread (a few lipids thick) that stretches and never snaps, in all 6
  runs. Real vesicles in the lab do the same; a shake breaks the thread. So:
  division almost, a "pair on a string". Footage `e15e_long_q0.5_s0.mp4`,
  `runs/e15e/neck.png` (distance up, thread thins to a floor), `cut_s2_300k.png`.
- **Footage:** `e15_o0.85_f200_sh0.1_s0.mp4` (tears), `e15b_*.mp4`,
  `e15c_box80_osm0.15_x5_sh0.02_s0.mp4`, `e15c_box50_osm0.3_x5_s0.mp4`,
  **`e15c_sh0.05_s1_400k.mp4` (split + "fusion")**, **`e15d_q0.5_s0.mp4`
  (charged heads: split, drift apart, meet across the box edge, split again)**,
  pictures `runs/e15c/long/check.png`, `runs/e15d/check200k.png` (charged vs
  neutral walls), `runs/e15d/drift.png` (the plot twist).

## Ep 12. A catalyst inside: can contents decide who grows? (2026-11-10/11)

- **Question:** so far every bag is the same; what's inside doesn't matter. For
  evolution it must. Lab result to test (Budin & Szostak 2011, Adamala & Szostak
  2013): a catalyst trapped in a vesicle joins fatty acids into two-tail lipids
  (like real cell membranes); those never leave, so that vesicle keeps its wall
  and grows while its neighbour shrinks.
- **Built:** two identical bags side by side; A holds 20 catalyst beads, B holds 20
  plain sugar beads (same fullness, no chemistry). Nothing says "grow".
- **Shortcut (said out loud):** the joining costs no energy yet. Next: a fuel
  molecule that gets used up.
- **First try (2026-11-11): the catalyst bag killed itself.** In all 6 runs the
  catalyst bag shrank, flattened and opened: 0% of its sugar and catalyst still
  inside at the end, while the plain bag kept its contents. Was it chemistry? No,
  the model: each join pulled two lipid heads to 0.89 apart while ordinary neighbours
  sit 1.15 apart, so every join shrank the inner wall a bit. Lipids from the outer
  layer flipped in to fill the gap, got joined too, and the inner wall ended up
  crowded and buckled inward. A real two-chain lipid takes about the room of its two fatty
  acids. Also, long-tail lipids (already "finished" membrane lipids in our model)
  were being joined into unrealistic 4-chain ones. Fixed both (joined pair keeps
  normal spacing; only fatty acids join), checked old runs still replay exactly,
  rerunning. Good on-camera moment: the frame where the green catalyst beads
  drift out into the open water (`e16_cats20_s0.mp4`, ~80k steps).
- **Second try: healthy, but no winner.** With the fix the catalyst bag stays round and
  sealed in all 6 runs and makes ~350-420 two-chain lipids. They gather into purple patches
  across the wall. But when both bags are fed, the catalyst bag ends up slightly *behind*
  (3 of 3 seeds), while the plain-vs-plain controls go either way (`fix_result.png`).
- **Why (the user's point):** nothing pulls lipid from one bag to the other. In the lab
  it's bookkeeping: fatty acids hop on and off walls at random, and two-chain lipids make a
  wall hold its fatty acids tighter, so that wall loses fewer than it gains. We counted
  hops: ours are rare (1-2% per 1000 steps; less oily fatty acids hop up to 10x more), and
  the two-chain lipids don't hold neighbours tighter. The extra hopping in the catalyst
  bag is from its inner layer into its own inside water (not lost). The real effect comes
  from two-chain lipids making the wall more ordered; our simple beads don't have that.
  Good on-camera moment: an honest null result, and finding which piece of physics is
  missing instead of forcing the answer.
- **Random oiliness (the user's idea, 2026-11-13):** in nature fatty acids aren't all one
  length; the chemistry that makes them gives a random mix, short ones commoner. Short ones are
  less oily. So we gave each new fatty acid 1, 2 or 3 tail beads at random. Within a few
  thousand steps the walls sorted themselves: short ones fell out into the water (300-600 loose,
  vs ~11 before), long ones stayed, so the walls end up with longer tails than the mix they
  started from. Nobody told them to choose. Hopping went up 10-30x, which is what the lab
  mechanism needs. With mostly short ones (50%) the plain bag opened for a while and resealed.
- **Full test with random oiliness (2026-11-14, 12 runs, 200k):** both bags stay sealed,
  the catalyst bag no longer falls behind, but it doesn't win either: its lead (-33 to +35)
  sits inside the spread of the "sugar in both" controls (-98 to +16). The catalyst turns
  over half its fatty acids into two-chain lipids, yet the rest leave its wall just as
  often as from the plain bag. In real cells, two-chain lipids stiffen the wall and grip
  their neighbours; our beads don't. Honest result: the trick needs that missing physics.
  Footage `e16mix_cats20_feed_s0.mp4`, plot `runs/e16/mix_result.png`.
- **Acid and not acid (the user's question, 2026-11-15):** soap and fatty acids switch
  between a neutral form (H+ on) and a charged one, with pH. Built it: every head can now
  gain or lose its H+, decided by energy, and a neutral head next to a charged one sticks a
  little (they share the H+). Test at 3 pH values, pairing off/on (6 runs, 30k): all bags
  sealed. Pairing keeps a third fewer fatty acids floating loose, but they still hop out as
  often, just come back twice as fast. No "best pH" for holding on. Caveat: our walls go
  half-charged at pH 5, real ones at 7-9 (our water is too simple). Plot
  `runs/e17/e17_result.png`.
- **Viewer ideas for later (from the user):** "ports" in the wall that light or heat
  open so charge can flow (real: light-, heat-, stretch-gated channels), and
  reactions that give off heat or gas and change what happens next. On the
  roadmap: fuel -> light-driven charge pump -> gated channels.
- **Footage:** `e16_cats20_s0.mp4` (first try: the bag collapses and leaks);
  `e16fix_cats20_s0.mp4` (fixed: stays sealed, purple patches grow), `runs/e16/fix_result.png`
  (wall A minus wall B, catalyst vs control).
- **Sudden pH swings, fuel and fatty alcohol (the user's idea, 2026-11-16):** ponds get
  sudden acid (volcanic gas, acid rain) or alkaline (spring water) events, like heat waves.
  Built: pH shocks that fade back; bags that keep their own H+ inside; fuel that a catalyst
  burns, giving off acid inside the bag; fatty alcohol (a never-charged cousin). Result: a
  shock to pH 2 or 8 never broke a bag. The outside of the wall changed at once, the inside a
  few thousand steps later: for a while the bag held a pH difference, then it leaked through.
  The fuel made only a tiny, short acid pulse inside, because the acid leaks out as fast as
  it's made. Fatty alcohol steadied the wall a little. Making the neutral heads a bit oily
  squeezed the water out of the bag into a thick oily shell (the lab sees vesicles turn to oil
  drops in acid). Footage idea: the shock plot (outside line jumps, inside line follows late).
  `runs/e18/e18_result.png`, `runs/e18/e18_slice.png`. Also the earlier pH test, rerun with
  3 seeds: same answer (`runs/e17/e17_result.png`).

## Ep 13. Light: a pump in the wall, and pigment made by the pond (2026-11-16)

- **More shocks (the user's ideas):** heat from burning fuel, a pH held at 2 or 8 for a long
  time, random shocks, rain that doubles the water then dries back. 18 runs, 3 seeds each:
  no bag broke. Long acid or alkaline spells shake loose more fatty acids, that's all.
- **A light-driven pump (pigment put in by hand):** 30 pigment heads in bag A's wall push H+
  into A when light hits them. Day and night (10k steps each): in the day A's inside goes acid,
  at night it leaks back within ~2000 steps, next day again. 6/6 runs. First time a bag keeps
  a real pH difference with its pond: light beats the leak that fuel (Ep 12) couldn't. It
  costs A some wall lipid. Plot `runs/e19/e19_pump.png` (the "day/night" wave is a nice shot).
- **"Did you code the pigment?" (the user's question):** yes, in that test. So next the pond
  made its own: HCN polymer, the brown-black tar from Ep 10. Saved a bag world and a puddle
  with lots of HCN, then poured them together like two waters meeting (new tool: save, carry
  on, join worlds). The puddle kept making tar; it drifted to the bags and touched the walls.
  In light it warms its neighbours (and, a guess, makes a little acid). Result: the bags
  didn't care. Heat spreads out at once; the acid lands outside. The pump only worked because
  it sat in the wall, pointing in.
- **The user's rain idea:** maybe life started in rainy times, so HCN never dried out. Tried
  it: the same HCN spread through all the water made ~50x less tar. Rain helps in other ways
  (fresh molecules, wet-dry cycles, Ep 4/5), but pigment needs HCN concentrated by drying or
  freezing, as chemists think.
- Next idea: real HCN tar is oily/insoluble, so it might sink into walls by itself (not done).
  Plot `runs/e20/e20_glow.png`.
- **Tar that sinks (2026-11-17, the user: "make the pigment, not a shortcut"):** real HCN tar
  doesn't dissolve, so long chains became oily "tar". Nothing placed it. It sank into the bag
  walls by itself, more into one bag than the other, by chance. Then in daylight that bag's inside
  went acid, and at night it leaked back, in all 6 lit runs; in the dark, nothing. No direction was
  coded: acid made on the outside of the wall washes into the pond, acid made on the inside is
  trapped. The first pump-like behaviour that nobody built. Weaker than the hand-made pump (about a
  third). Guesses: how oily the tar is, and that light on it makes acid. Plot `runs/e21/e21_tar.png`
  (the day/night dips are the shot for the episode).

## Ep 14. The habitable zone: our pond on other planets (2026-11-17/18)

- **The user's question:** Earth is at just the right distance from the sun; a bit closer or further
  and everything would have died. Does our pond agree? (Gentle correction: the sun doesn't really "go
  round" us, the planet spins; and Venus and Mars got their fates mostly from their air, which we don't
  model.)
- **Built:** a smooth sun. Light rises and falls through the day, the water warms a little after it.
  A closer planet gets more light (1/distance squared) and is warmer; chemistry runs faster when warm.
  Same pond (HCN puddle + two bags) on 5 planets, 3 seeds each.
- **What happened:** the bags held on every planet. Hot planets made ~4x the tar of cold ones, and more
  of it sank into walls. The tar pump (Ep 13) worked on Earth and a bit closer; further out it nearly
  vanished (too little tar, too little light). So the "zone" for our crude pump is Earth-ish.
- **Too close crashed:** at Venus distance the simulation broke on the first bright day. The tar
  soaked up so much light that its own heat outran the cooling, and the beads jumped. Checked: without
  photon heat the same world is fine. That's our model's limit, not the bags dying (real heat from one
  photon leaves in a trillionth of a second). Honest reading: we can't say yet what happens at Venus.
- Caveats: no air, no ice or boiling, bags moved from Earth, day hugely sped up.
  Plot `runs/e22/e22_zone.png` (four panels: tar, tar in walls, loose lipids, pump vs distance).
- **Fix attempt (2026-11-19/20):** we made the photon heat gentler: instead of a kick, the light warms
  the water around the tar for a short while. All 15 rerun. The same 4 hot runs still crashed: being
  gentle didn't matter, the light simply puts too much heat into the tar (it got ~5x hotter than the
  water). Lesson for the video: the bug was the *amount* of a guessed number, not the method.
- **Second look at the pump:** the tar amounts came out the same as the first run (robust), but the
  pump's best distance moved: now Earth, 1.2 and Mars-ish pump about equally. With 3 seeds per planet the
  pump is too noisy to say "Earth is best". Honest: tar-rich bags go acid by day wherever we could run.
  Plot `runs/e22/e22w_zone.png` (grey dashed = first run).
