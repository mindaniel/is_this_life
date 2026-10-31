"""Environment: the planet around the chemistry.

Three kinds of change act on a World:

- seasons: Earth orbits the sun. Temperature swings slowly up and down.
- day/night: Earth spins. Light turns on and off, fast. Sunlight warms
  matter (day_heat). `light` is kept as its own signal so pigments can
  absorb it later; for plain lipids light is only heat.
- currents: the sun heats unevenly, warm water rises, cool water sinks, and
  the water turns in swirls (eddies). The swirl pattern drifts like weather.
  Stronger by day, calmer at night. At our tiny scale this is how wind and
  waves show up: stirring.
- shocks: rare random events (asteroid impact). A shock heats everything
  (heat then fades away) and blasts beads near the impact point outward,
  mixing matter that was stuck in place.

Nothing here knows about membranes. It only changes kT, light, water flow,
and velocities.
"""
import numpy as np


class Environment:
    def __init__(self, base_kT=0.6, amp=0.0, period=20000,
                 day_period=0, day_heat=0.0,
                 current=0.0, eddies=2, drift_period=50000, night_current=0.3,
                 shock_every=0, shock_heat=0.5, shock_fade=3000,
                 blast_radius=15.0, blast_speed=3.0, calm_after=None, seed=0,
                 terrain=None, shock_local=False):
        """
        base_kT       average temperature
        amp           season size: kT swings between base-amp and base+amp
        period        steps per season cycle (one 'year')
        day_period    steps per day/night (0 = no sun cycle, constant light 0)
        day_heat      extra kT at noon from sunlight
        current       top water speed of the swirls (0 = still water)
        eddies        swirl pairs across the box
        drift_period  steps for the swirl pattern to drift one box length
        night_current fraction of the current left at night (with a day cycle)
        shock_every   mean steps between shocks (0 = never)
        shock_heat    extra kT right after a shock
        shock_fade    steps for the extra heat to fade to 37%
        blast_radius  beads closer than this to impact get kicked outward
        blast_speed   kick speed at impact point (falls to 0 at radius)
        calm_after    step after which cycle and shocks stop (fair comparison)
        terrain       Terrain: places that are warmer/colder, sunny/shaded
        shock_local   True: shock heat stays near the impact (big world).
                      False: the whole box heats (one small droplet)
        """
        self.base_kT, self.amp, self.period = base_kT, amp, period
        self.shock_every, self.shock_heat, self.shock_fade = shock_every, shock_heat, shock_fade
        self.blast_radius, self.blast_speed = blast_radius, blast_speed
        self.day_period, self.day_heat = day_period, day_heat
        self.current, self.eddies = current, eddies
        self.drift_period, self.night_current = drift_period, night_current
        self.light = 0.0
        self.calm_after = calm_after
        self.rng = np.random.default_rng(seed)
        self.extra_heat = 0.0
        self.terrain, self.shock_local = terrain, shock_local
        self.hot_spots = []   # [x, y, heat] of recent local impacts
        self.events = []   # (step, x, y)

    def calm(self, t):
        return self.calm_after is not None and t >= self.calm_after

    def kT_at(self, t):
        if self.calm(t):
            return self.base_kT + self.extra_heat
        season = self.amp * np.sin(2 * np.pi * t / self.period)
        sun = self.day_heat * self.light_at(t)
        if self.terrain is not None:
            sun *= self.terrain.mean_sun
        return self.base_kT + season + sun + self.extra_heat

    def light_at(self, t):
        """Sunlight 0..1: rises at dawn, peaks at noon, 0 all night."""
        if not self.day_period or self.calm(t):
            return 0.0
        return max(0.0, np.sin(2 * np.pi * t / self.day_period))

    def flow_at(self, t, L):
        """Water velocity field at time t: drifting swirls, or None if still.

        Built from a stream function, so water is never created or destroyed
        (what flows out of one spot flows into the next).
        """
        if not self.current or self.calm(t):
            return None
        strength = self.current
        if self.day_period:
            strength *= self.night_current + (1 - self.night_current) * self.light_at(t)
        k = 2 * np.pi * self.eddies / L
        px = 2 * np.pi * t / self.drift_period            # weather moves east...
        py = 2 * np.pi * t / (1.7 * self.drift_period)    # ...and slowly north

        def flow(pos):
            sx, cx = np.sin(k * pos[:, 0] + px), np.cos(k * pos[:, 0] + px)
            sy, cy = np.sin(k * pos[:, 1] + py), np.cos(k * pos[:, 1] + py)
            return strength * np.stack([sx * cy, -cx * sy], axis=1)

        return flow

    def shock(self, world, x=None, y=None):
        """Asteroid: heat spike + outward blast around (x, y)."""
        x = self.rng.uniform(0, world.L) if x is None else x
        y = self.rng.uniform(0, world.L) if y is None else y
        d = world._wrapped() - np.array([x, y])
        d -= world.L * np.round(d / world.L)
        dist = np.sqrt(np.einsum("ij,ij->i", d, d))
        hit = dist < self.blast_radius
        push = self.blast_speed * (1 - dist[hit] / self.blast_radius)
        world.vel[hit] += (push / np.maximum(dist[hit], 1e-9))[:, None] * d[hit]
        if self.shock_local:
            self.hot_spots.append([x, y, self.shock_heat])
        else:
            self.extra_heat = self.shock_heat
        self.events.append((world.t, x, y))

    def kT_field_at(self, t, L):
        """Temperature map: function pos -> kT per bead, or None if uniform."""
        if self.terrain is None and not self.hot_spots:
            return None
        mean = self.kT_at(t)
        sun = self.day_heat * self.light_at(t)
        terrain, spots, R = self.terrain, [s[:] for s in self.hot_spots], self.blast_radius

        def field(pos):
            kT = np.full(len(pos), mean)
            if terrain is not None:
                # kT_at already added the box-average sun; swap in the local amount
                kT += terrain.heat(pos) + sun * (terrain.sun(pos) - terrain.mean_sun)
            for x, y, heat in spots:
                d = pos - np.array([x, y])
                d -= L * np.round(d / L)
                kT += heat * np.exp(-np.einsum("ij,ij->i", d, d) / R ** 2)
            return kT

        return field

    def run(self, world, n, chunk=100):
        """Advance world n steps, updating temperature every chunk steps."""
        done = 0
        while done < n:
            k = min(chunk, n - done)
            if self.shock_every and not self.calm(world.t):
                if self.rng.random() < k / self.shock_every:
                    self.shock(world)
            self.light = self.light_at(world.t)
            world.kT = self.kT_at(world.t)
            world.flow = self.flow_at(world.t, world.L)
            world.kT_field = self.kT_field_at(world.t, world.L)
            world.step(k)
            fade = np.exp(-k / self.shock_fade)
            self.extra_heat *= fade
            for s in self.hot_spots:
                s[2] *= fade
            self.hot_spots = [s for s in self.hot_spots if s[2] > 0.01]
            done += k


class Terrain:
    """A map of places: each zone of an n x n grid has its own warmth and sun.

    heat[i][j]  kT added in zone (row i = y, column j = x). Hot vent: +0.2
    sun[i][j]   0..1 share of sunlight reaching the zone (1 open, 0.2 shaded)
    Values blend smoothly between zone centers, like real places do.
    """

    def __init__(self, L, heat, sun=None):
        self.L = float(L)
        self.heat_grid = np.asarray(heat, dtype=float)
        self.sun_grid = np.ones_like(self.heat_grid) if sun is None else np.asarray(sun, dtype=float)
        self.mean_sun = float(self.sun_grid.mean())

    def _blend(self, grid, pos):
        """Bilinear blend between zone centers, wrapping around the box."""
        n = grid.shape[0]
        u = (pos / self.L) * n - 0.5          # zone coordinates, centers at integers
        i0 = np.floor(u).astype(int)
        f = u - i0
        x0, y0 = i0[:, 0] % n, i0[:, 1] % n
        x1, y1 = (x0 + 1) % n, (y0 + 1) % n
        fx, fy = f[:, 0], f[:, 1]
        return ((1 - fx) * (1 - fy) * grid[y0, x0] + fx * (1 - fy) * grid[y0, x1]
                + (1 - fx) * fy * grid[y1, x0] + fx * fy * grid[y1, x1])

    def heat(self, pos):
        return self._blend(self.heat_grid, np.mod(pos, self.L))

    def sun(self, pos):
        return self._blend(self.sun_grid, np.mod(pos, self.L))


class Pond:
    """Round ponds that dry up and refill: wet-dry cycles.

    Water covers disks of radius R(t) around `centers` (default: one pond in
    the middle of the box). Outside is dry land. Sun evaporates a pond (R
    shrinks to r_dry), rain refills it (R grows back to r_wet), once per
    `period` steps. The run starts wet.

    Beads on dry land feel a spring pushing them back into the nearest pond:
    the retreating edge sweeps lipids together, so they get crowded while
    the pond is small and spread out again after rain. A bead inside any
    pond is wet. Neighbouring ponds that are both nearly full overlap: the
    water joins (a flood) and molecules can drift across, until drying
    splits them again. Each pond keeps its own timing (`phases` as a
    fraction of a cycle, `periods` per pond), so ponds are not in step.
    period = 0 keeps the ponds always wet, for comparison.

    rain > 0: refilling stirs the water. Rain drops and inflowing water make
    swirls whose speed follows how fast the pond is refilling (top speed
    `rain` at the fastest refill, nothing while it dries or stays full).

    weather: optional Weather (seasons, steady wind, storms). Wind drives a
    big circular current in each pond; summer dries ponds deepest, winter
    hardly at all. Without it the pond behaves exactly as before.

    Ponds need not be alike or round (all three default to plain circles):
      scales  per pond size: r_wet and the drying range times this
      depths  per pond: how much of the drying range it goes through (a
              shallow pond dries further, 1 = down to r_dry)
      shapes  per pond, a list of (m, amp, phase) wobbles of the shoreline:
              R(angle) = R (1 + sum amp cos(m angle + phase)). m=2 makes it
              oval, m=3 three-lobed, bigger m a ragged shore.

    streams: optional list of (upper, lower) pond pairs. When the upper pond
    is nearly full (its radius above `spill` of its fullest) it overflows
    and a stream `stream_width` wide runs downhill into the lower pond,
    carrying water and whatever floats in it (top speed `stream_speed`). In
    a dry spell the stream bed dries out and is land again.

    shore_max: cap on the shore spring's pull (None = no cap). Something
    left far out on land (a dried stream bed, a storm) is otherwise yanked
    back so hard that it smashes into frozen crust beads.
    """

    def __init__(self, L, r_wet, r_dry, period, k=5.0, rain=0.0, eddies=3, seed=0,
                 centers=None, phases=None, periods=None, weather=None, scales=None, depths=None, shapes=None,
                 streams=None, stream_width=6.0, stream_speed=1.0, spill=0.9, shore_max=None):
        self.weather = weather
        self.scales = None if scales is None else np.asarray(scales, dtype=float)
        self.depths = None if depths is None else np.asarray(depths, dtype=float)
        self.shapes = shapes
        self.streams = list(streams or [])
        self.stream_width, self.stream_speed, self.spill = stream_width, stream_speed, spill
        self.shore_max = shore_max
        self.L, self.r_wet, self.r_dry, self.period, self.k = float(L), r_wet, r_dry, period, k
        self.rain, self.eddies = rain, eddies
        self.rng = np.random.default_rng(seed)
        self.centers = np.array([[L / 2, L / 2]] if centers is None else centers, dtype=float)
        n = len(self.centers)
        self.phases = np.zeros(n) if phases is None else np.asarray(phases, dtype=float)
        self.periods = np.full(n, float(period)) if periods is None else np.asarray(periods, dtype=float)

    def _angle(self, t):
        # same arithmetic as the one-pond version, so old runs reproduce exactly
        return 2 * np.pi * t / np.maximum(self.periods, 1) + 2 * np.pi * self.phases

    def radii(self, t):
        """Radius of each pond at time t (its mean radius, if not round)."""
        R = self._radii(t)
        if self.depths is not None:
            R = self.r_wet - (self.r_wet - R) * self.depths
        if self.scales is not None:
            R = R * self.scales
        return R

    def _radii(self, t):
        if not self.period:
            return np.full(len(self.centers), float(self.r_wet))
        wet = 0.5 + 0.5 * np.cos(self._angle(t))   # 1 = full, 0 = driest
        if self.weather is not None:   # how deep this season's drying goes
            return self.r_wet - (self.r_wet - self.r_dry) * self.weather.season(t)["dry"] * (1 - wet)
        return self.r_dry + (self.r_wet - self.r_dry) * wet

    def wobble(self, angle):
        """Shoreline factor g and dg/dangle, shape (positions, ponds); angle has that shape too."""
        g, gp = np.ones_like(angle), np.zeros_like(angle)
        for i, waves in enumerate(self.shapes or []):
            for m, amp, phase in waves:
                g[:, i] += amp * np.cos(m * angle[:, i] + phase)
                gp[:, i] -= amp * m * np.sin(m * angle[:, i] + phase)
        return g, gp

    def _waves(self):
        """Shoreline waves as (ponds, waves) arrays m, amp, phase (amp 0 = no wave), and sum |amp|."""
        if getattr(self, "_wave_table", None) is None:
            n = max(len(w) for w in self.shapes)
            m, amp, ph = np.zeros((3, len(self.shapes), n))
            for i, waves in enumerate(self.shapes):
                for w, (mm, a, p) in enumerate(waves):
                    m[i, w], amp[i, w], ph[i, w] = mm, a, p
            self._wave_table = m, amp, ph, np.abs(amp).sum(1)
        return self._wave_table

    def gaps(self, d, r, R):
        """r - R g(angle) per (position, pond), with g and dg/dangle, like wobble() on every pair.

        Only pairs that matter are worked out: a position inside some pond's
        widest shoreline needs just the ponds it may be in (the rest are on
        dry land for it, so neither nearest nor wet); a position on dry land
        needs them all. Skipped pairs get gap inf, g 1, dg 0. R: radius per pond."""
        m, amp, ph, most = self._waves()
        gap, g, gp = np.full(r.shape, np.inf), np.ones(r.shape), np.zeros(r.shape)

        def fill(mask):
            i, j = np.nonzero(mask)
            a = np.arctan2(d[i, j, 1], d[i, j, 0])
            gg, ggp = np.ones(len(i)), np.zeros(len(i))
            for w in range(m.shape[1]):   # same sums, in the same order, as wobble()
                mm, aa, pp = m[j, w], amp[j, w], ph[j, w]
                gg += aa * np.cos(mm * a + pp)
                ggp -= aa * mm * np.sin(mm * a + pp)
            g[i, j], gp[i, j] = gg, ggp
            gap[i, j] = r[i, j] - R[j] * gg

        fill(r < (R * (1 + most))[None])
        fill((gap.min(1) > 0)[:, None] & np.isinf(gap))
        return gap, g, gp

    def geom(self, pos, t):
        """where() plus the nearest pond's radius R, wobble g and dg/dangle at each position."""
        d = pos[:, None, :] - self.centers[None]
        d -= self.L * np.round(d / self.L)
        r = np.sqrt(np.einsum("ijk,ijk->ij", d, d))
        R = self.radii(t)[None]
        if self.shapes is None:
            g = gp = None
            gap = r - R
        else:
            gap, g, gp = self.gaps(d, r, R[0])
        near = np.argmin(gap, axis=1)
        idx = np.arange(len(pos))
        out = (near, np.maximum(gap[idx, near], 0.0), d[idx, near], r[idx, near], R[0, near])
        if self.streams:   # a flowing stream bed is wet
            wet = np.zeros(len(pos), dtype=bool)
            for (a, b), flow in zip(self.streams, self.stream_flows(t)):
                if flow > 0:
                    wet |= np.abs(self._along(pos, a, b)[1]) < self.stream_width / 2
            out = (near, np.where(wet, 0.0, out[1])) + out[2:]
        if g is None:
            return out + (np.ones(len(pos)), np.zeros(len(pos)))
        return out + (g[idx, near], gp[idx, near])

    def stream_flows(self, t):
        """Per stream 0..1: how hard the upper pond overflows (0 = dry bed)."""
        full = self.r_wet * (1.0 if self.scales is None else self.scales)
        fill = self.radii(t) / full
        return [float(np.clip((fill[a] - self.spill) / (1 - self.spill), 0, 1)) for a, _ in self.streams]

    def _along(self, pos, a, b):
        """Position along the line centre a -> centre b (clipped to it) and distance across it."""
        ab = self.centers[b] - self.centers[a]
        ab -= self.L * np.round(ab / self.L)
        length = np.hypot(*ab)
        e = ab / length
        d = pos - self.centers[a]
        d -= self.L * np.round(d / self.L)
        u = d @ e
        across = d @ np.array([-e[1], e[0]])
        outside = (u < 0) | (u > length)
        across = np.where(outside, np.inf, across)
        return u, across, e, length

    def stream_flow_at(self, t):
        """Water velocity in the streams: pos -> velocity (0 away from them)."""
        flows = self.stream_flows(t)
        w = self.stream_width / 2

        def flow(pos):
            v = np.zeros_like(pos)
            for (a, b), f in zip(self.streams, flows):
                if f <= 0:
                    continue
                u, across, e, length = self._along(pos, a, b)
                # flat profile across, fastest in the middle; fades in over 8 units at both ends so
                # the stream draws water from the upper pond and spreads out in the lower one
                ramp = np.clip(np.minimum(u, length - u) / 8.0, 0, 1)
                prof = np.clip(1 - (across / w) ** 2, 0, 1)
                v += (self.stream_speed * f * prof * ramp)[:, None] * e
            return v

        return flow

    def stream_outline(self, i):
        """Stream i's bed as a 4-corner polygon (unwrapped coordinates)."""
        a, b = self.streams[i]
        ab = self.centers[b] - self.centers[a]
        ab -= self.L * np.round(ab / self.L)
        e = ab / np.hypot(*ab)
        side = np.array([-e[1], e[0]]) * self.stream_width / 2
        p, q = self.centers[a], self.centers[a] + ab
        return np.array([p + side, q + side, q - side, p - side])

    def outline(self, i, t, n=180):
        """Shoreline of pond i at time t as an (n, 2) polygon."""
        a = np.linspace(0, 2 * np.pi, n, endpoint=False)
        g = self.wobble(np.tile(a[:, None], (1, len(self.centers))))[0][:, i]
        return self.centers[i] + (self.radii(t)[i] * g)[:, None] * np.stack([np.cos(a), np.sin(a)], 1)

    def radius(self, t):
        """Mean pond radius (the radius, for one pond)."""
        return float(self.radii(t).mean())

    def refill_speeds(self, t):
        """Per pond, 0..1: 0 while drying or full, 1 at the fastest refill."""
        if not self.period:
            return np.zeros(len(self.centers))
        return np.maximum(0.0, -np.sin(self._angle(t)))   # R = cos: refills in the 2nd half

    def refill_speed(self, t):
        return float(self.refill_speeds(t).mean())

    def where(self, pos, t):
        """For each position: nearest pond, how far onto dry land (0 = wet), offset from its centre."""
        if self.shapes is not None:
            return self.geom(pos, t)[:4]
        d = pos[:, None, :] - self.centers[None]
        d -= self.L * np.round(d / self.L)
        r = np.sqrt(np.einsum("ijk,ijk->ij", d, d))
        gap = r - self.radii(t)[None]   # < 0 inside that pond
        near = np.argmin(gap, axis=1)
        idx = np.arange(len(pos))
        return near, np.maximum(gap[idx, near], 0.0), d[idx, near], r[idx, near]

    def random_point(self, rng, t, margin=1.0, pond=None):
        """Random spot in the water: in pond `pond`, or a pond picked by area. Returns (point, pond)."""
        R = np.maximum(self.radii(t) - margin, 0.0)
        if pond is None:
            pond = rng.choice(len(R), p=R ** 2 / (R ** 2).sum())
        if self.shapes is not None:   # throw darts at the bounding circle until one lands in water
            for _ in range(1000):
                r, a = 1.6 * R[pond] * np.sqrt(rng.uniform()), rng.uniform(0, 2 * np.pi)
                p = self.centers[pond] + r * np.array([np.cos(a), np.sin(a)])
                near, _, _, _, Rn, g, _ = self.geom(p[None], t)
                if near[0] == pond and r < Rn[0] * g[0] - margin:
                    return p % self.L, pond
        r, a = R[pond] * np.sqrt(rng.uniform()), rng.uniform(0, 2 * np.pi)
        return (self.centers[pond] + r * np.array([np.cos(a), np.sin(a)])) % self.L, pond

    def flow_at(self, t):
        """Rain-stirred water: random swirls (new pattern each call), or None.

        With weather: the wind current, plus rain splash scaled by the season."""
        speeds = self.rain * self.refill_speeds(t)
        if self.weather is not None:
            base = self.weather.flow_at(t, self, speeds * self.weather.season(t)["rain"])
            if not self.streams:
                return base
            streams = self.stream_flow_at(t)
            return lambda pos: base(pos) + streams(pos)
        if speeds.max() <= 0:
            return None
        k = 2 * np.pi * self.eddies / self.L
        px, py = self.rng.uniform(0, 2 * np.pi, 2)

        def flow(pos):
            sx, cx = np.sin(k * pos[:, 0] + px), np.cos(k * pos[:, 0] + px)
            sy, cy = np.sin(k * pos[:, 1] + py), np.cos(k * pos[:, 1] + py)
            speed = speeds[0] if len(speeds) == 1 else speeds[self.where(pos, t)[0]][:, None]
            return speed * np.stack([sx * cy, -cx * sy], axis=1)

        return flow

    def shore_at(self, t, frozen=None):
        """External force for World: pos -> (force, energy). Rock and crust are left alone."""
        k = self.k

        def shore(pos):
            _, out, d, r = self.where(pos, t)   # how far onto dry land
            if frozen is not None:
                out[frozen] = 0.0
            if self.shore_max is None:
                F = (-k * out / np.maximum(r, 1e-9))[:, None] * d
                return F, float(np.sum(0.5 * k * out ** 2))
            c = self.shore_max / k   # spring up to c onto land, constant pull beyond
            pull = np.minimum(k * out, self.shore_max)
            F = (-pull / np.maximum(r, 1e-9))[:, None] * d
            E = np.where(out < c, 0.5 * k * out ** 2, 0.5 * k * c * c + self.shore_max * (out - c))
            return F, float(np.sum(E))

        return shore


class Crust:
    """Shore crust: molecules stranded on the dry rim stick there until the water returns.

    When a pond dries, not everything follows the water inward: lipids,
    salt and sugar left on the wet mud dry into a crust (like the ring in a
    bathtub). Each check, a molecule fully on dry land sticks with chance
    `stick` (it is frozen in place: no moving, no reacting). When the water
    covers it again it dissolves back with chance `dissolve` per check.
    """

    def __init__(self, stick=0.02, dissolve=0.1):
        self.stick, self.dissolve = stick, dissolve
        self.mask = np.zeros(0, dtype=bool)   # beads that are crust (not rock)
        self.stuck = self.freed = 0

    def update(self, world, pond):
        if len(self.mask) < world.n:   # new beads (rain, meteorites) start dissolved
            self.mask = np.concatenate([self.mask, np.zeros(world.n - len(self.mask), dtype=bool)])
        _, out, _, _ = pond.where(world._wrapped(), world.t)
        mol = world.mol
        dry = np.bincount(mol, weights=out > 0, minlength=world.n_mol) == np.bincount(mol, minlength=world.n_mol)
        rock = np.zeros(world.n_mol, dtype=bool)
        rock[mol[world.frozen & ~self.mask]] = True   # real rock never changes
        crust = np.zeros(world.n_mol, dtype=bool)
        crust[mol[self.mask]] = True
        roll = world.rng.uniform(size=world.n_mol)
        stick = dry & ~crust & ~rock & (roll < self.stick)
        free = ~dry & crust & (roll < self.dissolve)
        if stick.any() or free.any():
            change = stick[mol] | free[mol]
            self.mask[change] = stick[mol][change]
            world.frozen[change] = stick[mol][change]
            world.vel[change] = 0.0
            world._pairs = world._charged_pairs = None   # frozen pairs are skipped: rebuild lists
            self.stuck += int(stick.sum())
            self.freed += int(free.sum())

    def count(self, world):
        """Molecules in the crust now."""
        return len(np.unique(world.mol[self.mask[:world.n]]))


class Weather:
    """Four seasons, steady wind and storms over the ponds.

    One `year` has four seasons that blend smoothly into each other:
      spring  rainy, moderate wind
      summer  warm, ponds dry deepest, calm
      autumn  windy, storms most likely
      winter  cool, ponds hardly dry, calm
    Each season sets: wind (x `wind`), rain splash (x Pond.rain), dry (how
    far ponds shrink, 1 = down to r_dry), kT shift, storm (storm chance).

    Wind blows from one direction (heading(t), which wanders slowly). It
    pushes the water downwind through the middle of a pond; the water piles
    up at the far shore and flows back along both sides, so the pond splits
    into two turning cells. Smaller swirls ride on top and drift slowly. At
    our tiny scale the wind itself is not felt; what reaches the molecules
    is this stirring and shear (water moving faster on one side of a vesicle
    than the other).

    Storms come at random (mean gap `storm_every` steps, more in autumn),
    last about `storm_length` steps, rise and fade smoothly and make the wind
    `storm_wind` times stronger. A storm is a whirl passing over: its wind
    turns all the way round while it passes, and near its middle the water
    turns in one big circle (whirl_at(t), signed: + anticlockwise).
    Call update(t) once per chunk of steps (it rolls the storm dice).

    current: how hard the wind stirs the whole pond (1 = full, 0.3 = the
    water mostly stays put). surface: wind drag on floating things. Wind
    mostly drags the surface: floating lipids slide downwind and pile up on
    the far shore while the water goes back underneath (our flat world has
    no underneath, so this is a push on floating beads only, surface x wind
    speed, downwind). It holds them there as long as the wind keeps its
    direction, like scum and foam in the corner of a lake. 0 = off.

    carry_at(): in a turning current, water pressure pushes inward and keeps
    anything floating (same density as water) on its circle. Our beads are
    heavy and slip against the water, so without this push a storm flings
    everything to the shore like a salad spinner; real molecules just ride
    along. The push is mass x the acceleration of the water following its
    own flow, (u . grad) u.
    """

    SEASONS = ("spring", "summer", "autumn", "winter")
    TABLE = {
        "spring": dict(wind=0.6, rain=1.0, dry=0.6, kT=0.0, storm=0.5),
        "summer": dict(wind=0.3, rain=0.3, dry=1.0, kT=0.05, storm=0.3),
        "autumn": dict(wind=1.0, rain=0.8, dry=0.5, kT=0.0, storm=2.0),
        "winter": dict(wind=0.5, rain=0.5, dry=0.2, kT=-0.05, storm=0.2),
    }

    def __init__(self, year=80000, wind=1.5, storm_every=40000, storm_length=4000, storm_wind=4.0,
                 eddies=3, drift_period=30000, seed=0, current=1.0, surface=0.0):
        self.year, self.wind, self.eddies, self.drift_period = year, wind, eddies, drift_period
        self.current, self.surface = current, surface
        self.storm_every, self.storm_length, self.storm_wind = storm_every, storm_length, storm_wind
        self.rng = np.random.default_rng(seed)
        self.storms = []   # (start, length, strength, spin)
        self.last = None
        self.heading0 = self.rng.uniform(0, 2 * np.pi)   # prevailing wind direction
        self.wander = self.rng.uniform(0, 2 * np.pi, 2)

    def season(self, t):
        """Values for time t, blended between this season and the next; name = nearest season."""
        x = (t / self.year * 4) % 4
        i = int(x)
        f = 0.5 - 0.5 * np.cos(np.pi * (x - i))   # smooth 0 -> 1 across the season
        a, b = self.TABLE[self.SEASONS[i]], self.TABLE[self.SEASONS[(i + 1) % 4]]
        out = {k: (1 - f) * a[k] + f * b[k] for k in a}
        out["name"] = self.SEASONS[i] if x - i < 0.5 else self.SEASONS[(i + 1) % 4]
        return out

    def update(self, t):
        """Maybe start a storm (call once per chunk)."""
        if self.last is not None and self.storm_every and not self.storm_at(t):
            chance = (t - self.last) / self.storm_every * self.season(t)["storm"]
            if self.rng.uniform() < chance:
                length = self.storm_length * self.rng.uniform(0.5, 1.5)
                strength = self.storm_wind * self.rng.uniform(0.7, 1.3)
                spin = 1.0 if self.rng.uniform() < 0.5 else -1.0
                self.storms.append((t, length, strength, spin))
        self.last = t

    def storm_at(self, t):
        """0 when calm, rising to the storm's strength and back."""
        if not self.storms:
            return 0.0
        start, length, strength = self.storms[-1][:3]
        if not start <= t < start + length:
            return 0.0
        return strength * np.sin(np.pi * (t - start) / length) ** 2

    def _middle(self, t):
        """0..1, how close to the storm's middle (1 = the whirl is right overhead)."""
        if not self.storms:
            return 0.0
        start, length = self.storms[-1][:2]
        if not start <= t < start + length:
            return 0.0
        return np.sin(np.pi * (t - start) / length) ** 8

    def wind_at(self, t):
        """Straight wind speed now (season's wind + the storm's, minus what goes into the whirl)."""
        return self.wind * (self.season(t)["wind"] + self.storm_at(t) * (1 - self._middle(t)))

    def whirl_at(self, t):
        """Speed of the storm's circular current, signed (+ anticlockwise); 0 outside storms."""
        if not self.storms:
            return 0.0
        return self.storms[-1][3] * self.wind * self.storm_at(t) * self._middle(t)

    def heading(self, t):
        """Direction the wind blows towards (radians). Wanders slowly; a storm turns it full circle."""
        a = self.heading0 + 0.5 * np.sin(2 * np.pi * t / self.year + self.wander[0]) \
            + 0.3 * np.sin(2 * np.pi * t / (0.37 * self.year) + self.wander[1])
        if self.storms:
            start, length, _, spin = self.storms[-1]
            if t >= start:
                x = min((t - start) / length, 1.0)
                a += spin * 2 * np.pi * (x - np.sin(2 * np.pi * x) / (2 * np.pi))   # smooth full turn
        return a

    def kT_shift(self, t):
        return self.season(t)["kT"]

    def _steady(self, t, pond):
        """Wind current + storm whirl + drifting swirls (no rain splash): pos -> velocity."""
        W, spin, a = self.wind_at(t), self.whirl_at(t), self.heading(t)
        if self.current != 1.0:
            W, spin = W * self.current, spin * self.current
        side = np.array([-np.sin(a), np.cos(a)])   # left of the wind
        k = 2 * np.pi * self.eddies / pond.L
        px = 2 * np.pi * t / self.drift_period
        py = 2 * np.pi * t / (1.7 * self.drift_period)
        radii = pond.radii(t)

        def u(pos):
            # every pond's current, added up: where full ponds overlap a bead feels both (each
            # fades to rest at its own shore, so the sum has no jumps; a jump where the nearest
            # pond switched made carry_at blow up)
            d = pos[:, None, :] - pond.centers[None]
            d -= pond.L * np.round(d / pond.L)
            r = np.sqrt(np.einsum("ijk,ijk->ij", d, d))
            if pond.shapes is None:
                gap, g, gp = r - radii[None], np.ones_like(r), np.zeros_like(r)
                reach = radii
            else:
                gap, g, gp = pond.gaps(d, r, radii)
                reach = radii * (1 + pond._waves()[3])
            near = np.argmin(gap, axis=1)
            # only (bead, pond) pairs inside that pond's widest shoreline can be in its water
            i, j = np.nonzero(r < reach[None])
            d, r, g, gp, R = d[i, j], r[i, j], g[i, j], gp[i, j], radii[j]
            # both currents come from stream functions psi of s = r / (R g(angle)) (0 centre,
            # 1 shore), so water follows the shoreline and is never created or lost
            s = r / (R * g)
            x = np.clip(s, 0.0, 1.0)
            r_ = np.maximum(r, 1e-9)[..., None]
            radial = d / r_
            turn = np.stack([d[..., 1], -d[..., 0]], axis=-1) / r_
            grad = (radial + (gp / g)[..., None] * turn) / (R * g)[..., None]   # grad s
            inside = (s < 1)[..., None]
            # wind: psi = W n (1 - s^2)^2, n = distance to the left of the wind line through the
            # centre: speed W downwind through the middle, back upwind along both sides, and
            # slowing to rest at the shore (the bank drags)
            n = d @ side
            q = (1 - x * x)[..., None]
            pair = W * q * (q * side - (4 * n * x)[..., None] * grad) * inside
            # whirl: psi = spin R (2 s^2 - 4/3 s^3), circular, speed spin 4s(1-s) like a gyre
            pair -= (spin * R * (4 * x - 4 * x * x))[..., None] * grad * inside
            gpsi = np.zeros_like(pos)
            np.add.at(gpsi, i, pair)
            v = np.stack([gpsi[:, 1], -gpsi[:, 0]], axis=1)
            sx, cx = np.sin(k * pos[:, 0] + px), np.cos(k * pos[:, 0] + px)
            sy, cy = np.sin(k * pos[:, 1] + py), np.cos(k * pos[:, 1] + py)
            return v + 0.4 * (W + abs(spin)) * np.stack([sx * cy, -cx * sy], axis=1), near

        return u

    def carry_at(self, t, pond, every=10, h=0.01):
        """External force pos -> (force, 0): the pressure push that keeps floating beads with the water.

        Recomputed every `every` calls (beads move ~0.1 in that time)."""
        u = self._steady(t, pond)
        cache = {"n": 0, "F": None}

        def carry(pos):
            if cache["F"] is None or cache["n"] % every == 0 or len(cache["F"]) != len(pos):
                v, _ = u(pos)
                ex, ey = np.array([h, 0.0]), np.array([0.0, h])
                dudx = (u(pos + ex)[0] - u(pos - ex)[0]) / (2 * h)
                dudy = (u(pos + ey)[0] - u(pos - ey)[0]) / (2 * h)
                cache["F"] = v[:, :1] * dudx + v[:, 1:] * dudy   # mass 1
            cache["n"] += 1
            return cache["F"], 0.0

        return carry

    def drift_at(self, t, pond, floating, every=10):
        """External force pos -> (force, 0): wind drag pushing floating beads downwind (in water only).

        floating: bool mask of beads that float (lipids). Recomputed every `every` calls."""
        a = self.heading(t)
        push = self.surface * (self.wind_at(t) + abs(self.whirl_at(t))) * np.array([np.cos(a), np.sin(a)])
        cache = {"n": 0, "F": None}

        def drift(pos):
            if cache["F"] is None or cache["n"] % every == 0 or len(cache["F"]) != len(pos):
                wet = pond.where(pos, t)[1] <= 0
                cache["F"] = (floating & wet)[:, None] * push[None]
            cache["n"] += 1
            return cache["F"], 0.0

        return drift

    def flow_at(self, t, pond, splash):
        """Water velocity: wind current (+ storm whirl) in each pond + drifting swirls + rain splash."""
        rain_k = 2 * np.pi * pond.eddies / pond.L
        rx, ry = pond.rng.uniform(0, 2 * np.pi, 2)   # rain drops: a new splash pattern each call

        def swirl(pos, k, px, py):
            sx, cx = np.sin(k * pos[:, 0] + px), np.cos(k * pos[:, 0] + px)
            sy, cy = np.sin(k * pos[:, 1] + py), np.cos(k * pos[:, 1] + py)
            return np.stack([sx * cy, -cx * sy], axis=1)

        steady = self._steady(t, pond)

        def flow(pos):
            v, near = steady(pos)
            if splash.max() > 0:
                v += splash[near][:, None] * swirl(pos, rain_k, rx, ry)
            return v

        return flow
