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
    """A round pond that dries up and refills: wet-dry cycles.

    Water covers a disk of radius R(t) in the middle of the box. Outside is
    dry land. Sun evaporates the pond (R shrinks to r_dry), rain refills it
    (R grows back to r_wet), once per `period` steps. The run starts wet.

    Beads beyond the shore feel a spring pushing them back into the water:
    the retreating edge sweeps lipids together, so they get crowded while
    the pond is small and spread out again after rain. (Real drying strands
    lipids as films on the mud; this is the simplest honest version.)
    period = 0 keeps the pond always wet, for comparison.

    rain > 0: refilling stirs the water. Rain drops and inflowing water make
    swirls whose speed follows how fast the pond is refilling (top speed
    `rain` at the fastest refill, nothing while it dries or stays full).
    """

    def __init__(self, L, r_wet, r_dry, period, k=5.0, rain=0.0, eddies=3, seed=0):
        self.L, self.r_wet, self.r_dry, self.period, self.k = float(L), r_wet, r_dry, period, k
        self.rain, self.eddies = rain, eddies
        self.rng = np.random.default_rng(seed)

    def radius(self, t):
        if not self.period:
            return self.r_wet
        wet = 0.5 + 0.5 * np.cos(2 * np.pi * t / self.period)   # 1 = full, 0 = driest
        return self.r_dry + (self.r_wet - self.r_dry) * wet

    def refill_speed(self, t):
        """0..1: 0 while drying or full, 1 at the fastest refill."""
        if not self.period:
            return 0.0
        return max(0.0, -np.sin(2 * np.pi * t / self.period))   # R = cos: refills in the 2nd half

    def flow_at(self, t):
        """Rain-stirred water: random swirls (new pattern each call), or None."""
        speed = self.rain * self.refill_speed(t)
        if speed <= 0:
            return None
        k = 2 * np.pi * self.eddies / self.L
        px, py = self.rng.uniform(0, 2 * np.pi, 2)

        def flow(pos):
            sx, cx = np.sin(k * pos[:, 0] + px), np.cos(k * pos[:, 0] + px)
            sy, cy = np.sin(k * pos[:, 1] + py), np.cos(k * pos[:, 1] + py)
            return speed * np.stack([sx * cy, -cx * sy], axis=1)

        return flow

    def shore_at(self, t, frozen=None):
        """External force for World: pos -> (force, energy). Rock is left alone."""
        R, k, L = self.radius(t), self.k, self.L

        def shore(pos):
            d = pos - L / 2
            d -= L * np.round(d / L)
            r = np.sqrt(np.einsum("ij,ij->i", d, d))
            out = np.maximum(r - R, 0.0)   # how far onto dry land
            if frozen is not None:
                out[frozen] = 0.0
            F = (-k * out / np.maximum(r, 1e-9))[:, None] * d
            return F, float(np.sum(0.5 * k * out ** 2))

        return shore

