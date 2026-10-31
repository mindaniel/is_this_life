"""Landscape level: a river basin of ponds in real units (metres, hours, moles).

The bead world (alife.chem) is about 100 nm across and lives for microseconds,
so it can't say where lipids go in a real landscape. This model can: each pond
is a few numbers (water volume, moles of lipid, formaldehyde, sugar) and an
hourly loop moves water and chemicals around:

- rain (wet and dry seasons, storms), runoff from the land around each pond
  (its catchment), evaporation (more by day, much more in hot ponds),
  seepage, hot springs that pulse on and off;
- a pond filled past its rim spills into its stream (downhill to the nearest
  lower pond) or, in a big storm or with no stream, floods over the land into
  the lower ponds nearby. The lowest pond (the lake) spills out of the basin;
- a pond that dries leaves its chemicals as a film on the mud; rain
  dissolves the film again (the wet-dry cycle).

Chemistry, per pond:
- lipid (fatty acid): arrives with hot-spring water and with dust/meteorite
  organics washed off the land; clay in rocky ponds slowly soaks it up.
  Vesicles can only exist above the critical vesicle concentration (CVC,
  ~10-40 mM for short fatty acids, lower with fatty alcohol mixed in).
- formaldehyde: made in the air by sunlight, falls with the rain.
- sugar (formose): 2+ formaldehyde -> sugars, catalysed by minerals, faster
  when hot. Sugars fall apart fast when hot (ribose half-life 73 min at
  100 C, 44 years at 0 C, Larralde 1995) and slowly when cool.

Rates nobody knows well (lipid supply, formose speed) are arguments; the
water numbers are typical for a warm volcanic island.
"""
import numpy as np

KINDS = ("plain", "hot", "rocky")


def ground(rng, size, relief=8.0, waves=6):
    """Smooth hills: height in metres at (x, y) points, a few long random waves."""
    angle = rng.uniform(0, 2 * np.pi, waves)
    wavelength = rng.uniform(0.5, 1.5, waves) * size
    k = 2 * np.pi / wavelength[:, None] * np.stack([np.cos(angle), np.sin(angle)], 1)
    amp = relief / 2 * rng.uniform(0.5, 1.0, waves) / np.sqrt(waves)
    phase = rng.uniform(0, 2 * np.pi, waves)
    return lambda p: (amp * np.cos(np.atleast_2d(p) @ k.T + phase)).sum(-1)


class Basin:
    """Ponds on hills; streams and floods run downhill; hourly water and chemistry."""

    def __init__(self, seed=3, ponds=9, size=300.0, lake=25.0, reach=40.0, kinds=None,
                 cvc=15.0, vent_lipid=0.05, dust=1e-9, clay=0.01, rain_hcho=0.01,
                 formose=0.05, hcho_life=30.0):
        rng = np.random.default_rng(seed)
        self.rng = rng
        self.size = size
        self.height = ground(rng, size)
        n = ponds
        # the lake sits in the lowest spot; the others are scattered, not touching
        g = np.linspace(lake + 10, size - lake - 10, 61)
        grid = np.stack(np.meshgrid(g, g), -1).reshape(-1, 2)
        low = np.quantile(self.height(grid), 0.5)        # ponds sit in hollows, not on hilltops
        lake_at = grid[np.argmin(self.height(grid))]   # lowest spot not at the map edge
        radius = np.r_[lake, rng.uniform(4.0, 12.0, n - 1)]
        centers = [lake_at]
        for i in range(1, n):
            for _ in range(5000):
                c = rng.uniform(radius[i] + 5, size - radius[i] - 5, 2)
                if self.height(c)[0] < low and all(np.hypot(*(c - q)) > radius[i] + radius[j] + 8 for j, q in enumerate(centers)):
                    centers.append(c)
                    break
            else:
                raise ValueError("ponds don't fit: fewer ponds or a bigger basin")
        self.centers = np.array(centers)
        self.radius = radius
        self.rim = self.height(self.centers)            # water spills at rim height
        self.depth = np.r_[3.0, rng.uniform(0.3, 1.5, n - 1)]
        self.area_max = np.pi * radius ** 2
        self.v_max = self.area_max * self.depth / 2      # bowl: area grows with depth
        self.kind = np.array(["plain"] + [KINDS[k] for k in rng.integers(0, 3, n - 1)]
                             if kinds is None else kinds.split(","))
        self.hot = self.kind == "hot"
        self.rocky = self.kind == "rocky"
        self.temp = np.where(self.hot, 70.0, 25.0)       # C; cool ponds follow the seasons
        self.streams = self._streams(reach)
        self.floods = self._flood_paths(2 * reach)
        self.catch = self._catchments()
        # chemistry settings
        self.cvc, self.vent_lipid, self.dust, self.clay = cvc, vent_lipid, dust, clay
        self.rain_hcho, self.formose, self.hcho_life = rain_hcho, formose, hcho_life
        # state: water (m^3) and moles in water / in the dried film on the mud
        self.v = self.v_max * 0.5
        self.lipid = self.v * 0.0
        self.hcho = self.v * 0.0
        self.sugar = self.v * 0.0                         # moles of sugar carbon
        self.film = np.zeros((3, n))                      # lipid, hcho, sugar left when dry
        self.land_lipid = np.zeros(n)                     # dust on each catchment, washed in by runoff
        self.spring_on = self.hot.copy()
        self.hour = 0
        self.spill = np.zeros((n, n))                     # m^3 moved i -> j this hour (for drawing)

    # ---- the map --------------------------------------------------------
    def _path_ok(self, i, j):
        """Can water from pond i's rim run to pond j? Not uphill, not through another pond."""
        a, b = self.centers[i], self.centers[j]
        d = np.hypot(*(b - a))
        e = (b - a) / d
        s = np.linspace(self.radius[i], d - self.radius[j], 12)
        if self.height(a + s[:, None] * e).max() > self.rim[i] + 0.05:
            return False
        for k in range(len(self.radius)):
            if k in (i, j):
                continue
            q = self.centers[k] - a
            u = np.clip(q @ e, 0, d)
            if np.hypot(*(q - u * e)) < self.radius[k]:
                return False
        return True

    def _gap(self, i, j):
        return np.hypot(*(self.centers[j] - self.centers[i])) - self.radius[i] - self.radius[j]

    def _streams(self, reach):
        """Each pond's stream: to the nearest lower pond it can reach (-1 = none; the lake drains out)."""
        n = len(self.radius)
        out = np.full(n, -1)
        for i in range(1, n):
            ok = [j for j in range(n) if self.rim[j] < self.rim[i] and self._gap(i, j) <= reach
                  and self._path_ok(i, j)]
            if ok:
                out[i] = min(ok, key=lambda j: self._gap(i, j))
        return out

    def _flood_paths(self, reach):
        """Where a flood over the land goes: lower ponds within reach, nearer ones get more."""
        n = len(self.radius)
        w = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if self.rim[j] < self.rim[i] and self._gap(i, j) <= reach and self._path_ok(i, j):
                    w[i, j] = 1 / max(self._gap(i, j), 1.0)
        rows = w.sum(1, keepdims=True)
        return np.divide(w, rows, out=np.zeros_like(w), where=rows > 0)

    def _catchments(self, cell=5.0):
        """Land area (m^2) whose rain runs into each pond: the nearest pond lower than that land."""
        g = np.arange(cell / 2, self.size, cell)
        p = np.stack(np.meshgrid(g, g), -1).reshape(-1, 2)
        h = self.height(p)
        d = np.hypot(*(p[:, None, :] - self.centers[None]).transpose(2, 0, 1))
        land = (d > self.radius[None]).all(1)
        d = np.where(self.rim[None] < h[:, None], d, np.inf)
        to = np.argmin(d, 1)
        ok = land & np.isfinite(d.min(1))
        return np.bincount(to[ok], minlength=len(self.radius)) * cell ** 2

    # ---- geometry of the water ------------------------------------------
    def area(self):
        """Water surface now: bowl-shaped ponds, so area ~ sqrt(volume)."""
        return self.area_max * np.sqrt(np.clip(self.v / self.v_max, 0, 1))

    def wet(self):
        return self.v > 1e-3 * self.v_max

    def conc(self, moles):
        """Concentration in mM (= mol/m^3); 0 where dry."""
        return np.where(self.wet(), moles / np.maximum(self.v, 1e-9), 0.0)

    # ---- weather --------------------------------------------------------
    def season(self):
        """0..1: 1 = middle of the wet season."""
        day = self.hour / 24 % 365
        return 0.5 + 0.5 * np.cos(2 * np.pi * (day - 60) / 365)

    def new_day(self):
        """Today's rain (m) and how many hours it falls."""
        wet = self.season()
        p = 0.1 + 0.5 * wet                               # chance of a rainy day
        if self.rng.random() > p:
            return 0.0, 0
        mm = self.rng.exponential(3 + 12 * wet)           # storms: a few days get 50-100 mm
        return mm / 1000, int(self.rng.integers(1, 6))

    # ---- one hour -------------------------------------------------------
    def step(self, rain):
        """Advance one hour. rain = metres of rain this hour."""
        h = self.hour % 24
        wet0 = self.wet()
        area = self.area()
        # water in
        runoff = 0.4 * max(rain - 0.002, 0.0)             # light rain soaks into the ground
        inflow = rain * area + runoff * self.catch
        if self.hour % 24 == 0:                            # hot springs pulse: on/off for days to weeks
            flip = self.rng.random(len(self.v)) < 1 / 10
            self.spring_on = np.where(self.hot & flip, ~self.spring_on, self.spring_on) & self.hot
        spring = np.where(self.spring_on, 0.15, 0.0)        # m^3/hour
        # water out: evaporation peaks at noon, x3 in hot ponds; seepage through the bottom
        sun = max(0.0, np.sin(np.pi * (h - 6) / 12))
        evap = (0.004 + 0.002 * (1 - self.season())) * np.pi / 24 * sun   # m/hour, 4-6 mm a day
        evap = evap * np.where(self.hot, 3.0, 1.0)
        seep = np.where(self.rocky, 0.0005, 0.002) / 24
        out = (evap + seep) * area
        # chemicals in: rain formaldehyde, spring lipid, dust on the land (washed in by runoff)
        self.hcho += rain * area * self.rain_hcho   # mM = mol/m^3
        self.lipid += spring * self.vent_lipid
        self.land_lipid += self.dust / 24 * self.catch
        washed = self.land_lipid * min(1.0, runoff / 0.02)
        self.land_lipid -= washed
        self.lipid += washed + self.dust / 24 * area
        self.v = np.maximum(self.v + inflow + spring - out, 0.0)
        # drying / rewetting: chemicals go to the film and back
        wet = self.wet()
        dried, rewet = wet0 & ~wet, ~wet0 & wet
        for k, name in enumerate(("lipid", "hcho", "sugar")):
            m = getattr(self, name)
            self.film[k] += np.where(dried, m, 0.0)
            m = np.where(dried, 0.0, m) + np.where(rewet, self.film[k], 0.0)
            self.film[k] = np.where(rewet, 0.0, self.film[k])
            setattr(self, name, m)
        self.film[1] *= np.where(wet, 1.0, 0.9)            # dry formaldehyde mostly evaporates
        self.rewet = rewet
        self._chemistry(wet)
        self._spill()
        self.hour += 1

    def _chemistry(self, wet):
        cool = 25 + 5 * (self.season() - 0.5)
        self.temp = np.where(self.hot, 70.0, cool)
        T = self.temp + 273.15
        # sugar half-life (hours) from Larralde 1995 (ribose): 73 min at 100 C, 44 y at 0 C
        sugar_life = 73 / 60 * np.exp(12900 * (1 / T - 1 / 373.15))
        sugar_life = np.where(wet, sugar_life, sugar_life * 10)   # dry sugar keeps longer
        c = self.conc(self.hcho)
        speed = np.where(self.rocky, 1.0, np.where(self.hot, 0.5, 0.1)) * np.exp(6000 * (1 / 298.15 - 1 / T))
        made = np.minimum(self.formose / 24 * speed * c ** 2 * self.v, self.hcho * 0.5)
        self.hcho -= made + self.hcho * (1 - 0.5 ** (1 / (self.hcho_life * 24)))
        self.sugar += made
        self.sugar *= 0.5 ** (1 / sugar_life)
        self.film[2] *= 0.5 ** (1 / (sugar_life * 10))
        self.lipid -= np.where(self.rocky & wet, self.clay / 24 * self.lipid, 0.0)

    def _spill(self):
        """Overfull ponds spill, highest first, so a cascade runs down in the same hour."""
        self.spill[:] = 0.0
        for i in np.argsort(-self.rim):
            extra = self.v[i] - self.v_max[i]
            if extra <= 0:
                continue
            share = extra / self.v[i]
            moved = np.array([self.lipid[i], self.hcho[i], self.sugar[i]]) * share
            self.v[i] = self.v_max[i]
            self.lipid[i] -= moved[0]
            self.hcho[i] -= moved[1]
            self.sugar[i] -= moved[2]
            # a stream carries up to 20 m^3/hour; the rest (or all, with no stream) floods over the land
            to = {}
            j = self.streams[i]
            if j >= 0:
                to[j] = min(extra, 20.0)
            rest = extra - sum(to.values())
            for k in np.flatnonzero(self.floods[i]):
                to[k] = to.get(k, 0.0) + rest * self.floods[i, k]
            for k, vol in to.items():
                f = vol / extra
                self.v[k] += vol
                self.lipid[k] += moved[0] * f
                self.hcho[k] += moved[1] * f
                self.sugar[k] += moved[2] * f
                self.spill[i, k] += vol

    # ---- a whole run -----------------------------------------------------
    def run(self, years=10, every=6, callback=None):
        """Run hourly; log every `every` hours. Returns a dict of arrays (time, pond)."""
        hours = int(years * 365 * 24)
        log = {k: [] for k in ("hour", "volume", "lipid_mM", "hcho_mM", "sugar_mM", "temp", "flow")}
        n = len(self.v)
        self.events = np.zeros(n, int)                     # rewetting a film rich enough to make vesicles
        self.rewets = np.zeros(n, int)
        self.ves_hours = np.zeros(n)
        self.sweet_hours = np.zeros(n)                     # vesicles + sugar (>1 mM C) together
        self.wet_hours = np.zeros(n)
        flow = np.zeros((n, n))
        rain_left, rain_mm = 0, 0.0
        for t in range(hours):
            if t % 24 == 0:
                rain_mm, rain_left = self.new_day()
                start = self.rng.integers(0, 24 - rain_left) if rain_left else 0
            raining = rain_left and start <= t % 24 < start + rain_left
            self.step(rain_mm / rain_left if raining else 0.0)
            flow += self.spill
            wet = self.wet()
            lip = self.conc(self.lipid)
            ves = wet & (lip > self.cvc)
            self.wet_hours += wet
            self.ves_hours += ves
            self.sweet_hours += ves & (self.conc(self.sugar) > 1.0)
            self.rewets += self.rewet
            self.events += self.rewet & ves
            if t % every == 0:
                log["hour"].append(self.hour)
                log["volume"].append(self.v / self.v_max)
                log["lipid_mM"].append(lip)
                log["hcho_mM"].append(self.conc(self.hcho))
                log["sugar_mM"].append(self.conc(self.sugar))
                log["temp"].append(self.temp.copy())
                log["flow"].append(flow / every)
                flow = np.zeros((n, n))
                if callback:
                    callback(self, log)
        return {k: np.array(v) for k, v in log.items()}
