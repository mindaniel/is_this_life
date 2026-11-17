"""The sun as seen from a pond: a smooth day, and a planet at some distance from its star.

Day and night come from the planet spinning: brightness rises and falls as max(0, sin) over a
`period` (day = first half). Scaled by pi/2 so a day gives as much light as the old on/off switch
(mean of max(0, sin) is 1/pi, of on/off 1/2). The water warms a little later than the light
(peak in the afternoon): kT(t) = kT_mean (1 + swing sin(2 pi t / period - lag)).

Distance d (Earth = 1): brightness 1 / d^2 (light spreads over a bigger sphere), mean temperature
kT_earth d^-0.5 (a planet's equilibrium temperature goes as d^-1/2; Earth's greenhouse kept, so
d = 1 is our usual kT 1.1 ~ 25 C and kT is proportional to absolute temperature).
Shortcuts (flagged): no atmosphere change with distance (real Venus / Mars fates were set by their
air), no freezing or boiling (water is implicit), one temperature for the whole box, the day is
~10^14 times shorter than a real one. In the bead world light comes from everywhere at once
(the box is ~100 nm, far thinner than anything that casts a shadow), so there is no direction.
"""
import numpy as np


class Sun:
    def __init__(self, d=1.0, period=20000, kT_earth=1.1, swing=0.03, lag=np.pi / 4):
        self.d, self.period, self.kT_earth, self.swing, self.lag = d, period, kT_earth, swing, lag
        self.kT_mean = kT_earth * d ** -0.5

    def brightness(self, t):
        """Light now, relative to Earth's old on/off day (0 at night)."""
        return 0.5 * np.pi * max(0.0, np.sin(2 * np.pi * t / self.period)) / self.d ** 2

    def kT(self, t):
        return self.kT_mean * (1 + self.swing * np.sin(2 * np.pi * t / self.period - self.lag))
