"""The 'periodic table' of the sim: every bead type and its properties.

A bead is a small group of atoms bundled together (coarse-grained). Each row:

    name     short id used in code
    charge   electric charge in units of e (+ attracts -, like repels)
    size     effective diameter in sigma. For ions this is the size WITH the
             shell of water molecules they hold around them (hydrated size)
    oily     0..1, how much it hides from water. Oily beads attract each
             other (hydrophobic effect); strength = sqrt(oily_a * oily_b)
    color    for drawing
    note     what it stands for

Adding an element = adding a row. No new code.
"""

ELEMENTS = [
    # name      charge size  oily  color      note
    ("head",     0,   0.95, 0.0, "#D85A30", "neutral lipid head (phosphocholine-like)"),
    ("tail",     0,   1.00, 1.0, "#EF9F27", "oily carbon chain, ~3-4 CH2 groups"),
    ("fa_head", -1,   0.95, 0.0, "#993556", "fatty acid head, COO- (charged, early-Earth lipid)"),
    ("Na",      +1,   0.90, 0.0, "#378ADD", "sodium ion, seawater's main + ion"),
    ("K",       +1,   0.80, 0.0, "#7F77DD", "potassium ion, smaller water shell than Na"),
    ("Cl",      -1,   0.85, 0.0, "#1D9E75", "chloride ion, seawater's main - ion"),
    ("Mg",      +2,   1.20, 0.0, "#97C459", "magnesium ion, holds a tight thick water shell"),
    ("Ca",      +2,   0.90, 0.0, "#5F5E5A", "calcium ion, sheds water easily, grabs COO- directly"),
    ("pigment",  0,   1.00, 0.8, "#3B6D11", "conjugated ring, oily; absorbs light (later)"),
    ("mineral", -1,   1.00, 0.0, "#888780", "clay/mineral surface site, frozen in place"),
    ("rock",     0,   1.00, 0.0, "#444441", "plain rock, frozen; builds walls and pools"),
    ("sugar",    0,   1.00, 0.0, "#D4537E", "dissolved molecule too big to cross a membrane (osmotic)"),
    ("carb",     0,   0.80, 0.0, "#8E5BB5", "one sugar carbon, CH2O. Alone = formaldehyde; chains = sugars"),
]

NAMES = [row[0] for row in ELEMENTS]
INDEX = {name: k for k, name in enumerate(NAMES)}
CHARGE = [row[1] for row in ELEMENTS]
SIZE = [row[2] for row in ELEMENTS]
OILY = [row[3] for row in ELEMENTS]
COLOR = [row[4] for row in ELEMENTS]

# Exceptions to the size mixing rule (mean of the two sizes).
# Cooke-Kremer-Deserno lipids use 0.95 between a head and a tail.
PAIR_SIZE = {
    ("head", "tail"): 0.95,
    ("fa_head", "tail"): 0.95,
}
