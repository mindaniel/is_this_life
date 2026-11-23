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
    ("stail",    0,   1.00, 1.0, "#F2CD7E", "short-chain tail (fatty-acid-like); Rules.set_oily tunes it"),
    ("hcn",      0,   0.80, 0.0, "#2A6FB0", "hydrogen cyanide, HCN. Alone = HCN; chains = HCN polymer (dark, absorbs light)"),
    ("nh3",      0,   0.75, 0.0, "#A3C9E8", "ammonia, NH3 (a base: speeds HCN polymer; N source for amino acids)"),
    ("cat",      0,   1.00, 0.0, "#2E8B57", "catalyst: a big folded molecule (short peptide / ribozyme), can't cross a membrane"),
    ("fa_acid",  0,   0.95, 0.0, "#C98BB0", "fatty acid head, COOH (H+ on: neutral). Protonation swaps it with fa_head"),
    ("fa_oh",    0,   0.95, 0.0, "#E8C3A0", "fatty alcohol head, OH (neutral, never charged): decanol-like co-surfactant"),
    ("fuel",     0,   1.00, 0.0, "#F2A93B", "fuel: an activated molecule (anhydride / thioester-like); burning it frees an acid H+"),
    ("waste",    0,   1.00, 0.0, "#8C8C8C", "burnt fuel"),
    ("tar",      0,   0.80, 0.5, "#4A3320", "HCN polymer grown long (Nitrogen tar_len): insoluble brown-black tar; oiliness a guess"),
    ("nx",       0,   1.00, 0.0, "#C0392B", "base X of an RNA-like monomer (pairs with Y); hangs off one nbb backbone bead"),
    ("ny",       0,   1.00, 0.0, "#2471A3", "base Y of an RNA-like monomer (pairs with X)"),
    ("nbb",      0,   1.00, 0.0, "#B9A27A", "RNA-like backbone (sugar-phosphate): monomers join into strands here; charge left out"),
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
    ("head", "stail"): 0.95,
    ("fa_head", "stail"): 0.95,
    ("fa_acid", "tail"): 0.95,
    ("fa_acid", "stail"): 0.95,
    ("fa_oh", "tail"): 0.95,
    ("fa_oh", "stail"): 0.95,
}

# Attractions that are not oiliness: a COOH head and a COO- head next to each other share
# the H+ (hydrogen bond), so half-charged fatty-acid walls hold together best. Strength
# in the same units as oily (tail-tail = 1); the value is a guess, tune with
# Rules.set_pair_eps. COOH-COOH pairs (dimers) are left out: in water they barely form.
PAIR_EPS = {
    ("fa_acid", "fa_head"): 0.5,
    ("fa_oh", "fa_head"): 0.5,   # alcohol OH gives a hydrogen bond to a COO- (why decanol steadies decanoate walls)
}
