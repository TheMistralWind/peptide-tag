"""Peptide chemistry for the name-to-peptide tool.

Everything here is meant to be literally true. The residue table below was
verified against RDKit's own peptide builder (Chem.MolFromSequence) for all 20
standard amino acids in four sequence contexts each, plus 40 random sequences,
and every free amino acid formula was checked against its reference formula.
See tests/test_chemistry.py, which re-runs those checks.

RDKit refuses B, J, O, U, X and Z, so we build peptide SMILES from residue
templates rather than calling MolFromSequence. That keeps pyrrolysine (O) and
selenocysteine (U) as real residues instead of silently swapping them out.
"""

from __future__ import annotations

import re
import unicodedata
from dataclasses import dataclass, field

from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Classes drive the colour coding in the artwork. "special" is the usual bucket
# for the three residues that break the normal rules: cysteine bridges, glycine
# has no side chain, proline kinks the backbone.
HYDROPHOBIC = "hydrophobic"
AROMATIC = "aromatic"
POLAR = "polar"
ACIDIC = "acidic"
BASIC = "basic"
SPECIAL = "special"
RARE = "rare"

CLASS_LABEL = {
    HYDROPHOBIC: "Hydrophobic",
    AROMATIC: "Aromatic",
    POLAR: "Polar",
    ACIDIC: "Acidic",
    BASIC: "Basic",
    SPECIAL: "Special",
    RARE: "Rare",
}


@dataclass(frozen=True)
class Residue:
    letter: str
    name: str
    code: str
    side_chain: str | None  # SMILES of the side chain hung off the alpha carbon
    klass: str
    hydropathy: float  # Kyte-Doolittle
    formula: str = ""  # the FREE amino acid, before any peptide bond
    note: str = ""


# side_chain is written so that N[C@@H](<side_chain>)C(=O) is the L residue.
# Glycine and proline are special-cased in _residue_smiles.
RESIDUES: dict[str, Residue] = {
    "A": Residue("A", "Alanine", "Ala", "C", HYDROPHOBIC, 1.8, "C3H7NO2"),
    "C": Residue("C", "Cysteine", "Cys", "CS", SPECIAL, 2.5, "C3H7NO2S",
                 note="Forms disulfide bridges, the staples that hold proteins shut."),
    "D": Residue("D", "Aspartic acid", "Asp", "CC(=O)O", ACIDIC, -3.5, "C4H7NO4"),
    "E": Residue("E", "Glutamic acid", "Glu", "CCC(=O)O", ACIDIC, -3.5, "C5H9NO4"),
    "F": Residue("F", "Phenylalanine", "Phe", "Cc1ccccc1", AROMATIC, 2.8, "C9H11NO2"),
    "G": Residue("G", "Glycine", "Gly", None, SPECIAL, -0.4, "C2H5NO2",
                 note="The only residue with no side chain, so it bends where others cannot."),
    "H": Residue("H", "Histidine", "His", "Cc1c[nH]cn1", BASIC, -3.2, "C6H9N3O2"),
    "I": Residue("I", "Isoleucine", "Ile", "[C@@H](C)CC", HYDROPHOBIC, 4.5, "C6H13NO2"),
    "K": Residue("K", "Lysine", "Lys", "CCCCN", BASIC, -3.9, "C6H14N2O2"),
    "L": Residue("L", "Leucine", "Leu", "CC(C)C", HYDROPHOBIC, 3.8, "C6H13NO2"),
    "M": Residue("M", "Methionine", "Met", "CCSC", HYDROPHOBIC, 1.9, "C5H11NO2S",
                 note="Every protein you have ever made started with this one."),
    "N": Residue("N", "Asparagine", "Asn", "CC(N)=O", POLAR, -3.5, "C4H8N2O3"),
    "O": Residue("O", "Pyrrolysine", "Pyl", "CCCCNC(=O)[C@@H]1CC=N[C@H]1C",
                 RARE, -3.9, "C12H21N3O3",
                 note="The 22nd amino acid. Methane-making microbes build it onto "
                      "a lysine and smuggle it in by reading a stop codon as a letter."),
    "P": Residue("P", "Proline", "Pro", None, SPECIAL, -1.6, "C5H9NO2",
                 note="Its side chain loops back onto the backbone, forcing a kink."),
    "Q": Residue("Q", "Glutamine", "Gln", "CCC(N)=O", POLAR, -3.5, "C5H10N2O3"),
    "R": Residue("R", "Arginine", "Arg", "CCCNC(N)=N", BASIC, -4.5, "C6H14N4O2"),
    "S": Residue("S", "Serine", "Ser", "CO", POLAR, -0.8, "C3H7NO3"),
    "T": Residue("T", "Threonine", "Thr", "[C@@H](C)O", POLAR, -0.7, "C4H9NO3"),
    "U": Residue("U", "Selenocysteine", "Sec", "C[SeH]", RARE, 2.5, "C3H7NO2Se",
                 note="The 21st amino acid. Cysteine with selenium swapped in for "
                      "sulfur. You have 25 proteins that use it."),
    "V": Residue("V", "Valine", "Val", "C(C)C", HYDROPHOBIC, 4.2, "C5H11NO2"),
    "W": Residue("W", "Tryptophan", "Trp", "Cc1c[nH]c2ccccc12", AROMATIC, -0.9,
                 "C11H12N2O2",
                 note="The biggest residue, and the one that makes proteins glow under UV."),
    "Y": Residue("Y", "Tyrosine", "Tyr", "Cc1ccc(O)cc1", AROMATIC, -1.3, "C9H11NO3"),
}

# Four letters of the alphabet are not amino acids at all. All four are real
# IUPAC ambiguity codes, so resolving each to a concrete residue is standard
# practice rather than a fudge, and saying which is more interesting than
# hiding it. O and U are deliberately NOT in here: they are genuinely the 22nd
# and 21st amino acids, and keeping them means names like ROOSA come out
# completely unchanged, which is the whole point of the toy.
SUBSTITUTIONS: dict[str, tuple[str, str]] = {
    "B": ("N", "B is Asx, the code used when a reading cannot distinguish "
               "asparagine from aspartic acid. It resolves here to asparagine."),
    "J": ("L", "J is Xle, the code for leucine or isoleucine. The two have "
               "identical mass, so instruments often cannot tell them apart. "
               "It resolves here to leucine."),
    "X": ("G", "X is the code for any amino acid at all. It resolves here to "
               "glycine, the smallest and simplest one."),
    "Z": ("Q", "Z is Glx, the code used when glutamine cannot be distinguished "
               "from glutamic acid. It resolves here to glutamine."),
}

# Pyrrolysine and selenocysteine are real, but no human protein contains
# pyrrolysine and only 25 contain selenocysteine, so a search for your name in
# the human proteome would always fail for any name with an O in it. For that
# search only, each is read as the residue it is built from. The peptide itself
# keeps the real residue.
PROTEOME_PROJECTION = {"O": "K", "U": "C"}


def to_standard(sequence: str) -> str:
    """Sequence rewritten using only the standard 20, for proteome search."""
    return "".join(PROTEOME_PROJECTION.get(a, a) for a in sequence)


# Side-chain charge at blood pH, used for the +/- marks in the artwork. This is
# deliberately not the same as the class: pyrrolysine looks like a lysine but
# its amine is capped by an amide, so it carries no charge, and selenocysteine
# looks like a cysteine but its selenol is acidic enough to be deprotonated.
CHARGE: dict[str, int] = {
    "K": 1, "R": 1, "H": 1,
    "D": -1, "E": -1, "U": -1,
}

# Histidine is only partly protonated at pH 7.4, so the mark overstates it
# slightly. The residue list says so rather than the artwork pretending.
PARTIAL_CHARGE = {"H"}

MAX_LENGTH = 60


@dataclass
class Swap:
    position: int  # 1-based position in the final sequence
    original: str
    replacement: str
    reason: str


@dataclass
class Peptide:
    raw: str
    sequence: str
    swaps: list[Swap] = field(default_factory=list)
    truncated: bool = False

    @property
    def length(self) -> int:
        return len(self.sequence)


# Nordic and German letters that decomposition alone does not resolve.
_TRANSLITERATE = str.maketrans({
    "Ø": "O", "ø": "o", "Æ": "AE", "æ": "ae", "Ð": "D", "ð": "d",
    "Þ": "TH", "þ": "th", "ß": "ss", "Ł": "L", "ł": "l", "Đ": "D", "đ": "d",
})


def _latinise(text: str) -> str:
    """Strip accents so Zoë, Jää and Håkan keep all of their letters."""
    text = (text or "").translate(_TRANSLITERATE)
    decomposed = unicodedata.normalize("NFKD", text)
    return "".join(c for c in decomposed if not unicodedata.combining(c))


def text_to_sequence(text: str) -> Peptide:
    """Turn arbitrary text into a peptide sequence, recording what changed."""
    letters = re.sub(r"[^A-Za-z]", "", _latinise(text)).upper()
    truncated = len(letters) > MAX_LENGTH
    letters = letters[:MAX_LENGTH]

    out: list[str] = []
    swaps: list[Swap] = []
    for i, ch in enumerate(letters, start=1):
        if ch in RESIDUES:
            out.append(ch)
        else:
            replacement, reason = SUBSTITUTIONS[ch]
            out.append(replacement)
            swaps.append(Swap(i, ch, replacement, reason))
    return Peptide(raw=text, sequence="".join(out), swaps=swaps, truncated=truncated)


def _residue_smiles(letter: str) -> str:
    if letter == "G":
        return "NCC(=O)"
    if letter == "P":
        # Proline's nitrogen is part of its own ring, so it cannot use the
        # generic backbone template.
        return "N1CCC[C@H]1C(=O)"
    return f"N[C@@H]({RESIDUES[letter].side_chain})C(=O)"


def build_smiles(sequence: str) -> str:
    """Correct SMILES for the free peptide, N-terminus to C-terminus."""
    if not sequence:
        return ""
    return "".join(_residue_smiles(a) for a in sequence) + "O"


_ATOM = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for element, number in _ATOM.findall(formula):
        counts[element] = counts.get(element, 0) + int(number or 1)
    return counts


def _hill(counts: dict[str, int]) -> str:
    """Carbon, then hydrogen, then everything else alphabetically. A count of
    one is written bare: the old code emitted S1, which no chemist writes."""
    order = (["C"] if counts.get("C") else []) + (["H"] if counts.get("H") else [])
    order += sorted(e for e in counts if e not in ("C", "H") and counts[e])
    return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}" for e in order)


def molecular_formula(sequence: str) -> str:
    """Hill-notation formula for the free peptide.

    Sum the free amino acids, then remove one water per peptide bond. Water is
    H2O, so each bond costs TWO hydrogens and one oxygen. The previous version
    removed one of each, which produced impossible formulas with a half-integer
    degree of unsaturation for every even-length name.

    tests/test_chemistry.py checks every result here against RDKit's own
    formula for the same molecule, so this arithmetic is machine-verified
    without RDKit being needed to serve a request.
    """
    if not sequence:
        return ""
    counts: dict[str, int] = {}
    for a in sequence:
        for element, n in parse_formula(RESIDUES[a].formula).items():
            counts[element] = counts.get(element, 0) + n
    bonds = len(sequence) - 1
    counts["H"] -= 2 * bonds
    counts["O"] -= bonds
    return _hill(counts)


def format_formula_html(formula: str) -> str:
    """C32H56N12O9 -> C<sub>32</sub>H<sub>56</sub>... Counts of 1 stay bare."""
    return re.sub(r"([A-Z][a-z]?)(\d*)",
                  lambda m: m.group(1) + (f"<sub>{m.group(2)}</sub>" if m.group(2) else ""),
                  formula)


_SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def format_formula_unicode(formula: str) -> str:
    """C32H56N12O9 -> C₃₂H₅₆N₁₂O₉, for images where markup is unavailable."""
    return formula.translate(_SUBSCRIPT)


def gravy(sequence: str) -> float:
    """Kyte-Doolittle average hydropathy.

    Biopython's own gravy() raises KeyError on pyrrolysine and selenocysteine,
    so we do it here. Those two borrow the value of the residue they are built
    from (lysine and cysteine respectively), which is an approximation and is
    labelled as one in the UI.
    """
    if not sequence:
        return 0.0
    return sum(RESIDUES[a].hydropathy for a in sequence) / len(sequence)


# pKa values, Bjellqvist-style, matching Biopython's defaults so our numbers
# agree with the standard tools. Selenocysteine is the one addition: its selenol
# is far more acidic than a thiol (pKa about 5.2 against cysteine's 8.5), and
# Biopython drops it silently, which would overstate pI for every name with a U.
_PKA_POSITIVE = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
_PKA_NEGATIVE = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0, "U": 5.2}
# The terminal groups shift depending on which residue carries them.
_PKA_NTERM = {"A": 7.59, "M": 7.00, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44, "E": 7.70}
_PKA_CTERM = {"D": 4.55, "E": 4.75}


def net_charge(sequence: str, ph: float) -> float:
    """Net charge of the peptide at a given pH."""
    if not sequence:
        return 0.0
    counts: dict[str, int] = {}
    for a in sequence:
        counts[a] = counts.get(a, 0) + 1

    charge = 0.0
    for group, pka in _PKA_POSITIVE.items():
        if group == "Nterm":
            n, pka = 1, _PKA_NTERM.get(sequence[0], pka)
        else:
            n = counts.get(group, 0)
        if n:
            charge += n * (1.0 / (1.0 + 10 ** (ph - pka)))
    for group, pka in _PKA_NEGATIVE.items():
        if group == "Cterm":
            n, pka = 1, _PKA_CTERM.get(sequence[-1], pka)
        else:
            n = counts.get(group, 0)
        if n:
            charge -= n * (1.0 / (1.0 + 10 ** (pka - ph)))
    return charge


def isoelectric_point(sequence: str) -> float:
    """pH at which the peptide carries no net charge, by bisection."""
    if not sequence:
        return 0.0
    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2
        if net_charge(sequence, mid) > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2, 2)


def composition(sequence: str) -> list[dict]:
    """Residues present, most frequent first, with their class and count."""
    counts: dict[str, int] = {}
    for a in sequence:
        counts[a] = counts.get(a, 0) + 1
    rows = []
    for letter, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        r = RESIDUES[letter]
        rows.append({
            "letter": letter, "name": r.name, "code": r.code,
            "klass": r.klass, "class_label": CLASS_LABEL[r.klass],
            "count": count, "note": r.note,
        })
    return rows


def analyse(sequence: str) -> dict:
    """Physicochemical properties. Safe on empty and exotic sequences."""
    if not sequence:
        return {}

    ana = ProteinAnalysis(sequence)

    def safe(fn, default=None):
        try:
            return fn()
        except Exception:
            return default

    charged_pos = sum(1 for a in sequence if RESIDUES[a].klass == BASIC)
    charged_neg = sum(1 for a in sequence if RESIDUES[a].klass == ACIDIC)

    props = {
        "length": len(sequence),
        "molecular_weight": safe(ana.molecular_weight, 0.0),
        "isoelectric_point": isoelectric_point(sequence),
        "charge_at_ph7": net_charge(sequence, 7.4),
        "aromaticity": safe(ana.aromaticity, 0.0),
        "gravy": gravy(sequence),
        "formula": molecular_formula(sequence),
        "positive": charged_pos,
        "negative": charged_neg,
    }
    props["formula_html"] = format_formula_html(props["formula"])
    props["water_soluble"] = props["gravy"] < 0
    return props
