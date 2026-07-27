"""Checks that the chemistry we display is actually true.

Run with: .venv/bin/python -m pytest tests/ -q
"""

import random
import string
import sys
from pathlib import Path

import pytest
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import chemistry as C  # noqa: E402

RDLogger.DisableLog("rdApp.*")

STANDARD_20 = "ACDEFGHIKLMNPQRSTVWY"

# Reference formulas for the free amino acids.
FREE_AA_FORMULA = {
    "A": "C3H7NO2", "R": "C6H14N4O2", "N": "C4H8N2O3", "D": "C4H7NO4",
    "C": "C3H7NO2S", "E": "C5H9NO4", "Q": "C5H10N2O3", "G": "C2H5NO2",
    "H": "C6H9N3O2", "I": "C6H13NO2", "L": "C6H13NO2", "K": "C6H14N2O2",
    "M": "C5H11NO2S", "F": "C9H11NO2", "P": "C5H9NO2", "S": "C3H7NO3",
    "T": "C4H9NO3", "W": "C11H12N2O2", "Y": "C9H11NO3", "V": "C5H11NO2",
    "U": "C3H7NO2Se", "O": "C12H21N3O3",
}


@pytest.mark.parametrize("letter", sorted(FREE_AA_FORMULA))
def test_free_amino_acid_formula(letter):
    """Each residue template, capped as a free amino acid, has the right formula."""
    assert C.molecular_formula(letter) == FREE_AA_FORMULA[letter]


@pytest.mark.parametrize("letter", list(STANDARD_20))
@pytest.mark.parametrize("context", ["{a}", "G{a}G", "{a}{a}{a}", "A{a}W"])
def test_matches_rdkit_peptide_builder(letter, context):
    """Our SMILES is identical to RDKit's own builder for the standard 20."""
    seq = context.format(a=letter)
    ours = Chem.MolFromSmiles(C.build_smiles(seq))
    theirs = Chem.MolFromSequence(seq)
    assert ours is not None and theirs is not None
    assert Chem.MolToSmiles(ours) == Chem.MolToSmiles(theirs)


def test_random_sequences_match_rdkit():
    rng = random.Random(20260727)
    for _ in range(60):
        seq = "".join(rng.choice(STANDARD_20) for _ in range(rng.randint(2, 40)))
        ours = Chem.MolFromSmiles(C.build_smiles(seq))
        theirs = Chem.MolFromSequence(seq)
        assert Chem.MolToSmiles(ours) == Chem.MolToSmiles(theirs), seq


@pytest.mark.parametrize("seq", ["RONIN", "RONIU", "OUO", "OU"])
def test_rare_residues_build_even_though_rdkit_rejects_them(seq):
    """Pyrrolysine and selenocysteine are why we do not use MolFromSequence."""
    assert Chem.MolFromSequence(seq) is None
    mol = Chem.MolFromSmiles(C.build_smiles(seq))
    assert mol is not None


def test_formula_matches_biopython_mass():
    """Formula and displayed weight must agree. The old code failed this."""
    for seq in ["RONIN", "ROOSA", "KLVFFAE", "GPGPG", "WYFAC"]:
        mol = Chem.MolFromSmiles(C.build_smiles(seq))
        from rdkit.Chem.Descriptors import MolWt
        assert MolWt(mol) == pytest.approx(C.analyse(seq)["molecular_weight"], abs=0.15), seq


def test_our_formula_arithmetic_matches_rdkit():
    """The app computes formulas without RDKit. RDKit is the check on that.

    This is the test that would have caught the old water-loss bug: it removed
    one hydrogen per peptide bond instead of two.
    """
    rng = random.Random(7)
    seqs = ["RONIN", "ROOSA", "ANNA", "GPGPG", "WYFAC", "AA", "M", "OU", "LAURA"]
    seqs += ["".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(rng.randint(2, 30)))
             for _ in range(40)]
    for seq in seqs:
        mol = Chem.MolFromSmiles(C.build_smiles(seq))
        assert mol is not None, seq
        assert C.molecular_formula(seq) == CalcMolFormula(mol), seq


def test_water_loss_is_two_hydrogens_per_bond():
    """Two glycines (C2H5NO2 each) minus one water is C4H8N2O3, not C4H9N2O4."""
    assert C.molecular_formula("G") == "C2H5NO2"
    assert C.molecular_formula("GG") == "C4H8N2O3"
    assert C.molecular_formula("GGG") == "C6H11N3O4"


def test_formula_never_writes_a_subscript_of_one():
    assert C.molecular_formula("C").endswith("S")
    assert "S1" not in C.molecular_formula("MMC")


def test_every_letter_of_the_alphabet_maps():
    for ch in string.ascii_uppercase:
        p = C.text_to_sequence(ch)
        assert len(p.sequence) == 1, ch
        assert p.sequence in C.RESIDUES


def test_substitutions_are_recorded_with_positions():
    p = C.text_to_sequence("Bob")
    assert p.sequence == "NON"
    assert [(s.position, s.original, s.replacement) for s in p.swaps] == [
        (1, "B", "N"), (3, "B", "N"),
    ]

    p = C.text_to_sequence("Jazz")
    assert p.sequence == "LAQQ"
    assert [(s.position, s.original) for s in p.swaps] == [(1, "J"), (3, "Z"), (4, "Z")]


def test_empty_and_junk_input_never_raises():
    for junk in ["", "123", "!!!", "   ", "🙂", None]:
        p = C.text_to_sequence(junk)
        assert p.sequence == ""
        assert C.analyse(p.sequence) == {}
        assert C.build_smiles(p.sequence) == ""


def test_accented_names_keep_their_letters():
    assert C.text_to_sequence("Zoë").sequence == "QOE"
    assert C.text_to_sequence("Håkan").sequence == "HAKAN"
    assert C.text_to_sequence("Jää").sequence == "LAA"
    assert C.text_to_sequence("Søren").sequence == "SOREN"


def test_length_is_capped():
    p = C.text_to_sequence("A" * 500)
    assert p.length == C.MAX_LENGTH
    assert p.truncated is True


def test_formula_html_omits_subscript_one():
    assert C.format_formula_html("C3H7NO2S") == (
        "C<sub>3</sub>H<sub>7</sub>NO<sub>2</sub>S"
    )


def test_gravy_handles_rare_residues():
    """Biopython's gravy() raises KeyError here; ours must not."""
    assert isinstance(C.gravy("RONIN"), float)
    assert C.gravy("IIII") == pytest.approx(4.5)


def test_proline_ring_closes_into_the_backbone():
    mol = Chem.MolFromSmiles(C.build_smiles("GPG"))
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(Chem.MolFromSequence("GPG"))


def test_stereochemistry_is_l_configuration():
    """All natural residues are L. A racemic or flat structure would be wrong."""
    mol = Chem.MolFromSmiles(C.build_smiles("AA"))
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    assert len(centres) == 2
    assert all(tag == "S" for _, tag in centres)
