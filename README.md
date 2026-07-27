# Peptide Tag

Type your name. Get the peptide it spells, and find out where that sequence
turns up inside a real human protein.

Live at <https://peptide-tag-production.up.railway.app>

## The idea

Amino acids each have a one-letter code, and 22 letters of the alphabet are one
of them. So most names are already a peptide: a real molecule with a real mass
that a lab could synthesise. Then you can go looking for it, because the human
proteome is a searchable string of 11.4 million residues.

About a third of first names appear in a human protein whole. For the rest, the
app finds the longest run of the name that does, which means everybody gets an
answer and longer names get more interesting rather than more hopeless.

## What it will not do

Everything shown is computed or looked up. There is no invented biology, no
predicted function, and no therapeutic claim. A five-residue peptide spelled
from somebody's name is not a drug and the page says so.

The 3D view is labelled as what it is: the sequence built in two textbook
backbone conformations, an alpha helix and an extended strand. A peptide this
short has no fixed structure, and being able to flip between the two shapes is
the honest way to show that.

## Running it

```bash
python -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
.venv/bin/python -m pytest tests/ -q
.venv/bin/python app.py            # http://127.0.0.1:5055
```

## How it is checked

`tests/test_chemistry.py` rebuilds every peptide with RDKit and compares
structures and molecular formulas against ours, residue by residue, in several
sequence contexts, plus random sequences. RDKit is a development dependency
only; the app computes formulas itself so it does not have to ship a 128 MB
chemistry toolkit to serve a page.

## Layout

| file | what it does |
|---|---|
| `app.py` | routes |
| `chemistry.py` | residues, name to sequence, SMILES, formula, pI, hydropathy |
| `proteome.py` | searching reviewed human Swiss-Prot |
| `artwork.py` | the residue chain: inline SVG, printable card, social image |
| `structure.py` | ideal-geometry 3D model |
| `printing.py` | watertight STL for 3D printing |

## Data

`data/human_swissprot.fasta.gz` is 20,431 reviewed human proteins from UniProt,
committed so the app has no runtime network dependency. Refresh with:

```bash
python -m proteome --refresh
```

Sequence data from [UniProt](https://www.uniprot.org/), CC BY 4.0.
Fonts under the SIL Open Font License, see `static/fonts/README.md`.
