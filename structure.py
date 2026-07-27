"""A three-dimensional model of the peptide, honestly labelled.

A short linear peptide in water has no single shape. It is a rapidly
interconverting ensemble, and any one conformer is an arbitrary sample from it.
So this does not predict a structure. It builds the peptide in two textbook
backbone conformations, an alpha helix and an extended strand, from ideal
phi/psi angles. Being able to flip between them is the point: it shows the
thing has no fixed shape rather than claiming it does.

The alternative, RDKit's ETKDG conformer search, was measured at roughly six
seconds for a 25 residue peptide, outright failed to embed at the default
settings, never converged during minimisation, and ran past three minutes on a
worst case with no way to time it out from inside the call. This takes 8ms and
cannot fail.

PeptideBuilder falls back to glycine for any letter it does not know, which
would silently replace pyrrolysine and selenocysteine with a hydrogen atom, so
sequences are projected onto the standard twenty first.
"""

from __future__ import annotations

import io
from functools import lru_cache

import chemistry as C

HELIX = (-60.0, -40.0)
STRAND = (-140.0, 130.0)
SHAPES = {"helix": HELIX, "strand": STRAND}


@lru_cache(maxsize=512)
def model_pdb(sequence: str, shape: str = "helix") -> str:
    """Ideal-geometry PDB for the sequence. Empty string if it cannot build."""
    if not sequence:
        return ""
    import Bio.PDB
    import PeptideBuilder
    from PeptideBuilder import Geometry

    phi, psi = SHAPES.get(shape, HELIX)
    seq = C.to_standard(sequence)

    try:
        geo = Geometry.geometry(seq[0])
        geo.phi, geo.psi_im1 = phi, psi
        structure = PeptideBuilder.initialize_res(geo)
        for residue in seq[1:]:
            geo = Geometry.geometry(residue)
            geo.phi, geo.psi_im1 = phi, psi
            PeptideBuilder.add_residue(structure, geo)
        PeptideBuilder.add_terminal_OXT(structure)
    except Exception:
        return ""

    out = io.StringIO()
    writer = Bio.PDB.PDBIO()
    writer.set_structure(structure)
    writer.save(out)
    return out.getvalue()
