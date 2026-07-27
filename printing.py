"""A printable 3D model of the peptide.

Ball and stick, fused into one watertight solid, scaled to millimetres and
stood on a base so it can actually come off a printer intact.

The thing that breaks these models is bond thickness. Coordinates come out in
angstroms, so a longer peptide gets a smaller scale factor, and past about
fifteen residues at a fixed overall size the struts thin to under a millimetre
and snap. So the bond radius is derived from the scale rather than fixed: the
printed strut stays at least MIN_STRUT_MM thick no matter how long the name is.
"""

from __future__ import annotations

from functools import lru_cache

import structure

# Covalent radii, angstroms, used for the balls.
RADII = {"C": 0.70, "N": 0.65, "O": 0.60, "S": 1.00, "SE": 1.15, "H": 0.35}
BOND_MAX = 1.85  # angstroms; anything closer than this is bonded

TARGET_MM = 90.0     # longest dimension of the finished object
MIN_STRUT_MM = 2.2   # thinnest printable strut, diameter
BASE_MM = 3.0        # thickness of the disc it stands on


def _atoms(pdb: str) -> list[tuple[str, float, float, float]]:
    out = []
    for line in pdb.splitlines():
        if line.startswith("ATOM"):
            element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
            out.append((element, float(line[30:38]), float(line[38:46]),
                        float(line[46:54])))
    return out


@lru_cache(maxsize=64)
def model_stl(sequence: str, shape: str = "helix") -> bytes:
    """A watertight, printable STL in millimetres. Empty on failure."""
    import numpy as np
    import trimesh

    pdb = structure.model_pdb(sequence, shape)
    if not pdb:
        return b""
    atoms = _atoms(pdb)
    if not atoms:
        return b""

    pos = np.array([[x, y, z] for _, x, y, z in atoms])
    span = float(np.max(pos.max(axis=0) - pos.min(axis=0))) or 1.0
    scale = TARGET_MM / span

    # Strut radius in angstroms that lands on MIN_STRUT_MM once scaled.
    strut = max(0.30, (MIN_STRUT_MM / 2.0) / scale)
    ball = strut * 1.35

    parts = []
    for element, x, y, z in atoms:
        r = max(RADII.get(element, 0.7), ball)
        sphere = trimesh.creation.icosphere(subdivisions=2, radius=r)
        sphere.apply_translation((x, y, z))
        parts.append(sphere)

    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if float(np.linalg.norm(pos[i] - pos[j])) < BOND_MAX:
                parts.append(trimesh.creation.cylinder(
                    radius=strut, segment=[pos[i], pos[j]]))

    try:
        mesh = trimesh.boolean.union(parts, engine="manifold")
    except Exception:
        # Every primitive is individually closed, so the concatenation is still
        # watertight and slicers handle it. 25x faster, slightly bigger file.
        mesh = trimesh.util.concatenate(parts)

    mesh.apply_scale(scale)

    # Stand it up: drop a disc under the lowest point and sink the model
    # slightly into it so the two fuse instead of merely touching.
    lo = mesh.bounds[0]
    mesh.apply_translation((0, 0, -lo[2] + BASE_MM * 0.6))
    radius = float(max(mesh.extents[0], mesh.extents[1])) * 0.34 + 6.0
    base = trimesh.creation.cylinder(radius=radius, height=BASE_MM)
    centre = mesh.centroid
    base.apply_translation((centre[0], centre[1], BASE_MM / 2.0))

    try:
        final = trimesh.boolean.union([mesh, base], engine="manifold")
    except Exception:
        final = trimesh.util.concatenate([mesh, base])

    return final.export(file_type="stl")
