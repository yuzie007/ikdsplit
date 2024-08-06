"""Sort atoms."""

import shutil

import ase.io
import numpy as np
import pandas as pd
from ase import Atoms


def sort_atoms(atoms: Atoms, atoms_ref: Atoms) -> Atoms:
    """Sort atoms.

    Parameters
    ----------
    atoms : Atoms
        Atoms to be sorted.
    atoms_ref : Atoms
        Atoms with the reference order.

    """
    scaled_positions_ref = atoms_ref.get_scaled_positions()
    scaled_positions = atoms.get_scaled_positions()

    symprec = 1e-12
    tmp = []
    for i in range(len(atoms)):
        diffs = scaled_positions[i] - scaled_positions_ref
        diffs -= np.round(diffs)
        indices = np.where(np.all(np.abs(diffs) < symprec, axis=1))[0]
        assert len(indices) in [0, 1]
        if len(indices) == 1:
            tmp.append(indices[0])
    assert len(tmp) == len(atoms)
    tmp = sorted(tmp)
    return atoms_ref[tmp]


def sort_all(reference: str | None) -> None:
    """Sort atoms in the order in the reference file.

    Only `scaled_positions` are referred to, and `cell` is kept.

    Original atomic positions should be stored in `SPOSCAR_regressed`.

    Parameters
    ----------
    reference : str | None, default = None
        Atoms with the reference positions.
        If None, each file is simply copied.

    """
    if reference is not None:
        atoms_ref = ase.io.read(reference)
    df = pd.read_csv("info_conventional.csv", skipinitialspace=True)
    for d in df.to_dict(orient="records"):
        index = d["configuration"]
        fin = f"RPOSCAR-{index:06d}"
        fout = f"SPOSCAR-{index:06d}"
        if reference is not None:
            atoms = ase.io.read(fin)
            atoms = sort_atoms(atoms, atoms_ref)
            atoms.write(fout, direct=True)
        else:
            shutil.copy2(fin, fout)


def add_arguments(parser):
    parser.add_argument(
        "--reference",
        help="reference atoms with the reference order",
    )


def run(args):
    sort_all(args.reference)
