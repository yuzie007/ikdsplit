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
    symprec = 1e-12
    tmp = []
    for i in range(len(atoms)):
        sp_ref = atoms.get_scaled_positions()[i]
        diffs = atoms_ref.get_scaled_positions() - sp_ref
        diffs -= np.round(diffs)
        indices = np.where(np.all(np.abs(diffs) < symprec, axis=1))[0]
        assert len(indices) in [0, 1]
        if len(indices) == 1:
            tmp.append(indices[0])
    tmp = sorted(tmp)
    return atoms_ref[tmp]


def add_arguments(parser):
    parser.add_argument(
        "--ref",
        help="atoms with the reference order",
    )


def run(args):
    """Sort atoms in the order in the reference file.

    Only `scaled_positions` are referred to, and `cell` is kept.

    Original atomic positions should be stored in `SPOSCAR_regressed`.

    Parameters
    ----------
    args : _type_
        _description_

    """
    atoms_ref = ase.io.read(args.ref)
    df = pd.read_csv("info_conventional.csv", skipinitialspace=True)
    for d in df.to_dict(orient="records"):
        index = d["index"]
        fin = f"RPOSCAR-{index:09d}"
        fout = f"SPOSCAR-{index:09d}"
        atoms = ase.io.read(fin)
        atoms = sort_atoms(atoms, atoms_ref)
        atoms.write(fout, direct=True)
