"""Sort atoms."""

import argparse
import pathlib
import shutil
import tomllib

import ase.io
import numpy as np
import pandas as pd
from ase import Atoms

from ikdsplit.runner import start_recursive
from ikdsplit.splitter import add_arguments


def get_reference() -> str:
    """Get reference.

    Parameters
    ----------
    reference : str | None, default = None
        Atoms with the reference positions.
        If None, each file is simply copied.

    """
    with pathlib.Path("ikdsplit.toml").open("rb") as f:
        config = tomllib.load(f)
    return config["sort"]["reference"]


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


def sort() -> None:
    """Sort atoms in the order in the reference file.

    Only `scaled_positions` are referred to, and `cell` is kept.

    Original atomic positions should be stored in `SPOSCAR_regressed`.
    """
    reference = get_reference()
    atoms_ref = ase.io.read(reference)
    df = pd.read_csv("info_conventional.csv", skipinitialspace=True)
    for d in df.to_dict(orient="records"):
        index = d["configuration"]
        fin = f"RPOSCAR-{index:06d}"
        fout = f"SPOSCAR-{index:06d}"
        if not pathlib.Path(fin).is_file():
            continue
        if reference is None:
            shutil.copy2(fin, fout)
        else:
            atoms = ase.io.read(fin)
            atoms = sort_atoms(atoms, atoms_ref)
            atoms.write(fout, direct=True)


def run(args: argparse.Namespace) -> None:
    """Run."""
    criteria = {
        "max_level": args.level,
        "min_order": args.order,
        "max_configurations": args.configurations,
    }
    start_recursive(sort, criteria)


def main() -> None:
    """Run as a script."""
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter_class)
    add_arguments(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
