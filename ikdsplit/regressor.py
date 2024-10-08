"""Regress to the original target supercell."""

import argparse
import functools
import pathlib
import tomllib

import ase.io
import numpy as np
import pandas as pd
from ase import Atoms
from ase.build import make_supercell
from ase.spacegroup.crystal_data import _lattice_centering

from ikdsplit.io import parse_config
from ikdsplit.runner import start_recursive
from ikdsplit.spacegroup import multiply
from ikdsplit.splitter import add_arguments
from ikdsplit.utils import format_df


def get_p2c(spacegroup: int) -> np.ndarray:
    """Get basis change from primitive to conventional."""
    centering = _lattice_centering[spacegroup]
    p2c = {
        "P": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        "C": np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]]),
        "F": np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]),
        "I": np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]),
        "R": np.array([[1, 0, 1], [-1, 1, 1], [0, -1, 1]]),
    }[centering]
    return np.array(p2c)


def change_coordinates(
    atoms: Atoms,
    basis_change: np.ndarray,
    origin_shift: np.ndarray,
) -> tuple[Atoms, pd.DataFrame]:
    """Change the coordinate system."""
    atoms.set_scaled_positions(atoms.get_scaled_positions() - origin_shift)
    return make_supercell(atoms, basis_change.T, order="atom-major")


def cumulate_coordinate_change() -> tuple[np.ndarray, np.ndarray]:
    """Cumulate changes of the coordinate system in the application order."""
    ops = []

    with pathlib.Path("ikdsplit.toml").open("rb") as f:
        config = tomllib.load(f)

    basis_change_first = get_p2c(config["space_group_number"])
    origin_shift_first = np.array([0.0, 0.0, 0.0])
    ops.append((basis_change_first, origin_shift_first))

    config = parse_config()

    ops_last = [
        (np.array(_["basis_change"]), np.array(_["origin_shift"]))
        for _ in config["regress"]["transformations"]
    ]
    ops.extend(ops_last)

    return functools.reduce(multiply, ops)


def regress() -> None:
    """Regress to the original target supercell."""
    basis_change, origin_shift = cumulate_coordinate_change()
    if not np.allclose(basis_change, np.rint(basis_change)):
        print("regressed cell not compatible ... skipped", end=" ")
        return

    # calculate atomic positions in the target supercell
    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)
    symbols = df["symbol"].unique()
    df[["x", "y", "z"]] -= origin_shift
    df[["x", "y", "z"]] @= np.linalg.inv(basis_change).T
    df = format_df(df)
    df.to_csv("atoms_regressed.csv", float_format="%24.18f", index=False)

    df = pd.read_csv("info_conventional.csv", skipinitialspace=True)
    ds = []
    for d in df.to_dict(orient="records"):
        index = d["configuration"]
        fin = f"PPOSCAR-{index:06d}"
        fout = f"RPOSCAR-{index:06d}"
        atoms = ase.io.read(fin)
        atoms = change_coordinates(atoms, basis_change, origin_shift)
        atoms.write(fout, direct=True)
        d.update({symbol: atoms.symbols.count(symbol) for symbol in symbols})
        ds.append(d)
    df = pd.DataFrame(ds)
    df.to_csv("info_regressed.csv", index=False)


def run(args: argparse.Namespace) -> None:
    """Run."""
    criteria = {
        "max_level": args.level,
        "min_order": args.order,
        "max_configurations": args.configurations,
    }
    if args.recursive:
        start_recursive(regress, criteria)
    else:
        regress()


def main() -> None:
    """Run as a script."""
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter_class)
    add_arguments(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
