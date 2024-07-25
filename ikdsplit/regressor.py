"""Regress to the original target supercell."""

import argparse
import functools
import os
import pathlib
import tomllib

import ase.io
import numpy as np
import pandas as pd
from ase import Atoms
from ase.build import make_supercell
from ase.spacegroup.crystal_data import _lattice_centering

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


def invert(
    basis_change: np.ndarray,
    origin_shift: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Invert a general change of the coordinate system.

    Parameters
    ----------
    basis_change : np.ndarray
        P
    origin_shift : np.ndarray
        p

    """
    inv_basis_change = np.linalg.inv(basis_change)
    inv_origin_shift = -1.0 * inv_basis_change @ origin_shift
    return inv_basis_change, inv_origin_shift


def multiply(
    op0: tuple[np.ndarray, np.ndarray],
    op1: tuple[np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Multiply two changes of the coordinate system.

    See Sec. (1.2.2.2) in ITA (2016).
    """
    basis_change_0, origin_shift_0 = op0
    basis_change_1, origin_shift_1 = op1
    basis_change = basis_change_0 @ basis_change_1
    origin_shift = basis_change_0 @ origin_shift_1 + origin_shift_0
    return basis_change, origin_shift


def cumulate_coordinate_change(
    transformations: list,
) -> tuple[np.ndarray, np.ndarray]:
    """Cumulate changes of the coordinate system in the application order."""
    ops = []

    wycksplit_toml = pathlib.Path("wycksplit.toml")
    if wycksplit_toml.is_file():
        with wycksplit_toml.open("rb") as f:
            d = tomllib.load(f)
        basis_change_first = get_p2c(d["space_group_number_sub"])
        origin_shift_first = np.array([0.0, 0.0, 0.0])
        ops.append((basis_change_first, origin_shift_first))

    fn = ""
    while True:
        fn = "wycksplit.toml" if not fn else os.path.join("..", fn)
        if not os.path.isfile(fn):
            break
        with pathlib.Path(fn).open("rb") as f:
            d = tomllib.load(f)
        ops.append(invert(d["basis_change"], d["origin_shift"]))

    if transformations:
        ops_last = [
            (np.array(_["basis_change"]), np.array(_["origin_shift"]))
            for _ in transformations
        ]
        ops.extend(ops_last)

    return functools.reduce(multiply, ops)


def regress(transformations: list) -> None:
    """Regress to the original target supercell."""
    basis_change, origin_shift = cumulate_coordinate_change(transformations)

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
        fin = f"PPOSCAR-{index:09d}"
        fout = f"RPOSCAR-{index:09d}"
        atoms = ase.io.read(fin)
        atoms = change_coordinates(atoms, basis_change, origin_shift)
        atoms.write(fout, direct=True)
        d.update({symbol: atoms.symbols.count(symbol) for symbol in symbols})
        ds.append(d)
    df = pd.DataFrame(ds)
    df.to_csv("info_regressed.csv", index=False)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments."""


def run(args: argparse.Namespace) -> None:
    """Run."""
    with pathlib.Path("ikdsplit.toml").open("rb") as f:
        d = tomllib.load(f)
    regress(d["regress"]["transformations"])
