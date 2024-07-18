import functools
import os

import ase.io
import numpy as np
import pandas as pd
import yaml
from ase import Atoms
from ase.build import make_supercell

from ikdsplit.utils import format_df


def change_coordinates(
    atoms: Atoms,
    basis_change: np.ndarray,
    origin_shift: np.ndarray,
) -> tuple[Atoms, pd.DataFrame]:
    """Change the coordinate system."""
    atoms.set_scaled_positions(atoms.get_scaled_positions() - origin_shift)
    atoms = make_supercell(atoms, basis_change.T, order="atom-major")
    return atoms


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


def cumulate_coordinate_change(ops_last) -> tuple[np.ndarray, np.ndarray]:
    """Cumulate changes of the coordinate system in the application order."""
    ops = []
    fn = ""
    while True:
        fn = "wycksplit.yaml" if not fn else os.path.join("..", fn)
        if not os.path.isfile(fn):
            break
        with open(fn, encoding="utf-8") as f:
            d = yaml.safe_load(f)
        ops.append(invert(d["rotation"], d["translation"]))
    ops.extend(ops_last)
    return functools.reduce(multiply, ops)


def add_arguments(parser):
    parser.add_argument("--transformations")


def run(args):
    ops = []
    if args.transformations is not None:
        with open(args.transformations, encoding="utf-8") as f:
            ops = yaml.safe_load(f)
        ops = [[np.array(_) for _ in op] for op in ops]
    basis_change, origin_shift = cumulate_coordinate_change(ops)

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
        index = d["index"]
        fin = f"POSCAR-{index:09d}"
        fout = f"RPOSCAR-{index:09d}"
        atoms = ase.io.read(fin)
        atoms = change_coordinates(atoms, basis_change, origin_shift)
        atoms.write(fout, direct=True)
        d.update({symbol: atoms.symbols.count(symbol) for symbol in symbols})
        ds.append(d)
    df = pd.DataFrame(ds)
    df.to_csv("info_regressed.csv", index=False)
