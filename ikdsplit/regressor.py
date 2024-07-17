import os

import ase.io
import numpy as np
import pandas as pd
import yaml
from ase import Atoms
from ase.build import make_supercell


def format_df(df):
    for k in df.columns:
        if df[k].dtype == object:
            df[k] = df[k].str.pad(12)
    return df


def change_coordinates(
    atoms: Atoms,
    df: pd.DataFrame,
    basis_change: np.ndarray,
    origin_shift: np.ndarray,
) -> tuple[Atoms, pd.DataFrame]:
    """Change the coordinate system."""
    atoms.set_scaled_positions(atoms.get_scaled_positions() - origin_shift)
    atoms = make_supercell(atoms, basis_change.T, order="atom-major")

    df[["x", "y", "z"]] -= origin_shift
    df[["x", "y", "z"]] @= basis_change.T

    return atoms, df


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


def add_arguments(parser):
    parser.add_argument(
        "--rotation",
        nargs=9,
        default=[1, 0, 0, 0, 1, 0, 0, 0, 1],
        type=int,
    )
    parser.add_argument(
        "--translation",
        nargs=3,
        default=[0.0, 0.0, 0.0],
        type=float,
    )


def run(args):
    atoms = ase.io.read("POSCAR")
    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)

    fn = ""
    while True:
        fn = "wycksplit.yaml" if not fn else os.path.join("..", fn)
        if not os.path.isfile(fn):
            break
        with open(fn) as f:
            d = yaml.safe_load(f)
        basis_change, origin_shift = invert(d["rotation"], d["translation"])
        atoms, df = change_coordinates(atoms, df, basis_change, origin_shift)

    basis_change = np.array(args.rotation).reshape(3, 3)
    origin_shift = np.array(args.translation)
    atoms, df = change_coordinates(atoms, df, basis_change, origin_shift)

    atoms.write("SPOSCAR_regressed", direct=True)

    df = format_df(df)
    df.to_csv("atoms_regressed.csv", float_format="%24.18f", index=False)
