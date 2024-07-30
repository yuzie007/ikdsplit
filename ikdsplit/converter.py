"""Convert `atoms_conventional.csv` for the subgroup."""

import argparse
import pathlib
import tomllib

import numpy as np
import pandas as pd
from ase.spacegroup import Spacegroup

import ikdsplit
from ikdsplit.spacegroup import (
    get_setting_for_origin_choice_2,
    parse_transformation,
)
from ikdsplit.utils import format_df


def expand_equivalent_positions(
    xyz0: np.ndarray,
    space_group_number: int,
) -> np.ndarray:
    """Expand symmetrically equivalent positions.

    Parameters
    ----------
    xyz0 : np.ndarray
        Original positions.
    space_group_number : int
        Space group number.

    Returns
    -------
    xyzs : np.ndarray
        Symmetrically equivalent positions.

    """
    setting = get_setting_for_origin_choice_2(space_group_number)
    symops = Spacegroup(space_group_number, setting).get_symop()
    xyzs = [xyz0]
    for r, t in symops:
        xyz = r @ xyz0 + t
        diff = xyz[None, :] - np.array(xyzs)
        diff -= np.rint(diff)
        if not np.any(np.all(np.isclose(diff, 0.0), axis=1)):
            xyzs.append(xyz)
    return np.array(xyzs)


def reduce_equivalent_positions(
    xyz0s: np.ndarray,
    space_group_number: int,
) -> np.ndarray:
    """Reduce symmetrically equivalent positions.

    Parameters
    ----------
    xyz0s : np.ndarray
        Positions to be reduced.
    space_group_number : int
        Space group number.

    Returns
    -------
    np.ndarray
        Reduced positions.

    """
    setting = get_setting_for_origin_choice_2(space_group_number)
    symops = Spacegroup(space_group_number, setting).get_symop()
    xyzs = [xyz0s[0]]
    for xyz0 in xyz0s:
        for r, t in symops:
            xyz = r @ xyz0 + t
            diff = xyz[None, :] - np.array(xyzs)
            diff -= np.rint(diff)
            if np.any(np.all(np.isclose(diff, 0.0), axis=1)):
                break
        else:
            xyzs.append(xyz0)
    return np.array(xyzs)


def find_wyckoff(xyz0: np.ndarray, space_group_number: int) -> str:
    """Find Wyckoff position.

    Parameters
    ----------
    xyz0 : np.ndarray
        Position to be checked.
    space_group_number : int
        Space group number.

    Returns
    -------
    str
        Wyckoff letter with multiplicity.

    """
    src = pathlib.Path(ikdsplit.__file__).parent / "database" / "wyckoff"
    wyckoffs = pd.read_csv(src / f"{space_group_number:03d}.csv")

    # get all symmetrically equivalent positions
    # one of them should be equal to the representative Wyckoff coordinates
    xyzs = expand_equivalent_positions(xyz0, space_group_number)

    for d in wyckoffs.to_dict(orient="records")[::-1]:
        # skip if #(symmetrically equivalent points) != multiplicity
        if len(xyzs) != d["multiplicity"]:
            continue
        for xyz in xyzs:
            r, t = parse_transformation(d["coordinates"])
            xyz_refined = r @ xyz + t
            diff = xyz - xyz_refined
            diff -= np.rint(diff)
            if np.allclose(diff, 0.0):
                return str(d["multiplicity"]) + d["wyckoff_letter"]
    raise RuntimeError(xyz0, space_group_number)


def convert() -> None:
    """Convert `atoms_conventional.csv` for the subgroup."""
    with pathlib.Path("wycksplit.toml").open("rb") as f:
        wycksplit = tomllib.load(f)

    spgno_sup = wycksplit["space_group_number_sup"]
    spgno_sub = wycksplit["space_group_number_sub"]

    df = pd.read_csv("../atoms_conventional.csv", skipinitialspace=True)

    ds = []
    for _, s in df.iterrows():
        xyz0 = s[["x", "y", "z"]].to_numpy(float)

        xyzs = expand_equivalent_positions(xyz0, spgno_sup)

        # change coordinates for subgroup
        xyzs -= wycksplit["origin_shift"]
        xyzs = (np.linalg.inv(wycksplit["basis_change"]) @ xyzs.T).T
        xyzs -= np.rint(xyzs)

        xyzs = reduce_equivalent_positions(xyzs, spgno_sub)

        for xyz in xyzs:
            d = {}
            d["symbol"] = s["symbol"]
            for k in s.index.to_numpy().tolist():
                if k.startswith("wyckoff_"):
                    d[k] = s[k]
            d[f"wyckoff_{spgno_sub:03d}"] = find_wyckoff(xyz, spgno_sub)
            d["x"], d["y"], d["z"] = xyz
            ds.append(d)

    df = pd.DataFrame(ds)
    df = format_df(df)
    filename = "atoms_conventional.csv"
    df.to_csv(filename, float_format="%24.18f", index=False)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    pass


def run(args: argparse.Namespace) -> None:
    convert()
