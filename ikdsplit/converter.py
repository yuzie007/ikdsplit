"""Convert `atoms_conventional.csv` for the subgroup."""

import argparse

import numpy as np
import pandas as pd

from ikdsplit.io import fetch_transformation
from ikdsplit.spacegroup import (
    expand_equivalent_positions,
    find_wyckoff,
    reduce_equivalent_positions,
)
from ikdsplit.utils import format_df


def convert(spg_sup: int, spg_sub: int) -> None:
    """Convert `atoms_conventional.csv` for the subgroup."""
    change_of_basis, origin_shift = fetch_transformation(spg_sup, spg_sub)

    df = pd.read_csv("../atoms_conventional.csv", skipinitialspace=True)

    ds = []
    for _, s in df.iterrows():
        xyz0 = s[["x", "y", "z"]].to_numpy(float)

        xyzs = expand_equivalent_positions(xyz0, spg_sup)

        # change coordinates for subgroup
        xyzs -= origin_shift
        xyzs = (np.linalg.inv(change_of_basis) @ xyzs.T).T
        xyzs -= np.rint(xyzs)

        xyzs = reduce_equivalent_positions(xyzs, spg_sub)

        for xyz in xyzs:
            d = {}
            d["symbol"] = s["symbol"]
            for k in s.index.to_numpy().tolist():
                if k.startswith("wyckoff_"):
                    d[k] = s[k]
            d[f"wyckoff_{spg_sub:03d}"] = find_wyckoff(xyz, spg_sub)
            d["x"], d["y"], d["z"] = xyz
            ds.append(d)

    df = pd.DataFrame(ds)
    df = format_df(df)
    filename = "atoms_conventional.csv"
    df.to_csv(filename, float_format="%24.18f", index=False)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments."""
    parser.add_argument("--sup", type=int)
    parser.add_argument("--sub", type=int)


def run(args: argparse.Namespace) -> None:
    """Run."""
    convert(args.sup, args.sub)
