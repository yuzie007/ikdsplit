"""Convert `atoms_conventional.csv` for the subgroup."""

import argparse
import pathlib
import tomllib

import numpy as np
import pandas as pd

from ikdsplit.spacegroup import (
    expand_equivalent_positions,
    find_wyckoff,
    reduce_equivalent_positions,
)
from ikdsplit.utils import format_df


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
