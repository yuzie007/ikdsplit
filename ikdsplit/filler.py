"""Fill."""

import argparse
import itertools

import numpy as np
import pandas as pd
from ase import Atoms
from ase.spacegroup import crystal, get_spacegroup

from ikdsplit.io import parse_config
from ikdsplit.spacegroup import get_setting_for_origin_choice_2


def index_wyckoff(df: pd.DataFrame, space_group_number: int) -> pd.DataFrame:
    """Index each site."""
    key = f"wyckoff_{space_group_number:03d}"
    tmp = df.groupby([key]).cumcount() + 1
    df[key] += tmp.astype(str)
    return df


def make_atoms(
    df: pd.DataFrame,
    spacegroup: int,
    cell: np.ndarray,
    *,
    primitive: bool = False,
) -> Atoms:
    """Make `Atoms` based on `atoms_conventional.csv` and `fill`."""
    # hexagonal cell (consistent with Bilbao) for rhombohedral space groups
    setting = get_setting_for_origin_choice_2(spacegroup)

    included = df["symbol"] != "X"  # remove vacancies
    symbols = df[included]["symbol"].str.strip().to_numpy()
    basis = df[included][["x", "y", "z"]]
    try:
        atoms = crystal(
            symbols,
            basis=basis,
            spacegroup=spacegroup,
            cell=cell,
            setting=setting,
            primitive_cell=primitive,
        )
    except Exception:
        atoms = crystal(
            symbols,
            basis=basis,
            spacegroup=spacegroup,
            cell=cell,
            primitive_cell=primitive,
        )
    return atoms


def make_images(
    df: pd.DataFrame,
    spacegroup: int,
    cell: np.ndarray,
    mapping: dict[str, list[str]],
    *,
    primitive: bool = False,
) -> pd.DataFrame:
    """Make all possible `Atoms`."""
    symbols = df["symbol"].unique()
    mappings = [mapping[_] for _ in df["symbol"]]

    ds = []
    filled: list[str]
    for i, filled in enumerate(itertools.product(*mappings)):
        df_included = df.copy()
        df_included["symbol"] = filled
        atoms = make_atoms(df_included, spacegroup, cell, primitive=primitive)

        fn = f"PPOSCAR-{i:09d}" if primitive else f"CPOSCAR-{i:09d}"
        atoms.write(fn, direct=True)

        d = {}
        d["configuration"] = i
        d["space_group_number"] = get_spacegroup(atoms).todict()["number"]
        d.update({symbol: atoms.symbols.count(symbol) for symbol in symbols})

        key = f"wyckoff_{spacegroup:03d}"
        d.update(dict(zip(df[key], filled, strict=True)))

        ds.append(d)

    return pd.DataFrame(ds)


def fill() -> None:
    """Fill atoms acoording to `atoms_conventional.csv`."""
    config = parse_config()

    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)
    df = index_wyckoff(df, config["space_group_number"])

    for primitive in [False, True]:
        info = make_images(
            df,
            config["space_group_number"],
            config["cell"],
            mapping=config["fill"],
            primitive=primitive,
        )
        fn = "info_primitive.csv" if primitive else "info_conventional.csv"
        info.to_csv(fn, index=False)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments."""


def run(args: argparse.Namespace) -> None:
    """Run."""
    fill()
