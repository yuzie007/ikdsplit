import itertools
import tomllib

import numpy as np
import pandas as pd
from ase import Atoms
from ase.spacegroup import crystal, get_spacegroup


def index_wyckoff(df: pd.DataFrame) -> pd.DataFrame:
    """Index each site."""
    tmp = df.groupby(["wyckoff"]).cumcount() + 1
    df["wyckoff"] += tmp.astype(str)
    return df


def make_atoms(
    df: pd.DataFrame,
    spacegroup: int,
    cell: np.ndarray,
    primitive: bool,
):
    """Make `Atoms` based on `atoms_conventional.csv` and `fill`."""
    # hexagonal cell (consistent with Bilbao) for rhombohedral space groups
    setting = 1 if (3 <= spacegroup <= 15 or 143 <= spacegroup <= 194) else 2

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
    primitive: bool,
    mapping: dict[str, list[str]],
) -> tuple[list[Atoms], pd.DataFrame]:
    """Make all possible `Atoms`."""
    symbols = df["symbol"].unique()
    mappings = [mapping[_] for _ in df["fill"]]

    images = []
    ds = []
    i = -1
    filled: list[str]
    for filled in itertools.product(*mappings):
        i += 1
        df_included = df.copy()
        df_included["symbol"] = filled
        atoms = make_atoms(df_included, spacegroup, cell, primitive)
        images.append(atoms)
        d = {}
        d["configuration"] = i
        d["space_group_number"] = get_spacegroup(atoms).todict()["number"]
        d.update({symbol: atoms.symbols.count(symbol) for symbol in symbols})
        d.update(dict(zip(df["wyckoff"], filled, strict=True)))
        ds.append(d)

    return images, pd.DataFrame(ds)


def fill(mapping: dict[str, list[str]]) -> None:
    cell = np.loadtxt("cell.dat")

    with open("wycksplit.toml", "rb") as f:
        wycksplit = tomllib.load(f)

    spacegroup = wycksplit["space_group_number_sub"]

    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)
    df = index_wyckoff(df)

    for primitive in [False, True]:
        images, df_tmp = make_images(
            df,
            spacegroup,
            cell,
            primitive=primitive,
            mapping=mapping,
        )
        for i, atoms in enumerate(images):
            fn = f"PPOSCAR-{i:09d}" if primitive else f"CPOSCAR-{i:09d}"
            atoms.write(fn, direct=True)
        fn = "info_primitive.csv" if primitive else "info_conventional.csv"
        df_tmp.to_csv(fn, index=False)


def add_arguments(parser):
    pass


def run(args):
    with open("ikdsplit.toml", "rb") as f:
        d = tomllib.load(f)
    mapping = d["fill"]
    fill(mapping)
