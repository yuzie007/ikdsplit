import itertools
import tomllib

import numpy as np
import pandas as pd
from ase.spacegroup import crystal, get_spacegroup
from ase import Atoms


def add_labels(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.groupby(["symbol", "wyckoff"]).cumcount() + 1
    df["label"] = df["symbol"] + "_" + df["wyckoff"] + tmp.astype(str)
    return df


def make_atoms(
    df: pd.DataFrame,
    spacegroup: int,
    cell: np.ndarray,
    primitive: bool,
):
    # hexagonal cell (consistent with Bilbao) for rhombohedral space groups
    setting = 1 if (3 <= spacegroup <= 15 or 143 <= spacegroup <= 194) else 2

    symbols = df["symbol"].str.strip().values
    basis = df[["x", "y", "z"]]
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
    choices: list[list[bool]],
) -> list[Atoms]:
    symbols = df["symbol"].unique()
    images = []
    ds = []
    i = -1
    for included in itertools.product(*choices):
        i += 1
        df_included = df[list(included)]
        atoms = make_atoms(df_included, spacegroup, cell, primitive)
        images.append(atoms)
        d = {}
        d["index"] = i
        d["space_group_number"] = get_spacegroup(atoms).todict()["number"]
        d.update({symbol: atoms.symbols.count(symbol) for symbol in symbols})
        d.update(dict(zip(df["label"], included, strict=True)))
        ds.append(d)

    return images, pd.DataFrame(ds)


def make_choices(
    df: pd.DataFrame,
    always: list[str],
    never: list[str],
    selected: list[str],
) -> list[list[bool]]:
    choices = []
    for _, row in df.iterrows():
        if row["symbol"] in always:
            choices.append([True])
        elif row["symbol"] in never:
            choices.append([False])
        elif row["symbol"] in selected:
            choices.append([True, False])
        else:
            raise RuntimeError
    return choices


def fill(always: list[str], never: list[str], selected: list[str]):
    cell = np.loadtxt("cell.dat")

    with open("wycksplit.toml", "rb") as f:
        mapping = tomllib.load(f)

    spacegroup = mapping["space_group_number_sub"]

    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)
    df = add_labels(df)

    choices = make_choices(df, always, never, selected)

    for primitive in [False, True]:
        images, df_tmp = make_images(
            df,
            spacegroup,
            cell,
            primitive=primitive,
            choices=choices,
        )
        for i, atoms in enumerate(images):
            fn = f"PPOSCAR-{i:09d}" if primitive else f"CPOSCAR-{i:09d}"
            atoms.write(fn, direct=True)
        fn = "info_primitive.csv" if primitive else "info_conventional.csv"
        df_tmp.to_csv(fn, index=False)


def add_arguments(parser):
    parser.add_argument(
        "--always",
        nargs="+",
        default=[],
        help="symbols always included",
    )
    parser.add_argument(
        "--never",
        nargs="+",
        default=[],
        help="symbols never included",
    )
    parser.add_argument(
        "--selected",
        nargs="+",
        default=[],
        help="symbols selected",
    )


def run(args):
    fill(args.always, args.never, args.selected)
