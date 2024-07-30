"""Utilities."""

import contextlib
import math
import os
import pathlib
import tomllib
import typing

import pandas as pd

import ikdsplit


@contextlib.contextmanager
def cd(path: str | pathlib.Path) -> typing.Generator:
    """Change directory temporalily.

    Parameters
    ----------
    path: Path
        Path to directory.

    """
    cwd = pathlib.Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


def format_df(df: pd.DataFrame) -> pd.DataFrame:
    """Format `atoms_XXX.csv`."""
    for k in df.columns:
        if df[k].dtype == object:
            df[k] = df[k].str.pad(6)
    return df


def get_subgroups(group: int) -> list[int]:
    """Get subgroup numbers of the given space group number.

    Parameters
    ----------
    group : int
        Space group number.

    Returns
    -------
    list[int]
        List of the numbers of subgroup.

    """
    src = pathlib.Path(ikdsplit.__file__).parent / "database" / "wycksplit"
    fns = os.listdir(src)
    fns = [_ for _ in fns if _.endswith(".toml")]
    fns = [_ for _ in fns if _.startswith(f"{group:03d}")]
    return sorted([int(fn.split(".")[0].split("_")[1]) for fn in fns])


def print_group(group: int, level: int, *args: tuple, **kwargs: dict) -> None:
    """Print group with indent.

    Parameters
    ----------
    group : int
        Space group number.
    level : int
        Level of indent.
    *args
        Positional arguments passed to `print`.
    **kwargs
        Keyword arguments passed to `print`.

    """
    s = ""
    if level > 1:
        s += "   " * (level - 1)
    if level > 0:
        s += "-> "
    s += f"{group:03d}"
    print(s, *args, **kwargs)


def count_configurations() -> int:
    """Count number of atomic configurations."""
    with pathlib.Path("ikdsplit.toml").open("rb") as f:
        config = tomllib.load(f)
    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)
    return math.prod([len(config["fill"][_]) for _ in df["symbol"]])
