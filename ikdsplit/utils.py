"""Utilities."""

import contextlib
import os
import pathlib
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
