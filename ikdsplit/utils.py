"""Utilities."""

import contextlib
import os
import pathlib
import typing

import pandas as pd


@contextlib.contextmanager
def cd(path: str | pathlib.Path) -> typing.Generator:
    """Change directory temporalily.

    Parameters
    ----------
    path: Path
        Path to directory.

    """
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


def format_df(df: pd.DataFrame) -> pd.DataFrame:
    for k in df.columns:
        if df[k].dtype == object:
            df[k] = df[k].str.pad(12)
    return df
