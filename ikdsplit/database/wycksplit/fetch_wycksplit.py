"""Fetch Wycksplit from the Bilbao Crystallographic Server."""

import fractions
import pathlib
import urllib

import numpy as np
import pandas as pd


def write_toml(
    spg_sup: int,
    spg_sub: int,
    change_of_basis: list,
    origin_shift: list,
) -> None:
    """Write `wycksplit.toml`."""
    fn = f"{spg_sup:03d}_{spg_sub:03d}.toml"
    with pathlib.Path(fn).open("w", encoding="utf-8") as f:
        f.write(f"space_group_number_sup = {spg_sup}\n")
        f.write(f"space_group_number_sub = {spg_sub}\n")
        r = repr(change_of_basis)
        f.write(f"basis_change = {r}\n")
        r = repr(origin_shift)
        f.write(f"origin_shift = {r}\n")


def parse_transformation(s: str) -> tuple[np.ndarray, np.ndarray]:
    """Parse transformation."""
    values = [float(fractions.Fraction(_)) for _ in s.split()]
    change_of_basis = [
        [values[0], values[1], values[2]],
        [values[4], values[5], values[6]],
        [values[8], values[9], values[10]],
    ]
    origin_shift = [values[3], values[8], values[11]]
    return change_of_basis, origin_shift


def fetch_transformation(spg_sup: int, spg_sub: int, index: int) -> None:
    """Get transformation."""
    url = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-tranmax?"
    query = {"super": spg_sup, "sub": spg_sub, "index": index}
    url += urllib.parse.urlencode(query)
    transformations = pd.read_html(url)[1]
    s = transformations.loc[0]["Transformation Matrix"]
    change_of_basis, origin_shift = parse_transformation(s)
    write_toml(spg_sup, spg_sub, change_of_basis, origin_shift)


def run() -> None:
    """Run."""
    for spg_sup in range(230, 231):
        print(spg_sup)
        url = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-lxi?"
        query = {"gnum": spg_sup}
        url += urllib.parse.urlencode(query)
        subgroups = pd.read_html(url)[1]
        subgroups = subgroups[subgroups["Type"] == "t"]
        for d in subgroups.to_dict(orient="records"):
            spg_sub = d["IT number"]
            index = d["Index"]
            print(f"-> {spg_sub:03d}")
            fetch_transformation(spg_sup, spg_sub, index)


if __name__ == "__main__":
    run()
