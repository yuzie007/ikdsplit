"""Fetch Wycksplit from the Bilbao Crystallographic Server."""

import time
import urllib

import pandas as pd


def fetch_transformation(spg_sup: int, spg_sub: int, index: int) -> None:
    """Get transformation."""
    url = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-tranmax?"
    query = {"super": spg_sup, "sub": spg_sub, "index": index}
    url += urllib.parse.urlencode(query)
    transformations = pd.read_html(url)[1]
    s = transformations.loc[0]["Transformation Matrix"]
    ts = s.split()
    d = {}
    d["supergroup"] = spg_sup
    d["subgroup"] = spg_sub
    d["Pxx"] = ts[0]
    d["Pxy"] = ts[1]
    d["Pxz"] = ts[2]
    d["Pyx"] = ts[4]
    d["Pyy"] = ts[5]
    d["Pyz"] = ts[6]
    d["Pzx"] = ts[8]
    d["Pzy"] = ts[9]
    d["Pzz"] = ts[10]
    d["px"] = ts[3]
    d["py"] = ts[7]
    d["pz"] = ts[11]
    return d


def run() -> None:
    """Run."""
    ds = []
    for supergroup in range(1, 231):
        print(f"{supergroup:03d}")
        url = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-lxi?"
        query = {"gnum": supergroup}
        url += urllib.parse.urlencode(query)
        subgroups = pd.read_html(url)[1]
        subgroups = subgroups[subgroups["Type"] == "t"]
        for d in subgroups.to_dict(orient="records"):
            subgroup = d["IT number"]
            index = d["Index"]
            print(f"-> {subgroup:03d}")
            ds.append(fetch_transformation(supergroup, subgroup, index))
            time.sleep(5)
    pd.DataFrame(ds).to_csv("transformations_t.csv", index=False)


if __name__ == "__main__":
    run()
