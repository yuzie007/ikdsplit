"""Fetch Wyckoff from the Bilbao Crystallographic Server."""

import urllib

import pandas as pd

mapper = {
    "Multiplicity": "multiplicity",
    "Wyckoff  letter": "wyckoff_letter",
    "Site symmetry": "site_symmetry",
    "Coordinates": "coordinates",
}
for spg_no in range(1, 231):
    print(spg_no)
    URL = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list?"
    query = {"gnum": spg_no, "what": "wpos"}
    URL += urllib.parse.urlencode(query)
    df = pd.read_html(URL)[1]
    # when centering translations exist
    if isinstance(df.columns, pd.MultiIndex):
        df = pd.read_html(URL, skiprows=[1])[1]
        df.columns = df.columns.get_level_values(0)
    df = df.rename(mapper, axis=1)
    df.loc[:, "coordinates"] = df["coordinates"].str.split().str[0]
    df.to_csv(f"{spg_no:03d}.csv", index=False)
