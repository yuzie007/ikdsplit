import pathlib
import shutil

import numpy as np
import tomllib

from ikdsplit.converter import convert
from ikdsplit.filler import fill
from ikdsplit.regressor import regress
from ikdsplit.sorter import sort_all
from ikdsplit.utils import cd


def run_each(config: dict):
    convert()

    d = config["fill"]
    fill(d["always"], d["never"], d["selected"])

    d = config["regress"]
    regress(d["transformations"])

    d = config["sort"]
    sort_all(d["reference"])


def run_all():
    with open("ikdsplit.toml", "rb") as f:
        config = tomllib.load(f)

    group = config["space_group_number"]

    src = pathlib.Path(__file__).parent / "database"
    fn = src / "subgroups.dat"
    subgroups = np.genfromtxt(fn, dtype=(int, int, "U1"), usecols=(0, 1, 2))
    g2h = subgroups[[_[0] == group for _ in subgroups]]
    subgroups = [_[1] for _ in g2h]
    for subgroup in subgroups:
        print(f"{group:03d} -> {subgroup:03d}")
        dn = pathlib.Path(f"{subgroup:03d}")
        dn.mkdir(parents=True, exist_ok=True)
        with cd(dn):
            fn = src / "wycksplit" / f"{group:03d}_{subgroup:03d}.yaml"
            shutil.copy(fn, "wycksplit.yaml")
            run_each(config)


def add_arguments(parser):
    pass


def run(args):
    run_all()
