import shutil
import os

import numpy as np


def add_arguments(parser):
    parser.add_argument("group", type=int)


def run(args):
    group = args.group
    src = os.path.join(os.path.dirname(__file__), "database")
    fn = os.path.join(src, "subgroups.dat")
    subgroups = np.genfromtxt(fn, dtype=(int, int, "U1"), usecols=(0, 1, 2))
    g2h = subgroups[[_[0] == group for _ in subgroups]]
    subgroups = [_[1] for _ in g2h]
    for subgroup in subgroups:
        print(f"{group:03d} -> {subgroup:03d}")
        dn = f"{subgroup:03d}"
        os.makedirs(dn, exist_ok=True)
        fn = os.path.join(src, "wycksplit", f"{group:03d}_{subgroup:03d}.yaml")
        shutil.copy(fn, os.path.join(dn, "wycksplit.yaml"))
