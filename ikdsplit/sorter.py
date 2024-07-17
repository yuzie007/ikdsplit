import ase.io
import numpy as np


def add_arguments(parser):
    parser.add_argument(
        "--ref",
        help="atoms with the reference order",
    )


def run(args):
    """Sort atoms in the order in the reference file.

    Only `scaled_positions` are referred to, and `cell` is kept.

    Original atomic positions should be stored in `SPOSCAR_regressed`.

    Parameters
    ----------
    args : _type_
        _description_

    """
    atoms = ase.io.read("SPOSCAR_regressed")
    atoms_ref = ase.io.read(args.ref)

    tmp = []
    for i in range(len(atoms)):
        sp_ref = atoms.get_scaled_positions()[i]
        diffs = atoms_ref.get_scaled_positions() - sp_ref
        diffs -= np.round(diffs)
        indices = np.where(np.all(np.abs(diffs) < 1e-12, axis=1))[0]
        assert len(indices) in [0, 1]
        if len(indices) == 1:
            tmp.append(indices[0])
    tmp = sorted(tmp)
    atoms_sym = atoms_ref[tmp]
    atoms_sym.write("SPOSCAR_symmetry", direct=True)
    atoms_met = atoms_sym[[atom.index for atom in atoms if atom.symbol != "H"]]
    atoms_met.write("SPOSCAR", direct=True)
