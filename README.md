# ikdsplit

## Installation

Published version:

```bash
pip install git+https://github.com/yuzie007/ikdsplit.git@0.2.0
```

Development version:

```bash
git clone https://github.com/yuzie007/ikdsplit.git
cd ikdsplit
pip install -e .
```

## Input files

3 input files are necessary.

`atoms_conventional.csv`

```csv
symbol,wyckoff,x,y,z
    Ti,    8a,    0.125000000000000000,    0.125000000000000000,    0.125000000000000000
    Cr,   16d,    0.500000000000000000,    0.500000000000000000,    0.500000000000000000
     H,   96g,    0.312500000000000000,    0.312500000000000000,    0.125000000000000000
```

`ikdsplit.toml`

```toml
space_group_number = 227

cell = 6.555923429801207014

[fill]
Ti = ['Ti']
Cr = ['Cr']
H = ['H', 'X']

[[regress.transformations]]  # origin choice 2 -> 1
basis_change = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
origin_shift = [-0.125, -0.125, -0.125]
[[regress.transformations]]  # conventional -> 2x2x2 primitive
basis_change = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
origin_shift = [0.0, 0.0, 0.0]

[sort]
reference = "FPOSCAR"  # atoms are sorted in this order
```

`FPOSCAR`

## Usage

Run

```bash
ikdsplit run -l 2  # 2 levels of maximal subgroups are checked
```

Check the number of obtained `SPOSCAR-*`

```bash
find . -name "SPOSCAR-*" | wc
```

## Output files

- `wycksplit.toml`: information of Wyckoff splitting
- `PPOSCAR-*`: primitive cell of the subgroup
- `CPOSCAR-*`: conventional cell of the subgroup
- `RPOSCAR-*`: atoms with the target supercell
- `SPOSCAR-*`: `RPOSCAR-*` sorted in the same order as `FPOSCAR`

## Release notes

### 0.1.0

- Initial release working only for `227` and up to level `2`
