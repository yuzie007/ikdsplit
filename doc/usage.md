# Usage

## Input files

3 input files are necessary.

`atoms_conventional.csv`

```csv
symbol,wyckoff,x,y,z
    Ti,    8a,    0.125000000000000,    0.125000000000000,    0.125000000000000
    Cr,   16d,    0.500000000000000,    0.500000000000000,    0.500000000000000
     H,   96g,    0.312500000000000,    0.312500000000000,    0.125000000000000
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

## Procedure

1. Fetch subgoup information using `ikdsplit split`
(Run recursively with the `-r` option)

```bash
ikdsplit split -r
```

2. Fill atoms in the primitive and the conventional cells
(Run up to the `l`-th maximal subgroups obtained from the original group)

```bash
ikdsplit fill -r -l 2
```

3. Regress atoms in the target supercell when possible

```bash
ikdsplit regress -r -l 2
```

4. Sort regressed atoms in the order specified in `FPOSCAR`

```bash
ikdsplit sort -r -l 2
```

## Tips

Check the number of obtained `SPOSCAR-*`

```bash
find . -name "SPOSCAR-*" | wc
```

## Output files

- `PPOSCAR-*`: primitive cell of the subgroup
- `CPOSCAR-*`: conventional cell of the subgroup
- `RPOSCAR-*`: atoms with the target supercell
- `SPOSCAR-*`: `RPOSCAR-*` sorted in the same order as `FPOSCAR`

### `atoms_{conventional,regressed}.csv`

This file contains Wyckoff positions for the subgroup.

### `info_{conventional,primitive,regressed}.csv`

This file contains the information on the actual space group number,
the numbers of atoms, and the elements filling Wyckoff positions for each
configuraiton.

```
configuration,space_group_number,Ti,Cr,H,4a1,8d1,16h1,32i1
0,227,4,8,48,Ti,Cr,H,H
1,141,4,8,16,Ti,Cr,H,X
2,141,4,8,32,Ti,Cr,X,H
3,227,4,8,0,Ti,Cr,X,X
```
