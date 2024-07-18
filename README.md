# ikdsplit

## Usage

`cell.dat`

```
    6.555923429801207014     0.000000000000000000     0.000000000000000000
    0.000000000000000000     6.555923429801207014     0.000000000000000000
    0.000000000000000000     0.000000000000000000     6.555923429801207014
```

`atoms_conventional.csv`

```csv
symbol,wyckoff,x,y,z
          Ti,          8a,    0.125000000000000000,    0.125000000000000000,    0.125000000000000000
          Cr,         16d,    0.500000000000000000,    0.500000000000000000,    0.500000000000000000
           H,         96g,    0.312500000000000000,    0.312500000000000000,    0.125000000000000000
```

`ikdsplit.toml`

```toml
space_group_number = 227

[fill]
always = ["Ti", "Cr"]
never = []
selected = ["H"]

[[regress.transformations]]  # origin choice 2 -> 1
"basis_change" = [[ 1,  0,  0], [ 0,  1,  0], [ 0,  0,  1]]
"origin_shift" = [-0.12500, -0.12500, -0.12500]
[[regress.transformations]]  # conventional -> 2x2x2 primitive
"basis_change" = [[ 0,  1,  1], [ 1,  0,  1], [ 1,  1,  0]]
"origin_shift" = [ 0.00000,  0.00000,  0.00000]

[sort]
reference = '/Users/ikeda/Documents/projects/2022/Laves_AFM/data/codes/C15p_sub/FPOSCAR'
```

`FPOSCAR`

## Release notes

### 0.1.0
