# PubChem Property Lookup for R

Vectorized R functions to retrieve chemical compound properties from the [PubChem PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest), with local RDS caching for fast repeated lookups.

## Features

- **Lookup by compound name** — `pcprop_from_name_vec()` resolves CAS, common names, IUPAC names, abbreviations, and synonyms via PubChem
- **Lookup by InChIKey** — `pcprop_from_inchikey_vec()`
- **Automatic caching** — results are saved to `.rds` files so repeated queries hit the local cache instead of the API
- **Tidyverse-friendly** — returns tibble-compatible output; works seamlessly inside `mutate()` and `unnest()`

## Properties returned

By default, both functions retrieve these properties (configurable via the `property` argument):

| Property | Description |
|---|---|
| `CID` | PubChem Compound ID |
| `XLogP` | Predicted octanol/water partition coefficient |
| `ConnectivitySMILES` / `CanonicalSMILES` | SMILES string |
| `MolecularFormula` | Molecular formula |
| `MolecularWeight` | Molecular weight |
| `ExactMass` | Exact monoisotopic mass |
| `TPSA` | Topological polar surface area |
| `HeavyAtomCount` | Number of non-hydrogen atoms |
| `Charge` | Formal charge |
| `InChI` | IUPAC InChI identifier |
| `InChIKey` | Hashed InChI key |
| `IUPACName` | IUPAC systematic name |

## Requirements

- R ≥ 4.1
- `tidyverse` (dplyr, purrr, stringr, tidyselect)
- `curl` available on the system PATH

## Usage

```r
source("ms_scripts.R")
library(tidyverse)

# Create a tibble of compound names
cmps <- c(
  "perfluorooctanoic acid",
  "pentadecafluorooctanoic acid",
  "PFNA"
) |>
  as_tibble() |>
  rename(compound = 1)

# Look up properties and unnest into columns
cmps_info <- cmps |>
  mutate(pcprops = pcprop_from_name_vec(compound)) |>
  unnest(pcprops)

cmps_info
```

### Lookup by InChIKey

```r
inchikeys <- c("SNGREZUHAYWORS-UHFFFAOYSA-N")

result <- pcprop_from_inchikey_vec(inchikeys)
```

### Request a single property

When a single property is requested, a character vector is returned instead of a data frame:

```r
cmps |>
  mutate(formula = pcprop_from_name_vec(compound, property = "MolecularFormula"))
```

## Caching

Results are cached in RDS files in the working directory:

- `pubchem_all_prop_cache.rds` — name-based lookups
- `pubchem_inchikey_prop_cache.rds` — InChIKey-based lookups

Delete these files to force fresh lookups from PubChem.

## Citation

If you use this software in your work, please cite it. A DOI will be issued shortly via [Zenodo](https://zenodo.org/) — pending Zenodo availability. This section will be updated with the full citation once the DOI is live.

## License

MIT
