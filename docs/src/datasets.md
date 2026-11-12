# Teaching datasets

```@example data
using IBist
available_datasets()
```

The datasets are `asah`, `depress`, `fluorosis`, `midbp`, `mihdl`, `nrs`,
and `rds`. `load_dataset(name)` returns a Dataset containing a DataFrame
and metadata. Each load is independent; `dataframe(dataset)` makes a deep
copy. Integer and numeric missing values use Julia's `missing`.
Categorical columns preserve their exported levels and ordered status.
Numeric value labels remain metadata and are not automatically recoded.

Prepared data, variable definitions, units, source references, category order,
and provenance come from the same R exporter used for Python. From ibist-R:

```sh
Rscript --vanilla data-raw/export-julia.R
```

The exporter accepts an optional Julia repository path. It invokes the
existing portable data exporter with a different output directory and
generates the comparison fixtures from R. It does not alter the Python
exports, R prepared objects, or cleaning rules. Generated files in `data/`
and `test/reference/` must not be edited manually.

Schema version 1 uses an index, JSON metadata, and gzip-compressed JSON
column arrays. Metadata retain the R revision, prepared-file and script
hashes, export hash, and raw-input availability. Some original NRS, RDS,
and depression inputs are absent from the workspace, so their prepared R
objects are exported as they stand; no observations are reconstructed.
The ASAH rows reconstruct separate published scale tables, with no valid
participant linkage across scales. Installed users need no network, R,
or Python to load any dataset.

R verifies every value on export with zero round-trip tolerance. Julia
checks every loaded value against the portable columns, all missing counts,
dimensions, category levels and order, and independent loads. Provenance
indicates uncommitted preparation changes when generated before those changes
are committed; regenerate from the final R revision before release.
