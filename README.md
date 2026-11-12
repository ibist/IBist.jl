# IBist.jl

Native Julia teaching data and biostatistical inference for proportions
and rates. This is the Julia sister package to
[ibist-R](https://github.com/ibist/ibist-R) and
[ibist-py](https://github.com/ibist/ibist-py). The package/module name is
`IBist`; the GitHub repository is `ibist/IBist.jl` and the local workspace
directory is `IBist.jl`.

This initial release provides seven datasets and
`or_2x2_ci()` with Wald, Gart, Cornfield, and Baptista–Pike intervals,
including conditional mid-p modifications and one-sided bounds.
It requires Julia 1.10 or later. R and Python are not required at runtime.

## Local development

The package is not yet registered. Install it from GitHub in Julia:

```julia
using Pkg
Pkg.add(url = "https://github.com/ibist/IBist.jl")
```

For development, clone the repository and run from its directory:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

To use it in another Julia environment, run
`Pkg.develop(path="/path/to/IBist.jl")`.
Dependencies and compatibility bounds live in `Project.toml`.

```julia
using IBist

result = or_2x2_ci([7 27; 1 33]; method = ["wald", "cornfield", "bp"])
result.estimate
dataframe(result)
or_2x2_ci([7 27; 1 33]; method = "bp", midp = true, side = "lower")

dataset = load_dataset("rds")
dataframe(dataset)
dataset.metadata["columns"]
```

DataFrames is a regular dependency, so the tabular teaching examples work
with a standard installation. Dataset columns preserve categorical level
order and missing values; definitions and provenance stay in `.metadata`.
Both datasets and interval results implement Tables.jl.

## Scope and validation

The statistical specification follows R. Shared references supplement
independent formulas, combinatorial probabilities, published endpoints,
conditional coverage enumeration, and invariance/boundary checks.
Tests also verify every exported data value against the portable source,
all missing-value counts, and factor order. See
[interval methods](docs/src/intervals.md) and
[datasets](docs/src/datasets.md).

Other interval families, tests, power, GLM methods, and approximation
diagnostics are not yet implemented. Shared C++ algorithms such as Wang's
will use bindings to the authoritative core rather than a new implementation.
The authoritative cross-language development policy is
[`CROSS_LANGUAGE_DEVELOPMENT.md`][policy] in the R repository.

[policy]: https://github.com/ibist/ibist-R/blob/main/CROSS_LANGUAGE_DEVELOPMENT.md

## Documentation

Public functions have Julia help strings (`?or_2x2_ci`, `?load_dataset`).
Documenter builds the manual from those strings and executable examples:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

CI checks Julia 1.10 and the current stable Julia on Linux, macOS, and
Windows, and builds the documentation. Registration is initiated manually;
TagBot creates release tags after General accepts a version. Documentation
deployment is not automatic. See [release instructions](docs/src/releases.md).

## License

Copyright (C) 2026 Elizabeth Schifano and Jun Yan.
GNU General Public License, version 3 or later; see [LICENSE](LICENSE).
