# IBist.jl

IBist provides Julia interfaces to the settled statistical specification
shared with the R and Python ibist packages. The first development version
contains teaching datasets and odds-ratio intervals; other families remain
future work. It has no R/Python runtime dependency.

```@example overview
using IBist
result = or_2x2_ci([7 27; 1 33]; method = ["wald", "cornfield", "bp"])
dataframe(result)
```

The estimate is the unadjusted sample odds ratio. Confidence limits retain
full precision internally, independent of their printed representation.
Method requests retain order and repetitions. This interval-only result
does not imply a hypothesis test.

```@example overview
dataset = load_dataset("rds")
dataframe(dataset)
```

Use `dataset.metadata` for definitions, categories, units and provenance.
DataFrames is installed with IBist. Both Dataset and CIResult implement the
Tables.jl interface for other compatible consumers.
