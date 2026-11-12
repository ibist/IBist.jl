# Odds-ratio intervals

Pass a matrix `[a b; c d]`, with independent groups in rows and events in
the first column. The sample odds ratio is `a*d/(b*c)`.

| Method | Construction |
| --- | --- |
| `wald` | Log-Wald interval; all cells must be positive |
| `gart` | Log-Wald interval with the Haldane–Anscombe 0.5 cell adjustment |
| `cornfield` | Inversion of the two one-sided conditional tails |
| `bp` | Inversion of a probability-ordered two-sided conditional test |

Conditional methods fix both margins. `midp=true` subtracts half the observed
mass from Cornfield's inclusive tails. For BP, half the mass of the complete
group tied with the observed table is included. Exact limits can be
conservative; mid-p intervals do not guarantee nominal coverage.

```@example sided
using IBist
result = or_2x2_ci([7 27; 1 33]; method = ["cornfield", "bp"], midp = true,
                   side = "lower")
dataframe(result)
```

The default `side="two.sided"` uses half the error probability in each
Cornfield tail. `side="lower"` and `side="upper"` use the full error
probability and leave the other endpoint at `Inf` or 0. One-sided BP uses
the Cornfield tail. The R/Python compatibility alias `alternative` is
available, but conflicting arguments are errors.

Calculations are native Julia. Conditional weights use stable adjacent
ratios and log normalization. Numerical inversion uses bisection on the
log odds-ratio scale, including boundaries at jumps in BP's ordered p-value.
Tolerance therefore differs in scale from R's direct odds-ratio inversion.
Conditional enumeration can be expensive for wide supports. At low
confidence levels some BP mid-p requests cannot bracket an accepted point
and raise an error. R/Python approximation diagnostics are not yet ported.

Validation includes 60 shared two-sided R reference calls and eight
additional one-sided calls covering low confidence levels, mid-p, and
repeated methods. Comparisons allow `2e-6` relative or `2e-8` absolute
endpoint error because R solves directly in odds-ratio units. These checks
supplement direct Wald formulas, closed-form small-table limits, independent
BigInt-combination probabilities, published Fagerland endpoints, table
rearrangements, and exact coverage enumeration for selected fixed margins.
Cross-language agreement alone does not establish coverage.

See [`or_2x2_ci`](@ref) for references and complete boundary behavior.
