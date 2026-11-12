"""
    or_2x2_ci(x; conf_level=0.95, method="wald", midp=false,
              side=nothing, alternative=nothing, tol=1e-8)

Estimate the odds ratio in a 2 × 2 matrix `[a b; c d]` of nonnegative
integer counts. The reported estimate is always `a*d/(b*c)`, including
when Gart's adjusted interval is requested. Return a [`CIResult`](@ref).

# Arguments
- `conf_level`: confidence level strictly between zero and one.
- `method`: a method name or vector/tuple of names. Order and repetitions
  are retained. `"wald"` uses the log-Wald interval and requires positive
  cells. `"gart"` adds 0.5 to every cell for interval construction.
  `"cornfield"` inverts the two conditional tails. `"bp"` inverts the
  probability-ordered, two-sided conditional test.
- `midp`: Boolean. For Cornfield subtract half the observed-table mass
  from each inclusive tail. BP assigns half the probability to the entire
  group tied with the observed table. Wald and Gart are unchanged.
- `side`: `"two.sided"` (default), `"lower"`, or `"upper"`. One-sided
  limits use the full error probability, with the other endpoint at 0 or
  `Inf`. One-sided BP uses the same ordered tail as Cornfield.
- `alternative`: compatibility alias: `"two.sided"`, `"greater"` (lower
  bound), or `"less"` (upper bound). Conflicting arguments are errors.
- `tol`: positive, finite absolute tolerance on the log odds-ratio scale.

# Notes
Rows are the two independent groups; the first column contains events.
Conditional methods fix both margins and enumerate Fisher's noncentral
hypergeometric distribution. They require positive row and column totals;
wide supports can be expensive. Probability weights use adjacent-count
ratios and normalization on the log scale. BP can have jumps in its p-value;
numerical boundaries need not solve a continuous equality exactly. At very
low confidence levels a mid-p acceptance set can be empty, in which case
`ArgumentError` is raised when a bound cannot be bracketed.

Exact conditional intervals can be conservative. Mid-p does not guarantee
nominal coverage. Wald and Gart are large-sample methods; computability
alone does not establish adequate approximation. This initial Julia release
does not yet provide R/Python's separate `asymcheck` diagnostics.

Gart constructs its interval with the Haldane–Anscombe 0.5 adjustment;
the reported sample estimate remains unadjusted. A zero numerator gives
estimate 0; a zero denominator gives `Inf`; both zero give `NaN`.
The latter can occur with empty margins, for which Gart remains available.
Counts must not exceed `2^53-1`, matching Python's exact-integer input bound.
Invalid tables or options raise `ArgumentError`.

# References
- Anscombe, F. J. (1956). On estimating binomial response relations.
  Biometrika 43, 461–464. https://doi.org/10.1093/biomet/43.3-4.461
- Baptista, J., and Pike, M. C. (1977). Algorithm AS 115: Exact two-sided
  confidence limits for the odds ratio in a 2 × 2 table.
  Applied Statistics 26, 214–220.
- Fagerland, M. W. (2012). Exact and mid-p confidence intervals for the
  odds ratio. The Stata Journal 12, 505–514.
- Gart, J. J. (1966). Alternative analyses of contingency tables.
  Journal of the Royal Statistical Society, Series B 28, 164–179.
- Haldane, J. B. S. (1956). The estimation and significance of the logarithm
  of a ratio of frequencies. Annals of Human Genetics 20, 309–311.
  https://doi.org/10.1111/j.1469-1809.1955.tb01285.x

# Examples
```julia
using IBist
result = or_2x2_ci([7 27; 1 33]; method=["wald", "cornfield"])
result.estimate
dataframe(result)
or_2x2_ci([7 27; 1 33]; method="bp", midp=true, side="lower")
```
"""
function or_2x2_ci(x; conf_level = 0.95, method = "wald", midp = false,
                   side = nothing, alternative = nothing, tol = 1e-8)
    tab = _validate_table(x)
    level = _probability(conf_level, "conf_level")
    midp isa Bool || throw(ArgumentError("midp must be true or false"))
    if !(tol isa Real) || tol isa Bool || !isfinite(tol) || tol <= 0
        throw(ArgumentError("tol must be a positive finite number"))
    end
    tolerance = Float64(tol)
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("tol is outside the Float64 range"))
    selected_side = _ci_side(side, alternative)
    methods = method isa AbstractString ? [method] :
        method isa Union{AbstractVector, Tuple} ? collect(method) : []
    valid = ("wald", "gart", "cornfield", "bp")
    if isempty(methods) || !all(m -> m isa AbstractString && m ∈ valid,
                                methods)
        throw(ArgumentError("method must contain wald, gart, cornfield or bp"))
    end
    estimate = Float64(tab[1, 1]) * tab[2, 2] /
        (Float64(tab[1, 2]) * tab[2, 1])
    error = (1 - level) / (selected_side == "two.sided" ? 2 : 1)
    z = cquantile(Normal(), error)
    intervals = Matrix{Float64}(undef, length(methods), 2)
    labels = String[]
    conditional = nothing
    for (i, name) in enumerate(methods)
        if name ∈ ("wald", "gart")
            cells = Float64.(tab) .+ (name == "gart" ? 0.5 : 0)
            all(>(0), cells) ||
                throw(ArgumentError("Wald interval requires positive cells"))
            eta = log(cells[1, 1]) + log(cells[2, 2]) -
                log(cells[1, 2]) - log(cells[2, 1])
            se = sqrt(sum(inv, cells))
            lower, upper = exp(eta - z * se), exp(eta + z * se)
            push!(labels, name)
        else
            if isnothing(conditional)
                conditional = _ConditionalTable(tab)
            end
            lower, upper = _conditional_interval(conditional, name, level,
                                                 midp, selected_side, tolerance)
            push!(labels, "$name-$(midp ? "midp" : "exact")")
        end
        intervals[i, 1] = selected_side == "upper" ? 0 : lower
        intervals[i, 2] = selected_side == "lower" ? Inf : upper
    end
    return CIResult(estimate, intervals, labels, level, selected_side,
                    "odds ratio", 1.0, midp)
end

function _probability(value, name)
    if !(value isa Real) || value isa Bool || !isfinite(value) ||
        !(0 < value < 1)
        throw(ArgumentError("$name must be a number strictly between 0 and 1"))
    end
    converted = Float64(value)
    0 < converted < 1 || throw(ArgumentError("$name exceeds Float64 precision"))
    return converted
end

function _ci_side(side, alternative)
    sides = ("two.sided", "lower", "upper")
    if !isnothing(side) && !(side isa AbstractString && side ∈ sides)
        throw(ArgumentError("side must be two.sided, lower or upper"))
    end
    if !isnothing(alternative)
        choices = ("two.sided", "greater", "less")
        index = alternative isa AbstractString ?
            findfirst(==(alternative), choices) : nothing
        if isnothing(index)
            throw(ArgumentError(
                "alternative must be two.sided, greater or less"))
        end
        legacy = sides[index]
        if !isnothing(side) && side != legacy
            throw(ArgumentError("side and alternative specify conflicting " *
                                "bounds"))
        end
        side = legacy
    end
    return isnothing(side) ? "two.sided" : String(side)
end

function _validate_table(x)
    if !(x isa AbstractMatrix) || size(x) != (2, 2) ||
        !all(v -> v isa Real && !(v isa Bool), x)
        throw(ArgumentError("x must be a 2 x 2 numeric matrix"))
    end
    if !all(v -> isfinite(v) && v >= 0 && isinteger(v), x)
        throw(ArgumentError("x must contain non-negative integer counts"))
    end
    all(<=(2^53 - 1), x) ||
        throw(ArgumentError("counts must not exceed 2^53-1"))
    return Matrix{Int}(x)
end

struct _ConditionalTable
    support::UnitRange{Int}
    base::Vector{Float64}
    index::Int
    center::Float64
end

function _ConditionalTable(tab::Matrix{Int})
    row1, row2 = sum(tab[1, :]), sum(tab[2, :])
    col1, col2 = sum(tab[:, 1]), sum(tab[:, 2])
    min(row1, row2, col1, col2) > 0 ||
        throw(ArgumentError("Need positive row and column totals"))
    support = max(0, col1 - row2):min(row1, col1)
    base = zeros(length(support))
    for i in 1:(length(support) - 1)
        s = support[i]
        base[i + 1] = base[i] + log(Float64(row1 - s)) +
            log(Float64(col1 - s)) - log(Float64(s + 1)) -
            log(Float64(row2 - col1 + s + 1))
    end
    index = tab[1, 1] - first(support) + 1
    base .-= base[index]
    center = all(>(0), tab) ? log(Float64(tab[1, 1])) +
        log(Float64(tab[2, 2])) - log(Float64(tab[1, 2])) -
        log(Float64(tab[2, 1])) : 0.0
    return _ConditionalTable(support, base, index, center)
end

function _log_probabilities(table::_ConditionalTable, eta)
    observed = table.support[table.index]
    weights = table.base .+ (table.support .- observed) .* eta
    maximum_weight = maximum(weights)
    shifted = weights .- maximum_weight
    return shifted .- log(sum(exp, shifted))
end

function _conditional_pvalue(table, eta, method, midp, lower)
    log_prob = _log_probabilities(table, eta)
    probabilities = exp.(log_prob)
    observed = table.index
    if method == "cornfield"
        tail = lower ? sum(@view probabilities[observed:end]) :
            sum(@view probabilities[1:observed])
        return tail - (midp ? probabilities[observed] / 2 : 0)
    end
    difference = log_prob .- log_prob[observed]
    tied = abs.(difference) .<= 1e-12
    smaller = (difference .< 0) .& .!tied
    return sum(probabilities[smaller]) +
        (midp ? 0.5 : 1) * sum(probabilities[tied])
end

function _bisect_boundary(objective, a, b, tol)
    fa, fb = objective(a), objective(b)
    fa == 0 && return a
    fb == 0 && return b
    signbit(fa) != signbit(fb) ||
        throw(ArgumentError("Cannot bracket conditional interval"))
    while b - a > tol
        middle = a + (b - a) / 2
        (middle == a || middle == b) && break
        fm = objective(middle)
        fm == 0 && return middle
        if signbit(fm) == signbit(fa)
            a, fa = middle, fm
        else
            b = middle
        end
    end
    return a + (b - a) / 2
end

function _conditional_interval(table, method, level, midp, side, tol)
    method = side == "two.sided" ? method : "cornfield"
    divisor = method == "cornfield" && side == "two.sided" ? 2 : 1
    target = (1 - level) / divisor
    limits = Float64[]
    step = log(10.0)
    for lower in (true, false)
        if lower && (side == "upper" || table.index == 1)
            push!(limits, 0.0)
            continue
        elseif !lower && (side == "lower" || table.index == length(table.base))
            push!(limits, Inf)
            continue
        end
        objective(eta) = _conditional_pvalue(table, eta, method, midp, lower) -
            target
        center = table.center
        if method == "cornfield"
            direction = lower ? 1 : -1
            a = b = center
            while direction * objective(a) > 0 && a > -750
                a -= step
            end
            while direction * objective(b) < 0 && b < 750
                b += step
            end
            root = _bisect_boundary(objective, a, b, tol)
        else
            if objective(center) < 0
                left = table.index > 1 ? table.base[table.index - 1] : -Inf
                right = table.index < length(table.base) ?
                    -table.base[table.index + 1] : Inf
                candidates = if isfinite(left) && isfinite(right)
                    [left + f * (right - left) for f in (0.5, 1e-6, 1 - 1e-6)]
                else
                    [isfinite(left) ? left + 1 : right - 1]
                end
                center = candidates[argmax(objective.(candidates))]
                objective(center) >= 0 ||
                    throw(ArgumentError("Cannot bracket BP at this level"))
            end
            direction = lower ? -1 : 1
            outside = center + direction * step
            while objective(outside) > 0 && abs(outside) < 750
                outside += direction * step
            end
            if objective(outside) > 0
                push!(limits, lower ? 0.0 : Inf)
                continue
            end
            a, b = minmax(center, outside)
            root = _bisect_boundary(objective, a, b, tol)
        end
        push!(limits, exp(root))
    end
    return limits[1], limits[2]
end
