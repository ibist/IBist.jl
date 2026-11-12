# Independent combinatorial probabilities avoid the implementation's recurrence.
function integer_probabilities(tab, theta)
    row1, row2 = sum(tab[1, :]), sum(tab[2, :])
    col1 = sum(tab[:, 1])
    support = max(0, col1 - row2):min(row1, col1)
    weights = [binomial(big(row1), k) * binomial(big(row2), col1 - k) *
               BigFloat(theta)^k for k in support]
    return support, Float64.(weights ./ sum(weights))
end

@testset "Wald, Gart and interval results" begin
    tab = [7 27; 1 33]
    for side in ("two.sided", "lower", "upper")
        result = or_2x2_ci(tab; method = ["gart", "wald", "gart"], side)
        @test result.estimate ≈ 7 * 33 / 27
        @test result.methods == ["gart", "wald", "gart"]
        @test result.side == side
        for (i, adjustment) in enumerate((0.5, 0.0, 0.5))
            cells = tab .+ adjustment
            eta = log(cells[1, 1] * cells[2, 2] /
                      (cells[1, 2] * cells[2, 1]))
            z = cquantile(Normal(), 0.05 / (side == "two.sided" ? 2 : 1))
            limits = exp.(eta .+ [-1, 1] .* z .* sqrt(sum(inv, cells)))
            side == "lower" && (limits[2] = Inf)
            side == "upper" && (limits[1] = 0)
            @test result.conf_int[i, :] ≈ limits
        end
    end
    result = or_2x2_ci(tab; method = ["wald", "gart"])
    @test Tables.istable(typeof(result))
    @test DataFrame(result).method == result.methods
    exported = dataframe(result)
    exported.lower[1] = -99
    @test result.conf_int[1, 1] > 0
    @test or_2x2_ci(tab; side = "lower").conf_int ==
        or_2x2_ci(tab; alternative = "greater").conf_int
    @test or_2x2_ci(tab; side = "upper").conf_int ==
        or_2x2_ci(tab; alternative = "less").conf_int
end

@testset "Published BP example and closed-form Cornfield limits" begin
    # Fagerland (2012), Stata Journal 12, 505–514, mercii example.
    exact = or_2x2_ci([7 27; 1 33]; method = "bp")
    midp = or_2x2_ci([7 27; 1 33]; method = "bp", midp = true)
    @test isapprox(exact.conf_int[1, :], [1.000, 195.495]; atol = 0.004)
    @test isapprox(midp.conf_int[1, :], [1.328, 98.838]; atol = 0.003)
    for midp in (false, true)
        # Weights for [0 2; 2 0] are 1, 4θ, θ²; solve its upper tail.
        target = midp ? 0.05 : 0.025
        expected = -2 + sqrt(3 + 1 / target)
        ci = or_2x2_ci([0 2; 2 0]; method = "cornfield", midp).conf_int
        @test ci[1, 1] == 0
        @test isapprox(ci[1, 2], expected; rtol = 1e-8)
        one = or_2x2_ci([0 2; 2 0]; method = "cornfield", midp,
                        side = "upper").conf_int
        target = midp ? 0.1 : 0.05
        @test isapprox(one[1, 2], -2 + sqrt(3 + 1 / target); rtol = 1e-8)
    end
end

@testset "Conditional probabilities and exact coverage" begin
    for tab in ([7 27; 1 33], [0 2; 2 0], [1 10^12; 2 10^12])
        model = IBist._ConditionalTable(tab)
        for theta in (0.001, 0.1, 1.0, 10.0, 1000.0)
            _, expected = integer_probabilities(tab, theta)
            actual = exp.(IBist._log_probabilities(model, log(theta)))
            @test isapprox(actual, expected; rtol = 1e-12, atol = 1e-15)
        end
    end
    # Ties in BP include the complete equal-probability group, also for mid-p.
    model = IBist._ConditionalTable([0 2; 2 0])
    @test IBist._conditional_pvalue(model, 0.0, "bp", false, true) ≈ 1 / 3
    @test IBist._conditional_pvalue(model, 0.0, "bp", true, true) ≈ 1 / 6
    for method in ("cornfield", "bp")
        limits = [or_2x2_ci([k 4-k; 4-k 1+k]; method).conf_int[1, :]
                  for k in 0:4]
        for theta in (0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0)
            _, probability = integer_probabilities([0 4; 4 1], theta)
            covered = [lo <= theta * (1 + 1e-7) &&
                       hi >= theta * (1 - 1e-7) for (lo, hi) in limits]
            @test sum(probability[covered]) >= 0.95 - 1e-12
        end
    end
    # Independently check inclusive tail equations for all small tables.
    for a in 0:2, b in 0:2, c in 0:2, d in 0:2
        tab = [a b; c d]
        min(a + b, c + d, a + c, b + d) == 0 && continue
        result = or_2x2_ci(tab; method = "cornfield")
        for (i, theta) in enumerate(result.conf_int[1, :])
            (theta == 0 || !isfinite(theta)) && continue
            support, probability = integer_probabilities(tab, theta)
            tail = i == 1 ? sum(probability[support .>= a]) :
                sum(probability[support .<= a])
            @test isapprox(tail, 0.025; rtol = 1e-8, atol = 1e-12)
        end
    end
end

@testset "Boundary values and invariance" begin
    for midp in (false, true), method in ("wald", "gart", "cornfield", "bp")
        tab = [7 27; 1 33]
        first = or_2x2_ci(tab; method, midp)
        transposed = or_2x2_ci(permutedims(tab); method, midp)
        reversed = or_2x2_ci(tab[2:-1:1, :]; method, midp)
        @test first.conf_int ≈ transposed.conf_int
        @test isapprox(first.conf_int, 1 ./ reversed.conf_int[:, 2:-1:1];
                       rtol = 1e-7)
    end
    @test or_2x2_ci([0 2; 3 4]; method = "bp").estimate == 0
    @test or_2x2_ci([1 0; 3 4]; method = "bp").estimate == Inf
    @test isnan(or_2x2_ci([0 0; 3 4]; method = "gart").estimate)
    @test all(isfinite, or_2x2_ci(zeros(Int, 2, 2); method = "gart").conf_int)
    small = or_2x2_ci([7 27; 1 33]; method = "cornfield", midp = true,
                      conf_level = 0.01)
    @test small.conf_int[1, 1] < small.conf_int[1, 2]
    large = or_2x2_ci([100000 200000; 300000 400000])
    @test large.conf_int[1, 1] < large.estimate < large.conf_int[1, 2]
end

@testset "Input validation" begin
    for x in ([], [1, 2, 3, 4], ones(2, 3), trues(2, 2), ["1" "2"; "3" "4"])
        @test_throws ArgumentError or_2x2_ci(x)
    end
    for value in (-1.0, 0.5, NaN, Inf, 2.0^53)
        @test_throws ArgumentError or_2x2_ci([1 2; 3 value])
    end
    for level in (0, 1, -1, NaN, Inf, true, [0.95], "0.95")
        @test_throws ArgumentError or_2x2_ci([1 2; 3 4]; conf_level = level)
    end
    invalid_methods = ("adjusted", "unknown", [], ["wald", 2], nothing)
    for method in invalid_methods
        @test_throws ArgumentError or_2x2_ci([1 2; 3 4]; method)
    end
    for tol in (0, -1, Inf, NaN, true, "0.01")
        @test_throws ArgumentError or_2x2_ci([1 2; 3 4]; tol)
    end
    for midp in (nothing, 0, 1, "yes", [true])
        @test_throws ArgumentError or_2x2_ci([1 2; 3 4]; midp)
    end
    for side in (true, "less", :lower, ["lower"])
        @test_throws ArgumentError or_2x2_ci([1 2; 3 4]; side)
    end
    @test_throws ArgumentError or_2x2_ci([1 2; 3 4]; alternative = "lower")
    @test_throws ArgumentError or_2x2_ci([1 2; 3 4]; side = "lower",
                                        alternative = "less")
    @test_throws ArgumentError or_2x2_ci([0 2; 3 4])
    for method in ("cornfield", "bp")
        @test_throws ArgumentError or_2x2_ci([0 0; 3 4]; method)
    end
end

reference_number(value) = value == "inf" ? Inf :
    value == "nan" ? NaN : Float64(value)
reference_table(rows) = permutedims(hcat((Int.(r) for r in rows)...))

@testset "R two-sided references" begin
    path = joinpath(@__DIR__, "reference", "or_2x2_r.json")
    reference = JSON3.read(read(path, String), Dict{String, Any})
    for case in reference["cases"]
        labels = [ci["method"] for ci in case["intervals"]]
        methods = [first(split(label, "-")) for label in labels]
        result = or_2x2_ci(reference_table(case["table"]);
                           method = methods, conf_level = case["conf_level"],
                           midp = case["midp"])
        @test result.methods == labels
        @test result.estimate ≈ reference_number(case["estimate"])
        for (i, ci) in enumerate(case["intervals"])
            expected = reference_number.([ci["lower"], ci["upper"]])
            @test isapprox(result.conf_int[i, :], expected;
                           rtol = 2e-6, atol = 2e-8)
        end
    end
end

@testset "R one-sided references" begin
    path = joinpath(@__DIR__, "reference", "or_sides_r.json")
    reference = JSON3.read(read(path, String), Dict{String, Any})
    for case in reference["cases"]
        result = or_2x2_ci(reference_table(case["table"]);
                           method = case["methods"], side = case["side"],
                           conf_level = case["conf_level"], midp = case["midp"])
        @test result.methods == case["labels"]
        @test result.side == case["side"]
        for (i, limits) in enumerate(case["conf_int"])
            @test isapprox(result.conf_int[i, :], reference_number.(limits);
                           rtol = 2e-6, atol = 2e-8)
        end
    end
end
