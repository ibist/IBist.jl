"""
    CIResult

Confidence-interval result with public fields `estimate`, `conf_int`
(one row per method, columns lower and upper), `methods`, `conf_level`,
`side`, `parameter`, `null_value`, and `midp`. Method order and repetitions
are retained. Stored numbers are never rounded for display.

Use `dataframe(result)` for a fresh method/lower/upper comparison table.
The result also implements Tables.jl. `side` is `"two.sided"`, `"lower"`,
or `"upper"`; the unused one-sided endpoint is the parameter boundary.
An interval-only result contains no hypothesis test or p-value.
"""
struct CIResult
    estimate::Float64
    conf_int::Matrix{Float64}
    methods::Vector{String}
    conf_level::Float64
    side::String
    parameter::String
    null_value::Float64
    midp::Bool
end

Tables.istable(::Type{CIResult}) = true
Tables.columnaccess(::Type{CIResult}) = true
Tables.columns(result::CIResult) = (
    method = copy(result.methods), lower = copy(result.conf_int[:, 1]),
    upper = copy(result.conf_int[:, 2]))

"""
    dataframe(x::Union{CIResult, Dataset})

Return an independent DataFrame of confidence limits or dataset columns.
Dataset variable definitions and provenance remain in `x.metadata`.
"""
dataframe(result::CIResult) = DataFrame(Tables.columns(result))

function Base.show(io::IO, result::CIResult)
    print(io, "CIResult(", result.parameter, ", estimate=",
          result.estimate, ", methods=", result.methods, ")")
end

function Base.show(io::IO, ::MIME"text/plain", result::CIResult)
    print(io, "Confidence intervals for ", result.parameter, "\n")
    @printf(io, "Sample estimate: %.4g\n", result.estimate)
    @printf(io, "%.4g%% confidence, %s\n", 100 * result.conf_level,
            result.side)
    width = max(6, maximum(length, result.methods))
    println(io, rpad("method", width), "       lower       upper")
    for (i, method) in enumerate(result.methods)
        print(io, rpad(method, width))
        @printf(io, "  %10.4g  %10.4g", result.conf_int[i, 1],
                result.conf_int[i, 2])
        i < length(result.methods) && println(io)
    end
end
