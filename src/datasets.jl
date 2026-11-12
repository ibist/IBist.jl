const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))

"""
    Dataset

Prepared teaching data with public fields `data::DataFrame` and `metadata`.
Metadata preserve variable definitions, units, category levels and order,
numeric value labels, source references, and preparation provenance.
`length(dataset)` is the number of rows. Each load is independent.

The Dataset implements Tables.jl; use `dataframe(dataset)` for an independent
DataFrame. Missing values use Julia's `missing`; R factors become categorical
columns with their exported level order. Numeric value labels remain metadata
and do not silently convert numeric columns into categories.
"""
struct Dataset
    data::DataFrame
    metadata::Dict{String, Any}
end

Base.length(dataset::Dataset) = dataset.metadata["nrows"]
Tables.istable(::Type{Dataset}) = true
Tables.columnaccess(::Type{Dataset}) = true
Tables.columns(dataset::Dataset) = Tables.columns(dataset.data)
dataframe(dataset::Dataset) = deepcopy(dataset.data)

function Base.show(io::IO, dataset::Dataset)
    print(io, "Dataset(\"", dataset.metadata["name"], "\", ",
          length(dataset), " rows, ", size(dataset.data, 2), " columns)")
end

"""
    available_datasets()

List bundled teaching dataset names in export order. No network or R/Python
installation is required. Pass a name to [`load_dataset`](@ref).
"""
function available_datasets()
    index = JSON3.read(read(joinpath(DATA_DIR, "index.json"), String),
                       Dict{String, Any})
    index["schema_version"] == 1 || error("Unsupported dataset schema")
    return String.(index["datasets"])
end

"""
    load_dataset(name)

Read a bundled teaching Dataset. For example:

```julia
using IBist
dataset = load_dataset("rds")
frame = dataframe(dataset)
dataset.metadata["columns"]
```

Unknown names raise `ArgumentError`. All columns and metadata are newly read
on each load. Data are generated from R's existing preparation workflow;
installed users require neither R nor network access. Some original raw
inputs are unavailable in the workspace: their metadata record this, and
the export uses the prepared R objects without reconstructing observations.
See `docs/src/datasets.md` for regeneration and schema details.
"""
function load_dataset(name)
    if !(name isa AbstractString) || name ∉ available_datasets()
        throw(ArgumentError("Unknown dataset; choose from " *
                            join(available_datasets(), ", ")))
    end
    metadata = JSON3.read(read(joinpath(DATA_DIR, "$name.json"), String),
                          Dict{String, Any})
    metadata["schema_version"] == 1 || error("Unsupported dataset schema")
    values = open(joinpath(DATA_DIR, "$name.json.gz")) do file
        stream = GzipDecompressorStream(file)
        try
            JSON3.read(read(stream, String), Dict{String, Any})
        finally
            close(stream)
        end
    end
    frame = DataFrame()
    for info in metadata["columns"]
        column = values[info["name"]]
        type = info["type"]
        T = type == "integer" ? Int : type == "number" ? Float64 :
            type == "boolean" ? Bool : String
        converted = if info["missing_count"] > 0
            Union{Missing, T}[isnothing(x) ? missing : T(x) for x in column]
        else
            T.(column)
        end
        if type == "category"
            converted = categorical(converted; ordered = info["ordered"])
            levels!(converted, String.(info["levels"]))
        end
        frame[!, Symbol(info["name"])] = converted
    end
    size(frame) == (metadata["nrows"], metadata["ncols"]) ||
        error("Dataset dimensions do not match metadata")
    return Dataset(frame, metadata)
end
