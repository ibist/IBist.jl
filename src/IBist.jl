"""
Teaching data and biostatistical inference for proportions and rates.

Load prepared data with [`load_dataset`](@ref) and estimate a two-by-two
odds ratio with [`or_2x2_ci`](@ref). This native Julia implementation
requires neither R nor Python at runtime.
"""
module IBist

using CategoricalArrays: categorical, levels!
using CodecZlib: GzipDecompressorStream
using DataFrames: DataFrame
using Distributions: Normal, cquantile
using Printf: @printf
import JSON3
import Tables

export CIResult, Dataset, available_datasets, dataframe, load_dataset,
       or_2x2_ci

include("results.jl")
include("datasets.jl")
include("odds_ratio.jl")

end
