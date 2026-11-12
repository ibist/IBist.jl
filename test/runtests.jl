using Test
using IBist
using CategoricalArrays: isordered, levels
using DataFrames: DataFrame, names
using Distributions: Normal, cquantile
using JSON3
using Tables

include("test_odds_ratio.jl")
include("test_datasets.jl")

@testset "Documented examples" begin
    result = or_2x2_ci([7 27; 1 33]; method = ["wald", "cornfield"])
    @test result.estimate ≈ 8.555555555555555
    @test size(dataframe(result)) == (2, 3)
    @test load_dataset("rds").metadata["name"] == "rds"
    @test occursin("cornfield-exact", sprint(show, MIME"text/plain"(), result))
    @test occursin("Dataset", sprint(show, load_dataset("rds")))
end
