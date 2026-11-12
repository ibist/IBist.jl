@testset "Prepared data and metadata" begin
    expected_rows = Dict("asah" => 309, "depress" => 1730, "fluorosis" => 24,
                         "midbp" => 16, "mihdl" => 24, "nrs" => 91795,
                         "rds" => 16)
    @test Set(available_datasets()) == Set(keys(expected_rows))
    for name in available_datasets()
        dataset = load_dataset(name)
        @test length(dataset) == expected_rows[name]
        @test dataset.metadata["name"] == name
        @test size(dataset.data, 2) == dataset.metadata["ncols"]
        @test Tables.istable(typeof(dataset))
        frame = dataframe(dataset)
        columns = dataset.metadata["columns"]
        @test names(frame) == [info["name"] for info in columns]
        @test isequal(DataFrame(dataset), frame)
        raw = open(joinpath(IBist.DATA_DIR, "$name.json.gz")) do file
            stream = IBist.GzipDecompressorStream(file)
            try
                JSON3.read(read(stream, String), Dict{String, Any})
            finally
                close(stream)
            end
        end
        for info in dataset.metadata["columns"]
            column = frame[!, info["name"]]
            @test count(ismissing, column) == info["missing_count"]
            @test length(column) == expected_rows[name]
            if info["type"] == "category"
                @test levels(column) == info["levels"]
                @test isordered(column) == info["ordered"]
                @test all(isequal(ismissing(x) ? nothing : String(x), y)
                          for (x, y) in zip(column, raw[info["name"]]))
            else
                @test all(isequal(ismissing(x) ? nothing : x, y)
                          for (x, y) in zip(column, raw[info["name"]]))
            end
        end
        @test haskey(dataset.metadata["provenance"], "r_revision")
    end
    # Generated IDs do not turn the ASAH scale tables into linked participants.
    asah = load_dataset("asah")
    @test count(==("WFNS"), asah.data.scale) == 103
    notes = asah.metadata["notes"]
    @test any(occursin("cannot be linked", note) for note in notes)
    fluorosis = load_dataset("fluorosis")
    description = fluorosis.metadata["columns"][1]["description"]
    @test occursin("otitis media", lowercase(description))
    weights = levels(load_dataset("rds").data.bwt)
    @test all(occursin("g", String(value)) for value in weights)
    first = load_dataset("nrs")
    fresh = load_dataset("nrs")
    copied = dataframe(first)
    copied.id[1] = -100
    first.metadata["name"] = "changed"
    @test fresh.data.id[1] != -100
    @test first.data.id[1] != -100
    @test fresh.metadata["name"] == "nrs"
    @test_throws ArgumentError load_dataset("../rds")
    @test_throws ArgumentError load_dataset(1)
end
