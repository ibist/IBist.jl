using Documenter
using IBist

makedocs(
    sitename = "IBist.jl",
    modules = [IBist],
    remotes = nothing,
    format = Documenter.HTML(prettyurls = false, repolink = nothing,
                            edit_link = nothing),
    pages = ["Overview" => "index.md", "Odds-ratio intervals" => "intervals.md",
             "Datasets" => "datasets.md", "API" => "api.md",
             "Releases" => "releases.md"],
    checkdocs = :exports,
)
