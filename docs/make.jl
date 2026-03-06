using Documenter, CO2BatchFill

push!(LOAD_PATH, "../src/")

makedocs(
    modules=[CO2BatchFill],
    sitename="CO2BatchFill",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        edit_link="main"
    ),
    warnonly=[:missing_docs, :docs_block, :cross_references],
    draft=false,
    pages=[
        "Introduction" => "index.md",
        "User Guide" => [
            "Core Types"         => "types.md",
            "Topography"         => "topography.md",
            "Layer Analysis"     => "layer_analysis.md",
            "Simulation"         => "simulation.md",
            "Analysis"           => "analysis.md",
            "Unit Conversion"    => "unit_conversion.md",
            "Visualization"      => "visualization.md",
            "Utilities"          => "utils.md",
        ],
        "Index" => "indexlist.md",
    ]
)

deploydocs(;
    repo="github.com/ellingsvee/CO2BatchFill.jl",
    devbranch="main",
)
