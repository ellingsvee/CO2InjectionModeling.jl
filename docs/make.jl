using Documenter, Literate, CO2BatchFill

push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "../examples/")

# Prepare example scripts
function set_to_cairo(content)
    content = replace(content, "GLMakie" => "CairoMakie")
    content = replace(content, "`CairoMakie`" => "GLMakie")
    return content
end

Literate.markdown("examples/multi_layer_filling.jl", "docs/src/"; execute=false, preprocess=set_to_cairo)
Literate.markdown("examples/multi_layer_ensemble.jl", "docs/src/"; execute=false, preprocess=set_to_cairo)
Literate.markdown("examples/sleipner_analysis.jl", "docs/src/"; execute=false, preprocess=set_to_cairo)


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
            "Core Types" => "types.md",
            "Topography" => "topography.md",
            "Layer Analysis" => "layer_analysis.md",
            "Simulation" => "simulation.md",
            "Analysis" => "analysis.md",
            "Unit Conversion" => "unit_conversion.md",
            "Visualization" => "visualization.md",
            "Utilities" => "utils.md",
        ],
        "Examples" => [
            "Multi-layer Filling" => "multi_layer_filling.md",
            "Multi-layer Ensemble" => "multi_layer_ensemble.md",
            "Sleipner Analysis" => "sleipner_analysis.md",
        ],
        "Index" => "indexlist.md",
    ]
)

deploydocs(;
    repo="github.com/ellingsvee/CO2BatchFill.jl",
    devbranch="main",
)
