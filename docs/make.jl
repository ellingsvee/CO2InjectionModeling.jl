using Documenter, CO2BatchFill

push!(LOAD_PATH, "../src/")

makedocs(
    modules=[CO2BatchFill],
    sitename="CO2BatchFill",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        edit_link="main"
    ),
    warnonly=[:missing_docs],
    pages=[
        "Utilities" => "utils.md",
    ]
)

deploydocs(;
    repo="github.com/ellingsvee/CO2BatchFill.jl",
    devbranch="main",
)
