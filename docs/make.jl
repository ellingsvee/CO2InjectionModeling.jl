using CO2InjectionModeling
using Documenter

DocMeta.setdocmeta!(CO2InjectionModeling, :DocTestSetup, :(using CO2InjectionModeling); recursive=true)

makedocs(;
    modules=[CO2InjectionModeling],
    authors="Elling Svee",
    sitename="CO2InjectionModeling.jl",
    format=Documenter.HTML(;
        canonical="https://ellingsvee.github.io/CO2InjectionModeling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ellingsvee/CO2InjectionModeling.jl",
    devbranch="main",
)
