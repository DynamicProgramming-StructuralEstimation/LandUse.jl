using Documenter, LandUse

makedocs(
    modules = [LandUse],
    checkdocs = :exports,
    sitename = "LandUse.jl",
    pages = Any["index.md"],
    repo = "https://gitlab.com/floswald/LandUse.jl/blob/{commit}{path}#{line}"
)
