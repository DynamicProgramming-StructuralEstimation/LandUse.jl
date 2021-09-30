# push!(LOAD_PATH,joinpath(@__DIR__,"..","src"))
using Documenter, LandUse, NNlib, Flux

makedocs(modules = [LandUse], 
         sitename = "LandUse.jl",
         pages = [
             "Home" => "index.md",
             "Single Region" => "model.md",
             "Country" => "country.md",
             "Interactive Views" => "interact.md",
             "Function Index" => "reference.md"
         ],
         format = Documenter.HTML(
         prettyurls = get(ENV, "CI", nothing) == "true")
    )

deploydocs(repo = "github.com/floswald/LandUse.jl.git")
