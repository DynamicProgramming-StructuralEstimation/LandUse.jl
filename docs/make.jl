push!(LOAD_PATH,"../src")
using Documenter, LandUse

makedocs(modules = [LandUse], 
         sitename = "LandUse.jl",
         pages = [
             "Home" => "index.md",
             "Single Region" => "model.md",
             "Country" => "country.md",
             "Interactive Views" => "interact.md",
             "Function Index" => "reference.md"
         ]
    )

deploydocs(repo = "github.com/floswald/LandUse.jl.git")
