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

cp(joinpath(@__DIR__,"..","LandUseR","docs"), joinpath(@__DIR__,"build","Rdocs"),force = true)


deploydocs(repo = "github.com/floswald/LandUse.jl.git")
