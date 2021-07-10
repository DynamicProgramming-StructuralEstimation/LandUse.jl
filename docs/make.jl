push!(LOAD_PATH,"../src")
using Documenter, LandUse

makedocs(modules = [LandUse], sitename = "LandUse.jl")

cp(joinpath(@__DIR__,"..","LandUseR","docs"), joinpath(@__DIR__,"build","Rdocs"),force = true)


deploydocs(repo = "github.com/floswald/LandUse.jl.git")