push!(LOAD_PATH,"../src")
using Documenter, LandUse

makedocs(modules = [LandUse], sitename = "LandUse.jl")

cp(joinpath(@__DIR__,"..","LandUseR","docs"), joinpath(@__DIR__,"build","Rdocs"),force = true)


<<<<<<< HEAD
deploydocs(repo = "github.com/floswald/LandUse.jl.git")
=======
# deploydocs(repo = "github.com/floswald/LandUse.jl.git")
>>>>>>> c788617d3058b4d9378b0a27b34fe0ba5f51b591
