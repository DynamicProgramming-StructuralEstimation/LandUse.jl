

function interact()
    # θrs = 0.32:0.02:5.5
    θrs = [0.32, 0.33, 0.34, 0.36, 0.38, 0.41, 0.48, 0.7, 1.35, 2.3, 3, 4.5, 5, 5.5]
    θus = 0.32:0.02:5.5
    ηs = [0.0, 0.05, 0.1]

    # mp = @manipulate for θr in slider(θrs, label = "θr"), θu in slider(θus, label = "θu")
    mp = @manipulate for θ in slider(θrs, label = "θ"), η in slider(ηs)
        LandUse.model(pars = Dict(:θr => θ, :θu => θ, :η => η))
    end
end
