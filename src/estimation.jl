
"""
returns a dict of empirical targets for the model
"""
function targets(p::Param)

    # vector-value moments: time varying stuff
    m = Dict()
    m[:rural_empl] = p.moments.Employment_rural

    # single value stuff: ratios of change etc
    # from average city over time
    m[:avg_density_fall] = 9.0
    m[:city_area] = 0.18
    m[:max_mode_increase] = 7.0

    # average city spatial moments in 2020
    m[:density_gradient_2020] = 6.0   # 1st tenth is 6 times denser than last 10-th

    m
end

"""
moment objective function
"""
function objective(x)
    # unpack X
    p = Param(par = Dict(:cbar => x[1],
                         :sbar => x[2],
                         :ηl => x[3],
                         :ηw => x[4],
                         :ηm => x[5],
                         :ηl => x[3]))

end