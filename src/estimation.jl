
"""
returns a dict of empirical targets for the model
"""
function targets(p::Param)

    # vector-value moments: time varying stuff
    m = OrderedDict()
    m[:rural_empl] = DataFrame(moment = ["rural_emp_" .* string(i) for i in p.moments.year], 
               data = p.moments.Employment_rural,
               model = zeros(length(p.moments.year)))

    # single value stuff: ratios of change etc
    # from average city over time
    m[:avg_density_fall] = DataFrame(moment = "avg_density_fall", data = 9.0, model = 0.0)
    m[:city_area] = DataFrame(moment = "city_area",data = 0.18, model = 0.0)
    m[:max_mode_increase] =  DataFrame(moment = "max_mode_increase", data = 7.0, model = 0.0)

    # average city spatial moments in 2020
    m[:density_gradient_2020] =  DataFrame(moment = "density_gradient_2020", data = 6.0 , model = 0.0)   # 1st tenth is 6 times denser than last 10-th

    # cross city moments
    # slope pop vs density
    m[:pop_vs_density_1876] = DataFrame(moment = "pop_vs_density_1876", data = 0.297, model = 0.0)
    m[:pop_vs_density_2015] = DataFrame(moment = "pop_vs_density_2015", data = 0.005, model = 0.0)

    return m
    
end

function dicts2df(d::OrderedDict)
    df = copy(d[:rural_empl])
    for (k,v) in d
        if k != :rural_empl
            append!(df,v)
        end
    end
    df
end

function x2dict(x)
    di = Dict(:cbar => x[1],
    :sbar => x[2],
    :ηl => x[3],
    :ηw => x[4],
    :ηm => x[5],
    :ϵs => x[6],
    :ϵsmax => x[6],
    :cτ => x[7],
    :ϕ1x => 0.15)
    di
end

"""
moment objective function for an optimizer
"""
function objective(x)
    # unpack X
    di = x2dict(x)

    p = Param(par = di, use_estimatedθ = false)
    try
        # find index of year 2020
        i2020 = argmin( abs.(p.moments.year .- 2020) )
        i2015 = argmin( abs.(p.moments.year .- 2015) )

        # run single city
        x1,M1,p1 = LandUse.run(p, estimateθ = false)
        d1 = dataframe(M1,p1)
        # multi country	
        xk,C,pk = runk(par = merge(di,Dict(:K => 2, :kshare => [0.5,0.5], :factors => [1.0,1.05])))

        # get data moments
        ta = targets(p)
        ta[:rural_empl].model = d1.rural_emp_model
        ta[:rural_empl].weights = ones(nrow(d1)) * 0.8

        # m = 0.0
        # m += sum(ta[:rural_empl].weights .* (ta[:rural_empl].data .- ta[:rural_empl].model).^2)

        ta[:avg_density_fall][!,:model] .= d1.citydensity[1] / d1.citydensity[i2020]
        ta[:avg_density_fall][!,:weights] .= 2.0

        # m += .(ta[:avg_density_fall].data .- ta[:avg_density_fall].model).^2

        ta[:city_area][!,:model] .= d1.cityarea[i2015]
        ta[:city_area][!,:weights] .= 1.0

        ta[:max_mode_increase][!,:model] .= maximum(d1.imode ./ d1.imode[1])
        ta[:max_mode_increase][!,:weights] .= 0.5

        ta[:density_gradient_2020][!,:model] .= M1[i2020].iDensity_q10 / M1[i2020].iDensity_q90
        ta[:density_gradient_2020][!,:weights] .= 1.0

        ta[:pop_vs_density_1876][!,:model]   .= (C[1].R[2].cityarea - C[1].R[1].cityarea) / (C[1].R[2].Lu - C[1].R[1].Lu)
        ta[:pop_vs_density_1876][!,:weights] .= 1.0
        ta[:pop_vs_density_2015][!,:model]   .= (C[end].R[2].cityarea - C[end].R[1].cityarea) / (C[end].R[2].Lu - C[end].R[1].Lu)
        ta[:pop_vs_density_2015][!,:weights] .= 1.0

        da = ta[:rural_empl]
        append!(da, ta[:avg_density_fall])
        append!(da, ta[:city_area])
        append!(da, ta[:max_mode_increase])
        append!(da, ta[:density_gradient_2020])
        append!(da, ta[:pop_vs_density_1876])
        append!(da, ta[:pop_vs_density_2015])

        sum(da.weights .* (da.data .- da.model).^2)
    catch 
        # @info "error at $(di)"
        return 999.9
    end
end

function runestim()
    x = bboptimize(
        objective, 
        SearchRange = [(0.7, 0.9),(0.1, 0.26), (0.0,0.5), (0.0, 0.5), (0.9, 2.0), (4.3, 8.0), (3.5, 4.5)], MaxFuncEvals = 500)

    LandUse.plot_i0k(x2dict(x))
end

