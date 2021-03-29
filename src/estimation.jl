function post_slack(job)
	txt = "payload={'text': '$job'}"
	# println(txt)
	# println(ENV["MIG_SLACK"])
	if haskey(ENV,"SLACK_HOOK")
        if LandUse.isflo
		    Base.run(`curl -X POST --data-urlencode $txt $(ENV["SLACK_HOOK"])`) 
        else
            Base.run(`curl -X POST --insecure --data-urlencode $txt $(ENV["SLACK_HOOK"])`) 
        end
		return nothing
	else
		error("you need a webhook into slack as environment variable SLACK_HOOK to post a message")
	end
end
function post_slack()
	haskey(ENV,"SLACK_HOOK") || error("you need a webhook into slack as environment variable SLACK_HOOK to post a message")
end

function post_file_slack(fname::String)
    if LandUse.isflo
        Base.run(`curl -F file=@$(fname) -F "initial_comment=Dashboard at best candidate" -F channels=CCW1NHS1K -H "Authorization: Bearer $(ENV["SLACK_FILES"])" https://slack.com/api/files.upload`)
    else
        Base.run(`curl --insecure -F file=@$(fname) -F "initial_comment=Dashboard at best candidate" -F channels=CCW1NHS1K -H "Authorization: Bearer $(ENV["SLACK_FILES"])" https://slack.com/api/files.upload`)
    end
end
function post_file_slack()
    haskey(ENV,"SLACK_HOOK") || error("you need a webhook into slack as environment variable SLACK_FILES to post files")
end



"""
returns a dict of empirical targets for the model
"""
function targets(p::Param)

    # vector-value moments: time varying stuff
    m = Dict()
    m[:rural_empl] = DataFrame(moment = ["rural_emp_" .* string(i) for i in p.moments.year], 
               data = copy(p.moments.Employment_rural),
               model = zeros(length(p.moments.year)))

    # single value stuff: ratios of change etc
    # from average city over time
    m[:avg_density_fall] = DataFrame(moment = "avg_density_fall", data = 7.9, model = 0.0)
    m[:city_area] = DataFrame(moment = "city_area",data = 0.18, model = 0.0)
    m[:max_mode_increase] =  DataFrame(moment = "max_mode_increase", data = 5.0, model = 0.0)

    # average city spatial moments in 2020
    # m[:density_9010_2020] =  DataFrame(moment = "density_gradient_2020", data = 6.0 , model = 0.0)   # 1st tenth is 6 times denser than last 10-th
    # exponential decay model. exponential coefficient:
    m[:density_decay_coef] =  DataFrame(moment = "density_decay_coef", data = -0.15 , model = 0.0)   # on 21 points
    m[:density_decay_MSE] =  DataFrame(moment = "density_decay_MSE", data = 0.0 , model = 0.0) 

    # housing spending Shares
    m[:housing_share_1900] = DataFrame(moment = "housing_share_1900", data = 0.23 , model = 0.0)
    m[:housing_share_2015] = DataFrame(moment = "housing_share_2015", data = 0.314 , model = 0.0)

    # 

    # cross city moments
    # slope pop vs density
    # m[:pop_vs_density_1876] = DataFrame(moment = "pop_vs_density_1876", data = 0.297, model = 0.0)
    # m[:pop_vs_density_2015] = DataFrame(moment = "pop_vs_density_2015", data = 0.005, model = 0.0)

    return m
    
end

function dicts2df(d::Dict)
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
    # :ηm => x[5],
    :cτ => x[5],
    :ϕ1x => 0.15)
    di
end

search_over() = OrderedDict(zip([:cbar,:sbar,:ηl , :ηw, :cτ], [(0.7, 0.9),(0.2, 0.26), (0.08,0.15), (0.45, 0.55), (2.0, 3.0)]))
# search_over() = OrderedDict(zip([:cbar], [(0.7, 0.88)]))
symrange(x,pad,n) = range(x-pad, stop = x+pad, length = n)





function p2x(p::Param)
    [ p.cbar,  p.sbar, p.ηl, p.ηw, p.cτ ]
end

"""
moment objective function for an optimizer
"""
function objective(x; moments = false, plot = false, save = false)
    # unpack X
    di = x2dict(x)

    p = Param(par = di, use_estimatedθ = false)
    try
        # find index of year 2020
        i1870 = argmin( abs.(p.moments.year .- 1870) )
        i1900 = argmin( abs.(p.moments.year .- 1900) )
        i2020 = argmin( abs.(p.moments.year .- 2020) )
        i2015 = argmin( abs.(p.moments.year .- 2015) )

        # run single city
        x1,M1,p1 = LandUse.run(p, estimateθ = false)
        d1 = dataframe(M1,p1)
        # multi country	
        # xk,C,pk = runk(par = merge(di,Dict(:K => 2, :kshare => [0.5,0.5], :factors => [1.0,1.05])))

        # get data moments
        ta = targets(p)
        ta[:rural_empl].model = copy(d1.rural_emp_model)
        ta[:rural_empl].weights = ones(nrow(d1)) * 0.01

        # m = 0.0
        # m += sum(ta[:rural_empl].weights .* (ta[:rural_empl].data .- ta[:rural_empl].model).^2)

        ta[:avg_density_fall][!,:model] .= d1.citydensity[i1870] / d1.citydensity[i2015]
        ta[:avg_density_fall][!,:weights] .= 0.0

        ta[:city_area][!,:model] .= d1.cityarea[i2015]
        ta[:city_area][!,:weights] .= 50000.0

        ta[:max_mode_increase][!,:model] .= maximum(d1.imode ./ d1.imode[1])
        ta[:max_mode_increase][!,:weights] .= 5.0

        # ta[:density_gradient_2020][!,:model] .= M1[i2020].iDensity_q10 / M1[i2020].iDensity_q90
        # ta[:density_gradient_2020][!,:weights] .= 1.0

        # exponential model
        ndensities = M1[i2020].iDensities ./ M1[i2020].iDensities[1]
        gradient,emod = expmodel(1:p.int_bins, ndensities)
        MSE = round(1000 * mse(emod),digits = 3)

        ta[:density_decay_coef][!,:model] .= gradient[2]
        ta[:density_decay_coef][!,:weights] .= 0.0
        ta[:density_decay_MSE][!,:model] .= MSE
        ta[:density_decay_MSE][!,:weights] .= 0.0
    
        # housing spending Shares
        ta[:housing_share_1900][!,:model] .= d1[i1900,:Ch] / d1[i1900,:C]
        ta[:housing_share_1900][!,:weights] .= 0.1
        ta[:housing_share_2015][!,:model] .= d1[i2015,:Ch] / d1[i2015,:C]
        ta[:housing_share_2015][!,:weights] .= 0.1
        
        # ta[:pop_vs_density_1876][!,:model]   .= (C[1].R[2].cityarea - C[1].R[1].cityarea) / (C[1].R[2].Lu - C[1].R[1].Lu)
        # ta[:pop_vs_density_1876][!,:weights] .= 10.0
        # ta[:pop_vs_density_2015][!,:model]   .= (C[end].R[2].cityarea - C[end].R[1].cityarea) / (C[end].R[2].Lu - C[end].R[1].Lu)
        # ta[:pop_vs_density_2015][!,:weights] .= 10.0

        da = ta[:rural_empl]
        append!(da, ta[:avg_density_fall])
        append!(da, ta[:city_area])
        append!(da, ta[:max_mode_increase])
        append!(da, ta[:density_decay_coef])
        append!(da, ta[:density_decay_MSE])
        append!(da, ta[:housing_share_1900])
        append!(da, ta[:housing_share_2015])
        # append!(da, ta[:density_gradient_2020])
        # append!(da, ta[:pop_vs_density_1876])
        # append!(da, ta[:pop_vs_density_2015])

        if moments
            vv = sum(da.weights .* (da.data .- da.model).^2)
            if plot
                po = both_plots(M1,p1,i2020, objvalue = vv)
                if save 
                    savefig(po, joinpath(@__DIR__,"..","out","bboptim_$(Dates.today())_plot.pdf")) 
                end
                return (vv , da, po)
            else
                return (vv , da)
            end
        else
            return sum(da.weights .* (da.data .- da.model).^2)
        end

        
    catch 
        # @info "error at $(di)"
        return Inf
    end
end

function runestim(;steps = 1000)
    # check slack
    post_slack()
    post_file_slack()

    # mm = :dxnes


    if steps > 500
        halfstep = floor(steps/2)
        optctrl = bbsetup(objective ; SearchRange = collect(values(search_over())),MaxSteps = halfstep)
        res100 = bboptimize(optctrl)
        best100  = best_candidate(res100)
        println("Best candidate after $halfstep steps: ", best100)
        println("saving")
        # Now serialize to a temp file:
        fp = joinpath(@__DIR__,"..","out","bboptim_$(Dates.today())_$halfstep.dat")
        fh = open(fp, "w")
        serialize(fh, (optctrl, res100))
        close(fh)

        # continue
        fh = open(fp, "r")
        optctrlb, res100b = deserialize(fh);
        close(fh)

        post_slack("[LandUse] successfully done save and reload after $halfstep steps")


        res2 = bboptimize(optctrlb; MaxSteps = steps - halfstep)
        best  = best_candidate(res2)


        # final save:
        rm(fp)
        fp = joinpath(@__DIR__,"..","out","bboptim_$(Dates.today()).dat")
        fh = open(fp, "w")
        serialize(fh, (optctrlb, res2))
        close(fh)
    else
        optctrl = bbsetup(objective ; SearchRange = collect(values(search_over())),MaxSteps = steps)
        res100 = bboptimize(optctrl)
        best  = best_candidate(res100)
        println("Best candidate after $steps steps: ", best)
        println("saving")
        # Now serialize to a temp file:
        fp = joinpath(@__DIR__,"..","out","bboptim_$(Dates.today()).dat")
        fh = open(fp, "w")
        serialize(fh, (optctrl, res100))
        close(fh)

    end
    
    try 
        x,m,pl = objective(best, moments = true, plot = true, save = true)
        txt = """
        [LandUse.jl] Estimation finished on $(gethostname()) after $steps steps:

        *Results:*
        ========

        *best candidate*: 
        ```
        $(x2dict(best))
        ```

        *best moments*:
        ```
        $m
        ```
        """

        println(txt)
        post_file_slack(joinpath(@__DIR__,"..","out","bboptim_$(Dates.today())_plot.pdf"))
        post_slack(txt)

        return res100
    catch
        println("final eval of objective failed")
        txt = """
        [LandUse.jl] Estimation finished on $(gethostname())

        *Results:*
        ========

        *best candidate*: 
        ```
        $(x2dict(best))
        ```
        """
    end

    
end

"Estimate exponential model on spatial density. "
function expmodel(xdata,ydata)
    mod2(x,par) = par[1] * exp.(par[2] .* x)
    out = curve_fit(mod2, xdata, ydata , [1.0, -0.3])
    (coef(out), out)
end

function comp_expmodels(b0,grad)
    xdata = range(0, stop=10, length=20)
    ydata = b0 .* exp.(grad .* xdata) + 0.01*randn(length(xdata))

    # model 1 normalizes data to first point
    m1(x,par) = exp.(par .* x)
    e1 = coef(curve_fit(m1, xdata, ydata ./ ydata[1], [-0.3]))

    # model 2 does nothing
    m2(x,par) = par[1] * exp.(par[2] .* x)
    e2 = coef(curve_fit(m2, xdata, ydata, [1.0,-0.3]))

    p1 = scatter(xdata, ydata ./ ydata[1], title = "normalized",leg = false)
    plot!(p1, x -> exp.(e1[1] * x), 0, 10, annotations = ([7],[0.7],["exp coef=$(round(e1[1],digits=1))"]))

    p2 = scatter(xdata, ydata , title = "true",leg = false)
    plot!(p2, x -> e2[1] * exp.(e2[2] * x), 0, 10, annotations = ([7],[0.7 * b0],["exp coef=$(round(e2[2],digits=1))"]))

    plot(p1,p2)
end