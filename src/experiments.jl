
function k20_dfexport(;overwrite = false)
    if overwrite
        df = FileIO.load(joinpath(intables, "d1d220.jld2"))
        d0 = df["d0"]
        d1 = df["d1"]
        C0 = df["C0"]
        C1 = df["C1"]
        p1 = df["p1"]
        p0 = df["p0"]

        # get exponential coefs for each city
        coefs0 = [Array(hcat(expmodel(c)[1]...)')[:,2] for c in C0]
        coefs1 = [Array(hcat(expmodel(c)[1]...)')[:,2] for c in C1]

        dc0 = DataFrame(gridmake(1:20, p0.T) , [:region, :year])
        dc0.exp_coef = vcat(coefs0...)
        dc0[!,:d1] .= 0.0
        dc0[!,:d2] .= 0.0
        dc0[!,:a ] .= p0.a

        dc1 = DataFrame(gridmake(1:20, p0.T) , [:region, :year])
        dc1.exp_coef = vcat(coefs1...)
        dc1[!,:d1] .= p1.d1
        dc1[!,:d2] .= p1.d2
        dc1[!,:a ] .= p1.a

        d0 = leftjoin(d0,dc0,on = [:region, :year])
        d1 = leftjoin(d1,dc1,on = [:region, :year])

        # build data export
        yearmap = popdata_mapyears(p0)
        # popdata_mapyears(p0)
        d0 = leftjoin(d0, yearmap, on = :year => :modelyears)
        d1 = leftjoin(d1, yearmap, on = :year => :modelyears)
        d0 = leftjoin(d0, p0.citylist, on = [:region => :rank, :datayears => :year])
        d1 = leftjoin(d1, p0.citylist, on = [:region => :rank, :datayears => :year])

        CSV.write(joinpath(dbtables, "k20-baseline.csv"), d0)
        CSV.write(joinpath(dbtables, "k20-d1d2.csv"), d1)
        FileIO.save(joinpath(intables, "k20-d1d2.jld2"), Dict("d0" => d0, "d1" => d1))

    else
		df = FileIO.load(joinpath(intables, "k20-d1d2.jld2"))
		d0 = df["d0"]
		d1 = df["d1"]
    end
    d0,d1
end


"""
Produces output from 20-city model. In particular comparison with single city case, and extension with d1, d2.
"""
function k20output(k;d1_ = 0.05,d2_ = 2.0, a = 2.67, overwrite = false)

    K = k
    if overwrite
        x,C0,p0 = LandUse.runk(par = Dict(:K => K,:kshare => [1/K for i in 1:K], :factors => ones(k), :gs => zeros(k)))
        println("baseline done")
        x,C1,p1 = LandUse.runk(par = Dict(:K => K,:kshare => [1/K for i in 1:K], :factors => ones(k), :gs => zeros(k), :d1 => d1_, :d2 => d2_, :a => a))
        d0 = dataframe(C0)
        d1 = dataframe(C1)

		FileIO.save(joinpath(intables, "d1d2$k.jld2"), Dict("d0" => d0, "d1" => d1, "p0" => p0, "p1" => p1, "C1" => C1, "C0" => C0))	
       
    else
        df = FileIO.load(joinpath(intables, "d1d2$k.jld2"))
		d0 = df["d0"]
		d1 = df["d1"]
        C0 = df["C0"]
		C1 = df["C1"]
		p1 = df["p1"]
		p0 = df["p0"]

    end

    # run a single city
    x,M,p = runm()

    p0x = select(subset(d0, :year => x->x.== 2020), :year, :Lu, :citydensity => LandUse.firstnorm => :fn, :region) 
    pl0 = @df p0x bar(:fn,xticks = ([1,2,3],["Paris","Lyon","Marseille"]), ylab = "rel density", title = "baseline")
    annotate!(pl0, [(2,0.5, Plots.text("$(round(p0x[2,:fn],digits = 6))"))])

    p1x = select(subset(d1, :year => x->x.== 2020), :year, :Lu, :citydensity => LandUse.firstnorm => :fn, :region)
    pl1 = @df p1x bar(:fn,xticks = ([1,2,3],["Paris","Lyon","Marseille"]), ylab = "rel density", title = "d1 = $d1_, d2 = $d2_")
    annotate!(pl1, [(2,0.5, Plots.text("$(round(p1x[2,:fn],digits = 6))"))])

    dd0 = select(subset(d0, :year => x->x.== 2020), :year, :Lu, :cityarea, :citydensity, :region)
    dd1 = select(subset(d1, :year => x->x.== 2020), :year, :Lu, :cityarea, :citydensity, :region)
    xx0 = lm(@formula( log(cityarea) ~ log(Lu) ), dd0)
    xx1 = lm(@formula( log(cityarea) ~ log(Lu) ), dd1)
    # return (d0, d1)
    # (xx0, xx1, d0, d1)
    b0 = bar([coef(xx0)[2]],ylims = (0,1), title = "baseline",annotations = (1.0, 0.8, Plots.text("coef = $(round(coef(xx0)[2],digits = 6))")))
    b1 = bar([coef(xx1)[2]],ylims = (0,1), title = "d1 = $d1_, d2 = $d2_",annotations = (1.0, 0.8, Plots.text("coef = $(round(coef(xx1)[2],digits = 6))")))

    s0 = @df dd0 scatter(:Lu, :cityarea, title = "baseline")
    s1 = @df dd1 scatter(:Lu, :cityarea, title = "d1 = $d1_, d2 = $d2_, a=$a")

    ts0 = ts_plots([C0[i].R[1] for i in 1:length(p0.T)], p0)
    ts1 = ts_plots([C1[i].R[1] for i in 1:length(p0.T)], p0)

    ts20 = ts_plots([C0[i].R[2] for i in 1:length(p0.T)], p0)
    ts21 = ts_plots([C1[i].R[2] for i in 1:length(p0.T)], p0)

    avg0 = select(d0, :year, :region, :citydensity, :cityarea, :d0)
    avg1 = select(d1, :year, :region, :citydensity, :cityarea, :d0)

    a0 = @df avg0 plot(:year, :citydensity, group = :region, title = "baseline av density" ,ylims = (0,maximum(avg0.citydensity)), leg = false)
    a1 = @df avg1 plot(:year, :citydensity, group = :region, title = "d1 = $d1_, d2 = $d2_, a=$a",ylims = (0,maximum(avg0.citydensity)), leg = false)

    phi0 = @df avg0 plot(:year, :cityarea, group = :region, title = "baseline cityarea", ylims = (0,maximum(avg1.cityarea)), leg = false)
    phi1 = @df avg1 plot(:year, :cityarea, group = :region, title = "d1 = $d1_, d2 = $d2_, a=$a", ylims = (0,maximum(avg1.cityarea)), leg = false)

    pout = plot(b0,b1, 
        #  plot(ts0[:n_densities],title = "baseline, k=1"), 
        #  plot(ts1[:n_densities], title = "d1 = $d1_, d2 = $d2_, k=1"),
        #  plot(ts20[:n_densities],title = "baseline, k=2"),
        #  plot(ts21[:n_densities],title = "d1 = $d1_, d2 = $d2_, k=2"),
         a0,a1,
         phi0,phi1, 
        #  layout = (5,2), size = (800,1100))
         layout = (3,2), size = (800,800))

    # compare single city to top 5 after 1880
    k5 = subset(d0, :region => x -> x .< 6, :year => x -> x .> 1870)

    # get single city plots
    tsp = ts_plots(M,p)

    # area 
    area1 = tsp[:phi]
    for ik in 1:5
        tk5 = subset(k5, :region => x -> x .== ik)
        plot!(area1, tk5.year, tk5.rel_cityarea, color = :grey, label = tk5)
    end

    area1

    # Dict(:C0 => C0,:p0=>p0, :C1 => C1, :p1=>p1, :d0=>d0 , :d1=>d1 ,:plot=>pout)
end




"""
run model with flat epsilon

https://github.com/floswald/LandUse.jl/issues/59
"""
function issue59(;save = false)
    p1 = Param()
    p2 = Param(par = Dict(:ϵflat => true, :ϵr => 4.0))

    x1,m1,p1 = run(p1)
    x2,m2,p2 = run(p2)

    d1 = dataframe(m1,p1); d2 = dataframe(m2,p2)
    d1[!,:type] .= "baseline"
    d2[!,:type] .= "ϵ(l) = 4"
    d = vcat(d1,d2)

    pl = Dict()
    pl[:density] = @df d plot(:year, :citydensity, group = :type, title = "Average Density")
    pl[:size] = @df d plot(:year, :cityarea, group = :type, title = "Urban Area", leg=false)
    pl[:qbar] = @df d plot(:year, :qbar_real, group = :type, title = "Avg House Price", leg=false)

    pl[:ndens0] = @df d plot(:year, :d0_n, group = :type, title = "Central n density", leg=false)
    pl[:ndensr] = @df d plot(:year, :dr_n, group = :type, title = "Fringe n density", leg=false)
    pl[:ndensa] = @df d plot(:year, :avgd_n, group = :type, title = "Avg n density", leg=false)
    o = plot(pl[:density], pl[:size], pl[:qbar], pl[:ndens0], pl[:ndensa], pl[:ndensr], layout = (2,3), size = (800,400))

    cs1 = cs_plots(m1[19],p1,19)
    cs2 = cs_plots(m2[19],p2,19)
    pl[:gradients] = plot(cs1[:D],cs2[:D], title = ["baseline" "flat ϵ"])

    if save 
        savefig(o,joinpath(dbplots,"issue59.pdf"))
        savefig(pl[:gradients],joinpath(dbplots,"issue59-gradients.pdf"))
    end
    pl
end


"""
run model with agglomeration forces

https://github.com/floswald/LandUse.jl/issues/60
"""
function issue60(;save = false)
    p1 = Param()
    p2 = Param(par = Dict(:η => 0.1))

    x1,m1,p1 = run(p1)
    x2,m2,p2 = run(p2)

    d1 = dataframe(m1,p1); d2 = dataframe(m2,p2)
    d1[!,:type] .= "baseline"
    d2[!,:type] .= "η = $(p2.η)"
    d = vcat(d1,d2)

    pl = Dict()
    pl[:density] = @df d plot(:year, :citydensity, group = :type, title = "Average Density")
    pl[:size] = @df d plot(:year, :cityarea, group = :type, title = "Urban Area", leg=false)
    pl[:qbar] = @df d plot(:year, :qbar_real, group = :type, title = "Avg House Price", leg=false)

    pl[:ndens0] = @df d plot(:year, :d0_n, group = :type, title = "Central n density", leg=false)
    pl[:ndensr] = @df d plot(:year, :dr_n, group = :type, title = "Fringe n density", leg=false)
    pl[:ndensa] = @df d plot(:year, :avgd_n, group = :type, title = "Avg n density", leg=false)
    o = plot(pl[:density], pl[:size], pl[:qbar], pl[:ndens0], pl[:ndensa], pl[:ndensr], layout = (2,3), size = (800,400))


    if save 
        savefig(o,joinpath(dbplots,"issue60.pdf"))
        # for (k,v) in pl
        #     savefig(v,joinpath(dbplots,"issue60-$(string(k)).pdf"))
        # end 
    end
    o
end


"""
run model with congestion forces

https://github.com/floswald/LandUse.jl/issues/64
"""
function issue64(;save = false)
    p1 = Param()
    p2 = Param(par = Dict(:ηa => 0.1))

    x1,m1,p1 = run(p1)
    x2,m2,p2 = run(p2)

    d1 = dataframe(m1,p1); d2 = dataframe(m2,p2)
    d1[!,:type] .= "baseline"
    d2[!,:type] .= "ηa = $(p2.ηa)"
    transform!(d1, :imode => (x -> x ./ x[1]) => :imode_n)
    transform!(d2, :imode => (x -> x ./ x[1]) => :imode_n)
    d = vcat(d1,d2)

    pl = Dict()
    pl[:density] = @df d plot(:year, :citydensity, group = :type, title = "Average Density")
    pl[:size] = @df d plot(:year, :cityarea, group = :type, title = "Urban Area", leg=false)
    pl[:mode] = @df d plot(:year, :imode_n, group = :type, title = "Mode increase", leg=false)

    pl[:ndens0] = @df d plot(:year, :d0_n, group = :type, title = "Central n density", leg=false)
    pl[:ndensr] = @df d plot(:year, :dr_n, group = :type, title = "Fringe n density", leg=false)
    pl[:ndensa] = @df d plot(:year, :avgd_n, group = :type, title = "Avg n density", leg=false)
    o = plot(pl[:density], pl[:size], pl[:mode], pl[:ndens0], pl[:ndensa], pl[:ndensr], layout = (2,3), size = (800,400))


    if save 
        savefig(o,joinpath(dbplots,"issue64.pdf"))
        # for (k,v) in pl
        #     savefig(v,joinpath(dbplots,"issue60-$(string(k)).pdf"))
        # end 
    end
    o
end


"""
run multi city with congestion

https://github.com/floswald/LandUse.jl/issues/66
"""
function issue66(;save = false)
    p1 = Param()
    p2 = Param(par = Dict(:ηa => 0.1))
    d1 = Dict(:K => 3, :kshare => [1/3,1/3,1/3], :factors => [1.0,1.01,1.1])
    d2 = Dict(:K => 3, :kshare => [1/3,1/3,1/3], :factors => [1.0,1.01,1.1], :ηa => 0.1)
    

    x1,m1,p1 = runk(par = d1)
    x2,m2,p2 = runk(par = d2)

    c1 = dashboard(m1, 19)
    c2 = dashboard(m2, 19)

    o = plot(c1, c2, size = (2200,700))


    if save 
        savefig(o,joinpath(dbplots,"issue64.pdf"))
        # for (k,v) in pl
        #     savefig(v,joinpath(dbplots,"issue60-$(string(k)).pdf"))
        # end 
    end
    o
end


"""
five city cross section

https://github.com/floswald/LandUse.jl/issues/67
"""
function issue67()
    x,C,p = k5()
    d = LandUse.dataframe(C)
    g = groupby(d, :region)

    dd = select(transform(g, :Lu => firstnorm => :Lun, 
                             :cityarea => firstnorm => :cityarean,
                             :citydensity => firstnorm => :citydensityn), 
                 :year, :Lu, :cityarea, :Lun, :cityarean,:citydensity,:citydensityn, :region)

    pl = Dict()

    # ratio of smallest to largest city's population over time
    gg = groupby(select(filter(x -> x.region .∈ Ref([1,5]), d), :region, :year, :Lu), :year)
    dg = combine(gg, :Lu => (x -> maximum(x) / minimum(x)) => :rel_Lu)

    # Relative Population and Area in final period
    # ============================================

    data = CSV.File(joinpath(LandUse.dboutdata,"top5poparea.csv")) |> DataFrame
    data[!,:region] = 1:5
    data[!,:group] .= "data"
    sort!(data, :region)
    
    m = select(subset(d, :year => x -> x .== 2020), :region, :Lu => (x -> x ./ maximum(x)) => :relative_pop, :cityarea => (x -> x ./ maximum(x)) => :relative_area)
    m[!,:group] .= "model"

    dm = vcat(select(data,Not(:LIBGEO)),m)
    dm = leftjoin(dm, select(data,:region, :LIBGEO => :city), on = :region)

    pl[:relpop] = @df dm groupedbar(:relative_pop, group = :group, legend = :left, xticks = (:region,:city), title = "2020 Population relative to Paris")
    pl[:relarea] = @df dm groupedbar(:relative_area, group = :group, legend = :left, xticks = (:region,:city), title = "2020 Area relative to Paris")
    savefig(pl[:relpop], joinpath(dbplots,"five-city-relpop.pdf"))
    savefig(pl[:relarea], joinpath(dbplots,"five-city-relarea.pdf"))


    # commuting time
    # ==============
    ct = combine(g, :ictime => last, :imode => last)
    rct = transform(ct, :ictime_last => maxnorm => :reltime, :imode_last => maxnorm => :relspeed)
    rct = leftjoin(rct, select(data,:region,:LIBGEO => :city), on = :region)
    sct = stack(select(rct,:city, :region, :reltime, :relspeed), Not([:region, :city]))

    pl[:relspeedtime] = @df sct groupedbar(:value, group = :variable, legend = :topleft, xticks = (:region,:city), title = "Average Speed and Commute Time rel to Paris")
    savefig(pl[:relspeedtime], joinpath(dbplots,"five-city-relspeedtime.pdf"))





    # Density
    # =======

    # cross section: bigger cities are denser, in all periods
    pl[:cross] = @df dd plot(:year, :citydensity, group = :region, yaxis = :log10, ylab = "log density", title = "Bigger cities are always denser.")
    savefig(pl[:cross], joinpath(dbplots,"five-city-cross.pdf"))
    
    # over time, the fall in density is more pronounced in large cities than in smaller ones
    pl[:time] = @df dd plot(:year, :citydensityn, group = :region, ylab = "normalized density (1840 = 1)", title = "...but they become less dense faster!")
    savefig(pl[:time], joinpath(dbplots,"five-city-time.pdf"))

    # show relative density and population over time
    pk = @df subset(dd, :region => x -> x .== 1) plot(:citydensityn, :Lun, series_annotation = Plots.text.(:year, 8, :right),xlab = "relative density (1840 = 1)",
     ylab = "relative Population (1840 = 1)", title = "Density vs Urban Population", label = "1")
    for ik in 2:5
        ds = subset(dd, :region => x -> x .== ik)
        if ik < 5
            plot!(pk, ds.citydensityn, ds.Lun, label = "$ik")
        else
            years = string.(collect(p.T))
            years[Not([1,5,10,15,19])] .= ""
            plot!(pk, ds.citydensityn, ds.Lun, label = "$ik", series_annotation = Plots.text.(years, 8, :right))
        end
    end
    pl[:rel] = pk
    savefig(pl[:rel], joinpath(dbplots,"five-city-rel.pdf"))
    (pl,d,g,dg)
end



"""
plot ``\\phi_k`` vs ``\\L{u,k}`` for all regions ``k``.
relates to https://github.com/floswald/LandUse.jl/issues/9
"""
function issue9()
    cpar = Dict(:S => 1.0, :L => 1.0,
                :K => 4,
                :θprop => [1.0,0.99,0.98,0.97],
                :kshare => [0.25 for i in 1:4])
    sols,C,cpar,pp = LandUse.runk(cpar = cpar)

    K = C[1].K

    open(joinpath(dbtables,"phi-vs-Lu.txt"),"w") do io
        @printf(io,"Year   k    θu      Lu       ϕ\n")
        @printf(io,"----------------------------------\n")

        anim = Animation()
        for (jt,it) in enumerate(C[1].T)
            ϕs = [C[jt].R[ik].ϕ for ik in 1:K]; Lus = [C[jt].R[ik].Lu for ik in 1:K]
            for ik in 1:K @printf(io,"%d   %d    %1.2f    %1.3f    %1.3f\n",it,ik,pp[ik].θus[jt]*pp[ik].θprop,Lus[ik],ϕs[ik]) end

            pl = plot(ϕs, Lus, title = "$it",
                      leg = false,
                      l = (:black,2),
                      m = (:circle, 5, :red),
                      xaxis = (L"\phi", (0.001,0.045)),
                      yaxis = (L"L_u", (0.01,0.3)))
            frame(anim)
        end
        g = gif(anim, joinpath(dbplots,"phi-vs-Lu.gif"),fps=0.5)
    end

end

"""
plot ``\\phi_k`` vs ``L_{u,k}`` for all regions ``k`` and different eps slopes.
relates to https://github.com/floswald/LandUse.jl/issues/10
"""
function issue10()
    cpar = Dict(:S => 1.0, :L => 1.0,
                :K => 4,
                :θprop => [1.0,0.97,0.95,0.94],
                :kshare => [0.25 for i in 1:4])
    sols,C,cp,p0 = LandUse.runk(cpar = cpar)

    K = C[1].K
    ϵ1 = 5.0
    ϵ2 = 1.0
    ϵ3 = 0.0
    ϵs = [p0[1].ϵs,ϵ1,ϵ2,ϵ3]

    # version with lower eps slope
    sols,C2,cp,pp = LandUse.runk(cpar = cpar,
                                 par = Dict(ik => Dict(:ϵsmax => ϵs[2]) for ik in 1:K))

    # version with lower eps slope
    sols,C3,cp,pp = LandUse.runk(cpar = cpar,
                                 par = Dict(ik => Dict(:ϵsmax => ϵs[3]) for ik in 1:K))

    # version with lower eps slope
    sols,C4,cp,pp = LandUse.runk(cpar = cpar,
                              par = Dict(ik => Dict(:ϵsmax => ϵs[4]) for ik in 1:K))

    df = DataFrame()
    fmt = FormatExpr("{1:d}    {2:d}   {3:>4.1f}    {4:>1.2f}    {5:>1.2f}    {6:>1.2f}")
    # open(joinpath(dbtables,"phi-vs-Lu-eps.txt"),"w") do io
        # print(io,"Year    k   ϵs      θu      Lu      ϕ\n")
        # print(io,"----------------------------------------\n")

        anim = Animation()
        for (jt,it) in enumerate(C[1].T)
            # ϕs =  log.(1000 .* [C[jt].R[ik].ϕ for ik in 1:K]);   Lus = log.(1000 .* [C[jt].R[ik].Lu for ik in 1:K])
            # ϕs2 = log.(1000 .* [C2[jt].R[ik].ϕ for ik in 1:K]); Lus2 = log.(1000 .* [C2[jt].R[ik].Lu for ik in 1:K])
            # ϕs3 = log.(1000 .* [C3[jt].R[ik].ϕ for ik in 1:K]); Lus3 = log.(1000 .* [C3[jt].R[ik].Lu for ik in 1:K])
            # ϕs4 = log.(1000 .* [C4[jt].R[ik].ϕ for ik in 1:K]); Lus4 = log.(1000 .* [C4[jt].R[ik].Lu for ik in 1:K])

            ϕs =  [C[jt].R[ik].ϕ for ik in 1:K];   Lus= [C[jt].R[ik].Lu for ik in 1:K]
            ϕs2 = [C2[jt].R[ik].ϕ for ik in 1:K]; Lus2 =[C2[jt].R[ik].Lu for ik in 1:K]
            ϕs3 = [C3[jt].R[ik].ϕ for ik in 1:K]; Lus3 =[C3[jt].R[ik].Lu for ik in 1:K]
            ϕs4 = [C4[jt].R[ik].ϕ for ik in 1:K]; Lus4 =[C4[jt].R[ik].Lu for ik in 1:K]

            df1 = DataFrame(year = it, k = 1:K , ϵs = ϵs[1], θu = [pp[ik].θus[jt]*pp[ik].θprop for ik in 1:K] , Lu = Lus, ϕ= ϕs)
            df2 = DataFrame(year = it, k = 1:K , ϵs = ϵs[2], θu = [pp[ik].θus[jt]*pp[ik].θprop for ik in 1:K] , Lu = Lus2, ϕ= ϕs2)
            df3 = DataFrame(year = it, k = 1:K , ϵs = ϵs[3], θu = [pp[ik].θus[jt]*pp[ik].θprop for ik in 1:K] , Lu = Lus3, ϕ= ϕs3)
            df4 = DataFrame(year = it, k = 1:K , ϵs = ϵs[4], θu = [pp[ik].θus[jt]*pp[ik].θprop for ik in 1:K] , Lu = Lus4, ϕ= ϕs4)
            append!(df,df1)
            append!(df,df2)
            append!(df,df3)
            append!(df,df4)
            pl = plot(Lus, ϕs, title = "$it",
                      label = latexstring("\\epsilon = $(ϵs[1])"),
                      legend = :bottomright,
                      l = (:black,2),
                      m = (:circle, 5, :red),
                      yaxis = (L"\log \phi", :log10, (0.0001, 0.06)),
                      xaxis = (L"\log L_u", :log10, (0.001, 0.6)))
                      # yaxis = (L"\phi", (-1.,6)),
                      # xaxis = (L"L_u", (1.01,7)))
            plot!(pl,Lus2, ϕs2,
                      label = latexstring("\\epsilon = $(ϵs[2])"),
                      l = (:green,2),
                      m = (:circle, 5, :red))
            plot!(pl,  Lus3,ϕs3,
                     label = latexstring("\\epsilon = $(ϵs[3])"),
                     l = (:blue,2),
                     m = (:circle, 5, :red))
            plot!(pl,Lus4, ϕs4,
                   label = latexstring("\\epsilon = $(ϵs[4])"),
                   l = (:orange,2),
                   m = (:circle, 5, :red))

            frame(anim)
            # for ik in 1:K printfmtln(io,fmt,it,ik,ϵs[ik],pp[ik].θus[jt]*pp[ik].θprop,Lus[ik],ϕs[ik]) end

    end
    g = gif(anim, joinpath(dbplots,"phi-vs-Lu-eps.gif"),fps=0.5)
    df
end

"""
measure increase in land value at fringe. by default starts in 1960 and measures
increase up to 2020.

https://github.com/floswald/LandUse.jl/issues/11
"""
function issue11(;istart = 8, istop = 12, discount::Float64 = 0.03)
    (x,M,p) = run(Param())  # run standard single region

    # get location of fringe in start
    ϕstart = M[istart].ϕ
    ϕstop = M[istop].ϕ

    # land value in year `stop` at that distance
    setperiod!(p,istop)
    ρstart = ρ(ϕstart,p,M[istop])
    ρrstop = M[istop].ρr  # ρ at fringe in year stop

    # difference discounted back to 1990
    dd = (ρstart - ρrstop) / ((1 + discount)^5)^(istop-istart)
    out1 = dd / M[istart].ρr
    println("land value at fringe increased by $(round(out1,digits = 2)) percent \n from $(p.T[istart]) to $(p.T[istop]), discounted at $(100*discount) percent p.a.")

    # with flat epsilon
    (x,M,p) = run(Param(par = Dict(:ϵsmax => 0.0)))

    # get location of fringe in start
    ϕstart = M[istart].ϕ
    ϕstop = M[istop].ϕ

    # land value in year `stop` at that distance
    setperiod!(p,istop)
    ρstart = ρ(ϕstart,p,M[istop])
    ρrstop = M[istop].ρr  # ρ at fringe in year stop

    # difference discounted back to 1990
    dd = (ρstart - ρrstop) / ((1 + discount)^5)^(istop-istart)
    out2 = dd / M[istart].ρr
    println("with ϵ flat:")
    println("land value at fringe increased by $(round(out2,digits = 2)) percent \n from $(p.T[istart]) to $(p.T[istop]), discounted at $(100*discount) percent p.a.")
    out1,out2
end

function issue12_gif(n)
    us = range(1.01,1.04,length = n)
    ϵs = round.(range(1.0,5,length=8),digits=1)
    anim = Animation()
    for ie in ϵs
        x = issue12(ie,gu = us)
        frame(anim, x[6])
    end
    g = gif(anim, joinpath(dbplots,"issue12-$n.gif"),fps=0.5)
    g
end


# LandUse.runk(par =
#     Dict(1 => Dict(:θut=> 1.0, :θrt=> 1.0,:θu_g => 1.2, :θr_g => 1.2),
#          2 => Dict(:θut=> 1.0, :θrt=> 1.0,:θu_g => 1.19, :θr_g => 1.2)))


# multik(Dict(1 => Dict(:θut=> 1.0, :θrt=> 1.0,:θu_g => 1.2, :θr_g => 1.2), 2 => Dict(:θut=> 1.0, :θrt=> 1.0,:θu_g => 1.19, :θr_g => 1.2)))


"""
Runs K regions as a country and produces 2 plots
"""
function multik(ppar::Dict; save = false)
    K = length(ppar)

    cpar = Dict(:S => 1.0, :L => 1.0,
                :K => K,
                :kshare => [1/K for i in 1:K])

    sols,C,cpar,pp = LandUse.runk(cpar = cpar,par = ppar)
    pl = plot(impl_plot_ts_all(C)..., layout = (2,2))

    if save savefig(pl,joinpath(dbplots,"multi-$K-TS.pdf")) end

    pl2 = impl_plot_slopes(C)
                    # title = latexstring("City size vs pop \$\\epsilon\$ = $ϵ"))
    if save savefig(pl2,joinpath(dbplots,"multi-$K-phi-Lu.pdf")) end
    sols,C,cpar,pp,pl,pl2
end

# function issue12_1(ϵ; gu = [1.06,1.05], gr = 1.06)
#     cpar = Dict(:S => 1.0, :L => 1.0,
#                 :K => 2,
#                 :kshare => [1/2 for i in 1:2])
#
#     ppar = Dict(i => Dict(:ϵsmax => 0.0,:θagg => [0.32],
#                            :θrt => [0.32],:θut => [0.32], :θu_g => gu[i] ,
#                            :θr_g => gr , :ϵr => ϵ ,:ϵs => 0.0)
#                            for i in 1:2)
#     sols,C,cpar,pp = LandUse.runk(cpar = cpar,par = ppar)
#     pl = LandUse.plot_ts_all(C)
#     savefig(pl,joinpath(dbplots,"multi-2-TS.pdf"))
#
#     d = dataframe(C)
#     d.lϕ = log.(d.ϕ)
#     d.lu = log.(d.Lu)
#     gd = groupby(d,:year)
#     gd = combine(gd, AsTable([:lϕ, :lu]) => (x -> round.(diff(vcat(extrema(x.lu)...)) ./ diff(vcat(extrema(x.lϕ)...)),digits = 1)) => :slope)
#     d  = innerjoin(d,gd,on = :year)
#     transform!(d, AsTable([:year, :slope]) => (x -> string.(x.year) .* ": slope=" .* string.(x.slope) ) => :year_s)
#
#     cols = range(colorant"red",colorant"blue",length = length(pp[1].T))
#     dd = select(d, :year_s, :region, :lϕ, :lu)
#     pl2 = @df dd plot(:lϕ,:lu,group = :year_s,
#                     xlab = L"\log \phi",
#                     ylab = L"\log L_u",
#                     marker = (:circle, 4),
#                     colour = cols',
#                     legend = :topleft)
#     savefig(pl2,joinpath(dbplots,"multi-2-phi-Lu.pdf"))
#     sols,C,cpar,pp,pl,pl2
# end

"""
https://github.com/floswald/LandUse.jl/issues/22
"""
function issue_22()
    x,M,p = run(Region,Param())
    d = dataframe(M,p)
    df = @linq d |>
         select(:year,:Ch,:Cu,:Cr,:C ) |>
         transform(h = :Ch ./ :C,u = :Cu ./ :C,r = :Cr ./ :C) |>
         select(:year, :h, :u , :r)
    ds = stack(df, Not(:year))
    pl = @df ds plot(:year,:value, group = :variable,
               linewidth = 2, title = "Spending Shares")
    savefig(pl, joinpath(dbplots,"spending-shares.pdf"))
    pl
end


"""
https://github.com/floswald/LandUse.jl/issues/21
"""
function issue_21(n=8)

    es = range(0.5,8,length = n)
    e = es[1]
    x,M,p = run(Region,Param(par = Dict(:ϵsmax => 0.0, :ϵr => e)))
    # d = select(dataframe(M,p),:year,:Hr,:H0,:hr,:h0,:ρr,:qr,:ρ0,:q0, :ϕ, :Sr, :Srh)
    d = dataframe(M,p)
    d[!,:ϵr] .= e

    for (i,e) in enumerate(es[2:end])
        # println("e = $e")
        x,M,p = run(Region,Param(par = Dict(:ϵsmax => 0.0, :ϵr => e)))
        d0 = dataframe(M,p)
        d0[!,:ϵr] .= e
        append!(d,d0)
    end

    ds = stack(d,Not([:year,:ϵr]))
    # sims = [[:ϕ],[:q0, :ρ0] , [:qr, :ρr], [:Hr , :hr],[:H0 , :h0],[:Sr , :Srh]]
    sims = [:ϕ,:Sr,:Srh,:d0,
            :dr,:q0 ,:qr,:ρ0 ,
            :ρr, :H0,:Hr,:h0,
            :hr , :Lu, :Lr, :r]
    titles = ["phi"; "Sr"; "Srh"; "Central Density";
               "Fringe Density";"Central price q0"; "Fringe prices qr"; "Central land rho 0";
              "Fringe land rho r"; "Central H supply"; "Fringe H supply"; "Central H demand";
              "Fringe H demand"; "Lu"; "Lr"; "r"]

    plt = Any[]
    colors = setup_color(n)
    # colors = colormap("blues",n)
    for i in 1:length(sims)
        x = @linq ds |>
            where((:variable .== sims[i]))
            # where((:variable .∈ Ref(sims[i])))
        px = @df x plot(:year, :value, group = :ϵr,
                        title = titles[i],
                        titlefontsize=10,
                        color = colors',
                        # label=nms[i],
                        legend = false,
                        linewidth=2,marker = (:circle,3))
        push!(plt, px)
    end
    # plot(plt...,layout = (2,3))
    # ds

    plt_l = scatter(zeros(n), 1e-32 .* ones(n), zcolor = es,
            color = :inferno, colorbar = true, colorbar_title = L"\epsilon_r", legend = false,
            grid = false, axis = false, markerstrokealpha = 0, markersize = 0.1)
    psize=(1200,700)
    l = @layout [ [a b c d
                  e f g h
                  i j k m
                  q r s t] u{0.05w}]
    p = plot(plt..., plt_l, size=psize, layout = l)
    savefig(p,joinpath(dbplots,"varyepsr.pdf"))
    p

    # @df plot(:year,)
end


"""
https://github.com/floswald/LandUse.jl/issues/15
"""
function issue15()
    cpar = Dict(:S => 1.0, :L => 1.0,
                :K => 4,
                :θg => [1.1,1.0,0.98,0.96],
                :kshare => [0.25 for i in 1:4])
    sols,C,cpar,pp = LandUse.runk(cpar = cpar)
end

function fixed_ρ()
    p0=LandUse.Param(par = Dict(:ϵsmax =>0.0))
    LandUse.fixed_rho(p0, fi = "fixed-rho")
    p0=LandUse.Param(par = Dict(:θug => [1.2 for i in 1:14],
                  :θrg => [1.2 for i in 1:14],
                   :ϵsmax =>0.0))
    LandUse.fixed_rho(p0, fi = "fixed-rho-constg")
    p0=LandUse.Param(par = Dict(:θug => [1.2 for i in 1:14],
                  :θrg => [1.15 for i in 1:14],
                   :ϵsmax =>0.0))
    LandUse.fixed_rho(p0, fi = "fixed-rho-highu-g")
    p0=LandUse.Param(par = Dict(:θug => [1.15 for i in 1:14],
                  :θrg => [1.2 for i in 1:14],
                   :ϵsmax =>0.0))
    LandUse.fixed_rho(p0, fi = "fixed-rho-highr-g")

    p0=LandUse.Param(par = Dict(:ϵsmax =>0.0,:σ => 0.4))
    LandUse.fixed_rho(p0, fi = "fixed-rho-lowsig")
    p0=LandUse.Param(par = Dict(:θug => [1.2 for i in 1:14],
                  :θrg => [1.2 for i in 1:14],
                   :ϵsmax =>0.0,:σ => 0.4))
    LandUse.fixed_rho(p0, fi = "fixed-rho-constg-lowsig")
    p0=LandUse.Param(par = Dict(:θug => [1.2 for i in 1:14],
                  :θrg => [1.15 for i in 1:14],
                   :ϵsmax =>0.0,:σ => 0.4))
    LandUse.fixed_rho(p0, fi = "fixed-rho-highu-g-lowsig")
    p0=LandUse.Param(par = Dict(:θug => [1.15 for i in 1:14],
                  :θrg => [1.2 for i in 1:14],
                   :ϵsmax =>0.0,:σ => 0.4))
    LandUse.fixed_rho(p0, fi = "fixed-rho-highr-g-lowsig")

end

function fixed_rho(p::Param; fi = nothing)
    pyplot()

    x,Mr,p = run(Region,p)
    p.ρrbar = Mr[1].ρr
    x,Mu,p = run(Urban,p)

    dr = dataframe(Mr,p)
    du = dataframe(Mu,p)
    dr.model = ["baseline" for i in 1:nrow(dr)]
    du.model = ["urban" for i in 1:nrow(du)]

    # nms = [:year,:model,:area,:Lu,:ϕ, :ρr]
    # nms = [:year,:model,:area, :ρ0 ,:ϕ, :ρr]
    sims = [:Lu, :ϕ, :r, :q0, :ρ0 ,:qr, :ρr, :wr , :wu0]
    lat1 = latexstring("\$L_u\$; \$\\sigma\$=$(p.σ)")
    tis = [lat1, L"\phi", "r", L"q_0", L"\rho_0" ,L"q_r", L"\rho_r", "wr, gr=$(p.θrg[1])" ,"wu0, gu=$(p.θug[1])"]
    nms = [:year , :model , sims...]

    drs = stack(select(dr, nms, Not([:year,:model])))
    dus = stack(select(du, nms, Not([:year,:model])))
    dd = [drs; dus]

    pnms = sims
    plt = Any[]
    # tis = [:area :Lu :phi :rho_r]
    # tis = [:area :ρ0 :phi :rho_r]
    for i in 1:length(pnms)
        x = @linq dd |>
            where(:variable .== pnms[i])
        px = @df x plot(:year, :value, group = :model,
                        title = tis[i],
                        titlefontsize=12,
                        label = i == 1 ? ["baseline" "urban"] : "",
                        # legend = :topleft,
                        linewidth=2,marker = (:circle,4))
        push!(plt, px)
    end
    pl = plot(plt...)
    if isnothing(fi)
        fin = "fixed-rho.pdf"
    else
        fin = "$fi.pdf"
    end
    savefig(pl,joinpath(dbplots,fin))
    pl
end


function output_3Ms()

    # single region
    x,M,p = LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:ζ => 0.5, :τ1 => 0.98)))
    dd = LandUse.ts_plots(M,p)

    si = Dict()
    wihe = (700,340)
    si[:alloc] = plot(dd[:pop],dd[:spending], size = wihe)
    si[:space] = plot(dd[:Sr],dd[:phi], size = wihe)
    si[:h]     = plot(dd[:hr100],dd[:Hr100], size = wihe)
    si[:dens]  = plot(dd[:avdensity],dd[:densities], size = wihe)
    si[:r_rho] = plot(dd[:r_rho], size = wihe)
    si[:r_y] = plot(dd[:r_y], size = wihe)

    for (k,v) in si
        savefig(v, joinpath(dbplots,"$k.pdf"))
    end

    return si





    # multi regions
    mm = LandUse.multik(
               Dict(1 => Dict(:θut=> 1.01, :θrt=> 1.0,:θu_g => 1.31, :θr_g => 1.2, :ζ => 0.5),
                    2 => Dict(:θut=> 1.0, :θrt=> 1.0,:θu_g => 1.3, :θr_g => 1.2, :ζ => 0.5)))
    LandUse.savefig(mm[5],joinpath(LandUse.dbplots,"multi2.pdf"))
    # Labor Alloc and City size
    # p1 = plot(dd[:pop])

end

"""
    https://github.com/floswald/LandUse.jl/issues/36

1. same growth in sectors
    i. high cbar vs low sbar: show implied city density time series to see that only that config works
    ii. show falling housing spending share as well
2. implications of growth in either sector only
3. identify commuting cost params by matching time series data
"""
function issue36( ; save=false)
    r = Dict()

    emptyplot = plot(legend=false,grid=false,foreground_color_subplot=:white)
    size1by2 = (800,400)

    # 1. same growth in sectors
    #     i. high cbar vs low sbar: show implied city density time series to see that only that config works
    #     ii. show falling housing spending share as well
    r[1] = Dict()

    p1  = Param() # baseline param: high cbar and low sbar
    x,M,p0  = run(p1)
    pl1 = LandUse.ts_plots(M,p1)

    r[1][:baselinestats] = Dict(:ruralpop => Dict(:from => M[1].Lr , :to => M[end].Lr / p1.Lt[end]),
                                :rspend => Dict(:from => M[1].Cr / M[1].C, :to => M[end].Cr / M[end].C),
                                :phi => Dict(:from => M[1].ϕ, :to => M[end].ϕ ))
    r[1][:baseline1] = plot(pl1[:Lr_data],pl1[:pop],pl1[:n_densities],pl1[:r_y],
                                layout = (2,2))
    r[1][:baseline2] = plot(pl1[:avdensity_pop],pl1[:n_densities],layout = (1,2),
                            size = size1by2)

    r[1][:prices] = LandUse.plot(pl1[:r_real],pl1[:r_y], layout = (1,2))
    r[1][:baseline_fit] = plot(pl1[:Lr_data],pl1[:spending],pl1[:n_densities],pl1[:avdensity],
                                layout = (2,2)) #




    p2 = Param(par = Dict(:cbar => 0.4, :sbar => 0.7)) # low cbar and high sbar
    x,M,p0  = run(p2)
    pl2 = LandUse.ts_plots(M,p2)
    r[1][:low_cbar] = plot(pl2[:Lr_data],pl2[:spending],pl2[:avdensity],pl2[:r_y],
                                layout = (2,2),link = :x)


    # 2. implications of growth in either sector only
    # i. u grows faster than r
    r[2] = Dict()
    p3 = LandUse.Param(par = Dict(:θrt => 1.0,:θr_g => 1.05))
    # p2 = LandUse.Param(par = Dict(:θu_g => 1.08,:θut => 1.0, :θrt => 1.0,:θr_g => 1.01))
    x,M,p0  = LandUse.run(p3)
    pl3 = LandUse.ts_plots(M,p3)
    r[2][:u_fast1] = LandUse.plot(pl3[:Lr_data],pl3[:spending],pl3[:r_y],emptyplot,
                                layout = (2,2))
    r[2][:u_fast2] = LandUse.plot(pl3[:avdensity_pop],plot(pl3[:n_densities],leg = :bottomleft),layout = (1,2),
                            size = size1by2)

    # p3 = LandUse.Param(par = Dict(:θu_g => 1.01,:θut => 1.0, :θrt => 1.0,:θr_g => 1.09))
    p4 = LandUse.Param(par = Dict(:θu_g => 1.0,:θut => 1.0))
    x,M,p0  = LandUse.run(p4)
    pl4 = LandUse.ts_plots(M,p4)
    r[2][:r_fast1] = LandUse.plot(pl4[:Lr_data],pl4[:spending],pl4[:r_y],emptyplot,
                                layout = (2,2))
    r[2][:r_fast2] = LandUse.plot(pl4[:avdensity_pop],plot(pl4[:n_densities],leg = :bottomleft),layout = (1,2),
                            size = size1by2)

    r[3] = plot(pl1[:mode],pl1[:ctime],layout = (2,1), link = :x)

    # 4. constant growth in both sectors
    r[4] = Dict()
    p = LandUse.Param(par = Dict(:θu_g => 1.25,:θut => 1.0, :θrt => 1.0,:θr_g => 1.25))
    x,M,p0  = LandUse.run(p)
    pl5 = LandUse.ts_plots(M,p)
    r[4] = LandUse.plot(pl5[:Lr_data],pl5[:spending],pl5[:avdensity],pl5[:r_y],
                                layout = (2,2),link = :x)



    # save plots
    if save
        savefig(r[1][:baseline1], joinpath(dbplots,"issue36-baseline.pdf"))
        savefig(r[1][:baseline_fit], joinpath(dbplots,"issue36-baseline-fit.pdf"))

        savefig(r[1][:low_cbar], joinpath(dbplots,"issue36-low-cbar.pdf"))
        savefig(r[1][:prices], joinpath(dbplots,"issue36-prices.pdf"))
        savefig(r[2][:u_fast1], joinpath(dbplots,"issue36-u-fast.pdf"))
        savefig(r[2][:r_fast1], joinpath(dbplots,"issue36-r-fast.pdf"))
        savefig(r[4], joinpath(dbplots,"issue36-constant.pdf"))
        savefig(r[3], joinpath(dbplots,"issue36-commute.pdf"))

        # changing size
        savefig(r[1][:baseline2], joinpath(dbplots,"issue36-baseline2.pdf"))
        savefig(r[2][:u_fast2], joinpath(dbplots,"issue36-u-fast2.pdf"))
        savefig(r[2][:r_fast2], joinpath(dbplots,"issue36-r-fast2.pdf"))



    end
    return r

end


"""
produces all output from the model needed to compile the paper.

This uses current baseline parameters defined in `params.json`
"""
function output_paper(;save = false)

    p1  = Param() # baseline param: high cbar and low sbar
    x,M,p0  = run(Region,p1)
    pl1 = LandUse.ts_plots(M,p1)

    # names of plots we want
    key = [:pop,:spending,:avdensity,:phi,:densities,:r_y,:r_rho,:q0_qr,:mode,:ctime]

    # subset to those
    r = filter( x -> x.first ∈ key , pl1)

    # save if wanted
    if save
        for (k,v) in r
            savefig(v, joinpath(dbplots,"$k.pdf"))
        end
    end
    r
end
