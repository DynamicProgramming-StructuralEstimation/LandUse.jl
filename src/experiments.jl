

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
plot ``\\phi_k`` vs ``\\L{u,k}`` for all regions ``k`` and different eps slopes.
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
                      yaxis = (L"\log \phi", :log, (0.0001, 0.06)),
                      xaxis = (L"\log L_u", :log, (0.001, 0.6)))
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



function issue12(; gf = [1.0,1.02,1.04])
    cpar = Dict(:S => 1.0, :L => 1.0,
                :K => 3,
                :kshare => [1/3 for i in 1:3])

    ppar = Dict(i => Dict(:ϵsmax => 0.0, :θut => [gf[i] for ti in 1:14],:ϵr => 2.0) for i in 1:3)
    sols,C,cpar,pp = LandUse.runk(cpar = cpar,par = ppar)
    pl = LandUse.plot_ts_all(C)
    savefig(pl,joinpath(dbplots,"multi-3-TS.pdf"))

    d = dataframe(C)
    d.lϕ = log.(d.ϕ)
    d.lu = log.(d.Lu)
    dd = select(d, :year, :region, :lϕ, :lu)
    pl2 = @df dd plot(:lϕ,:lu,group = :year,
                    xlab = L"\log \phi",
                    ylab = L"\log L_u",
                    marker = (:circle, 4),
                    legend = :bottomright)
    savefig(pl2,joinpath(dbplots,"multi-3-phi-Lu.pdf"))
    sols,C,cpar,pp,pl,pl2
end

function issue12_1(ϵ; gf = [1.0,1.06])
    cpar = Dict(:S => 1.0, :L => 1.0,
                :K => 2,
                :kshare => [1/2 for i in 1:2])

    ppar = Dict(i => Dict(:ϵsmax => 0.0, :θut => [gf[i] for ti in 1:14], :ϵr => ϵ) for i in 1:2)
    sols,C,cpar,pp = LandUse.runk(cpar = cpar,par = ppar)
    pl = LandUse.plot_ts_all(C)
    savefig(pl,joinpath(dbplots,"multi-2-TS.pdf"))

    d = dataframe(C)
    d.lϕ = log.(d.ϕ)
    d.lu = log.(d.Lu)
    dd = select(d, :year, :region, :lϕ, :lu)
    pl2 = @df dd plot(:lϕ,:lu,group = :year,
                    xlab = L"\log \phi",
                    ylab = L"\log L_u",
                    marker = (:circle, 4),
                    legend = :bottomright)
    savefig(pl2,joinpath(dbplots,"multi-2-phi-Lu.pdf"))
    sols,C,cpar,pp,pl,pl2
end

"""
https://github.com/floswald/LandUse.jl/issues/21
"""
function issue_21(n=8)

    es = range(0.5,8,length = n)
    e = es[1]
    x,M,p = run(Region,Param(par = Dict(:ϵsmax => 0.0, :ϵr => e)))
    # d = select(dataframe(M,p.T),:year,:Hr,:H0,:hr,:h0,:ρr,:qr,:ρ0,:q0, :ϕ, :Sr, :Srh)
    d = dataframe(M,p.T)
    d[!,:ϵr] .= e

    for (i,e) in enumerate(es[2:end])
        # println("e = $e")
        x,M,p = run(Region,Param(par = Dict(:ϵsmax => 0.0, :ϵr => e)))
        d0 = dataframe(M,p.T)
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
