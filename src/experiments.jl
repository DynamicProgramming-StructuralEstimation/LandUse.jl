

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
