

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

    anim = Animation()
    for (jt,it) in enumerate(C[1].T)
        ϕs = [C[jt].R[ik].ϕ for ik in 1:K]; Lus = [C[jt].R[ik].Lu for ik in 1:K]
        pl = plot(ϕs, Lus, title = "$it",
                  leg = false,
                  l = (:black,2),
                  m = (:circle, 5, :red),
                  xaxis = (L"\phi", (0.001,0.045)),
                  yaxis = (L"L_u", (0.01,0.3)))
        frame(anim)
    end
    g = gif(anim, joinpath(dbpath,"phi-vs-Lu.gif"),fps=0.5)

end

"""
plot ``\\phi_k`` vs ``\\L{u,k}`` for all regions ``k`` and different eps slopes.
relates to https://github.com/floswald/LandUse.jl/issues/10
"""
function issue10()
    cpar = Dict(:S => 1.0, :L => 1.0,
                :K => 4,
                :θprop => [1.0,0.99,0.98,0.97],
                :kshare => [0.25 for i in 1:4])
    sols,C,cp,p0 = LandUse.runk(cpar = cpar)

    K = C[1].K
    ϵ1 = 5.0
    ϵ2 = 1.0
    ϵ3 = 0.0

    # version with lower eps slope
    sols,C2,cp,pp = LandUse.runk(cpar = cpar,
                                 par = Dict(ik => Dict(:ϵsmax => ϵ1) for ik in 1:K))

    # version with lower eps slope
    sols,C3,cp,pp = LandUse.runk(cpar = cpar,
                                 par = Dict(ik => Dict(:ϵsmax => ϵ2) for ik in 1:K))

    # version with lower eps slope
    sols,C4,cp,pp = LandUse.runk(cpar = cpar,
                              par = Dict(ik => Dict(:ϵsmax => ϵ3) for ik in 1:K))



    anim = Animation()
    for (jt,it) in enumerate(C[1].T)
        ϕs = [C[jt].R[ik].ϕ for ik in 1:K]; Lus = [C[jt].R[ik].Lu for ik in 1:K]
        ϕs2 = [C2[jt].R[ik].ϕ for ik in 1:K]; Lus2 = [C2[jt].R[ik].Lu for ik in 1:K]
        ϕs3 = [C3[jt].R[ik].ϕ for ik in 1:K]; Lus3 = [C3[jt].R[ik].Lu for ik in 1:K]
        ϕs4 = [C4[jt].R[ik].ϕ for ik in 1:K]; Lus4 = [C4[jt].R[ik].Lu for ik in 1:K]
        pl = plot(ϕs, Lus, title = "$it",
                  label = latexstring("\\epsilon = $(p0[1].ϵs)"),
                  legend = :bottomright,
                  l = (:black,2),
                  m = (:circle, 5, :red),
                  xaxis = (L"\phi", (0.001,0.045)),
                  yaxis = (L"L_u", (0.01,0.33)))
        plot!(pl, ϕs2, Lus2,
                  label = latexstring("\\epsilon = $(ϵ1)"),
                  l = (:green,2),
                  m = (:circle, 5, :red))
        plot!(pl, ϕs3, Lus3,
                 label = latexstring("\\epsilon = $(ϵ2)"),
                 l = (:blue,2),
                 m = (:circle, 5, :red))
        plot!(pl, ϕs4, Lus4,
               label = latexstring("\\epsilon = $(ϵ3)"),
               l = (:orange,2),
               m = (:circle, 5, :red))

        frame(anim)
    end
    g = gif(anim, joinpath(dbpath,"phi-vs-Lu-eps$ϵ.gif"),fps=0.5)

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
