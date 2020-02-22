

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
        pl = plot(ϕs, Lus, title = "$it",xlab = L"\phi", ylab = L"L_u",
                  leg = false,
                  l = (:black,2),
                  m = (:circle, 5, :red),
                  xlims = (0,0.05),
                  ylims = (0,0.3))
        frame(anim)
    end
    g = gif(anim, joinpath(dbpath,"phi-vs-Lu.gif"),fps=0.5)

end
