

using Interact

function i1()

	ti = 1:14
	epsimax = 3:10.0
	gfac = 1.03:0.05:1.5

	@manipulate for growthtype = Dict("orig" => 1, "linear" => 2),
		            time in ti,
					growth in gfac,
					epsim in slider(epsimax, value = 10.0, label = "ϵ-slope")

					if growthtype == 2
						p0 = LandUse.Param(par = Dict(:θug => [growth for i in 1:14],
											 :θrg => [growth for i in 1:14],
						                     :ϵsmax => epsim))
					    x,M,p = LandUse.run(p0)
						LandUse.plot_ts_xsect(M,p0,time)
					else
						p0 = LandUse.Param(par = Dict(:ϵsmax => epsim))
						x,M,p = LandUse.run(p0)
						LandUse.plot_ts_xsect(M,p0,time)
					end

	end
end

function interact2()
	K = 2
	gs = 1.0:0.01:1.5
	ks = 0.0:0.1:1.0

    # mp = @manipulate for θr in slider(θrs, label = "θr"), θu in slider(θus, label = "θu")
    mp = @manipulate for g in slider(gs, label = "growth2", value =1.0 ),
		                 k in slider(ks, label = "k-share2", value = 0.5)
		cpar = Dict(:S => 1.0, :L => 1.0,
					:K => 2,
					:θg => [1.0,g],
					:kshare => [(1-k), k])
		sols,C,cpar,pp = LandUse.runk(cpar = cpar)
		try
			LandUse.plot_ts_all(C)
		catch
			println("infeasible")
		end
		# catch
		# 	println("infeasible")
		# end
    end
end
