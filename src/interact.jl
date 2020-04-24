

using Interact

function i0()
	epsimax = 0.0:10.0
	gfac = 1.00:0.05:1.25
	sigmas = 0.1:0.1:0.99

	@manipulate for growthtype = Dict("orig" => 1, "orig-u" => 2, "orig-r" => 3, "u-r const" => 4),
					ugrowth in slider(gfac,value = 1.01, label = "u-growthrate"),
					rgrowth in slider(gfac,value = 1.01, label = "r-growthrate"),
					epsim in slider(epsimax, value = 0.0, label = "ϵ-slope"),
					sigma in slider(sigmas, value = 0.9, label = "sigma")

					if growthtype == 1
						p0 = LandUse.Param(par = Dict(:ϵsmax => epsim,:σ => sigma ))
						x,M,p = LandUse.run(LandUse.Region,p0)
						LandUse.plot_ts(M,p0)

					elseif growthtype == 2
						println("does not work")
						# p0 = LandUse.Param(par = Dict(:θug => [ugrowth for i in 1:14],
						# 					 :ϵsmax => epsim,
						# 					 :σ => sigma ))
						# x,M,p = LandUse.run(LandUse.Region,p0)
						# LandUse.plot_ts(M,p0)
					#
					elseif growthtype == 3
						println("does not work")
					# 		p0 = LandUse.Param(par = Dict(:θrg => [rgrowth for i in 1:14],
					# 							 :ϵsmax => epsim,
					# 							 :σ => sigma ))
					# 		x,M,p = LandUse.run(LandUse.Region,p0)
					# 		LandUse.plot_ts(M,p0,time)
					#
					elseif growthtype == 4
						p0 = LandUse.Param(par = Dict(:θug => [ugrowth for i in 1:14],
											 :θrg => [rgrowth for i in 1:14],
						                     :ϵsmax => epsim,
											 :σ => sigma ))
					    x,M,p = LandUse.run(LandUse.Region,p0)
						LandUse.plot_ts(M,p0)
					end
	end
end


function i1()

	ti = 1:14
	epsimax = 0.0:10.0
	gfac = 1.00:0.05:1.25
	sigmas = 0.1:0.1:0.99

	@manipulate for growthtype = Dict("orig" => 1, "orig-u" => 2, "orig-r" => 3, "u-r const" => 4),
		            time in slider(ti, value = 14, label = "time"),
					ugrowth in slider(gfac,value = 1.24, label = "u-growthrate"),
					rgrowth in slider(gfac,value = 1.24, label = "r-growthrate"),
					epsim in slider(epsimax, value = 0.0, label = "ϵ-slope"),
					sigma in slider(sigmas, value = 0.9, label = "sigma")

					if growthtype == 1
						p0 = LandUse.Param(par = Dict(:ϵsmax => epsim,:σ => sigma ))
						x,M,p = LandUse.run(LandUse.Region,p0)
						LandUse.plot_ts(M,p0,time)

					elseif growthtype == 2
						p0 = LandUse.Param(par = Dict(:θug => [ugrowth for i in 1:14],
											 :ϵsmax => epsim,
											 :σ => sigma ))
						x,M,p = LandUse.run(LandUse.Region,p0)
						LandUse.plot_ts(M,p0,time)
					#
					# elseif growthtype == 3
					# 		p0 = LandUse.Param(par = Dict(:θrg => [rgrowth for i in 1:14],
					# 							 :ϵsmax => epsim,
					# 							 :σ => sigma ))
					# 		x,M,p = LandUse.run(LandUse.Region,p0)
					# 		LandUse.plot_ts(M,p0,time)
					#
					# elseif growthtype == 4
					# 	p0 = LandUse.Param(par = Dict(:θug => [ugrowth for i in 1:14],
					# 						 :θrg => [rgrowth for i in 1:14],
					# 	                     :ϵsmax => epsim,
					# 						 :σ => sigma ))
					#     x,M,p = LandUse.run(LandUse.Region,p0)
					# 	LandUse.plot_ts(M,p0,time)


					# elseif growthtype == 3
					# 	gg = (LandUse.originalθ[end] - LandUse.originalθ[1]) / 13
					# 	g0 = LandUse.originalθ[1] .+ gg .* [i for i in 0:13]
					# 	g = g0[2:end] ./ g0[1:end-1]
					# 	p0 = LandUse.Param(par = Dict(:θug => [growth for i in 1:14],
					# 						 :θrg => [growth for i in 1:14],
					# 	                     :ϵsmax => epsim))
					#     x,M,p = LandUse.run(p0)
					# 	LandUse.plot_ts_xsect(M,p0,time)

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
