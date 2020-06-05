









function works()
	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.99, :ζ => 0.1)))
	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.99, :ζ => 0.2)))
	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.99, :ζ => 0.3)))
	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.99, :ζ => 0.4)))
	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.98, :ζ => 0.5)))


	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.4, :ζ => 0.2)))
	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.5, :ζ => 0.2)))
	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.6, :ζ => 0.2)))
	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.7, :ζ => 0.2)))
	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.8, :ζ => 0.2)))


	LandUse.run(LandUse.Region,LandUse.Param(par = Dict(:τ1 => 0.9, :ζ => 0.3)))
end

function iccost()
	sbs = 0.0:0.05:0.5
	cbs = 0.0:0.05:0.8
	sis = 0.1:0.1:0.99
	θus = 1.0:0.05:1.3
	θrs = 1.0:0.05:1.3
	esl = 0.0:10:100.0
	eps = 0.0:0.1:6.0
	t1s = 0.5:0.1:1.0
	zetas = 0.0:0.1:0.5

	p = Param()
	@manipulate for sb in slider(sbs, label = "sbar", value = p.sbar),
					cb in slider(cbs, label = "cbar", value = p.cbar),
					# sig in slider(sis, label = "σ", value = p.σ),
					# # agg in slider(ags, label = "agg_g", value = p.θagg_g),
					# θug in slider(θus, label = "θu_g", value = p.θu_g),
					# esl in slider(esl, label = "ϵ-slope", value = p.ϵsmax),
					# epl in slider(eps, label = "ϵ", value = p.ϵr),
					# θrg in slider(θrs, label = "θr_g", value = p.θr_g),
					t1 in slider(t1s, label = "τ1", value = 0.9),
					# t1 in slider(range(0.99,0.99,length = 1)),
					zeta in slider(zetas, label = "ζ", value = p.ζ)
					# zeta in slider(range(0.1,0.1,length = 1))
			# t1 = 0.99
			x,M,p = try
				run(Region,
				             # Param(par = Dict(:ζ => zeta,:τ1 => t1,:ϵs => 0.0, :ϵsmax => esl,:ϵr => epl,
							 #                  :sbar => sb, :cbar => cb, :σ => sig)))
							 Param(par = Dict(:ζ => zeta, :τ1 => t1, :sbar => sb, :cbar => cb)))


			catch
				# println("error with ζ = $zeta, τ1 = $t1")
				(0,0,0)
			end

		if isa(M,Vector)
		   	dd = ts_plots(M,p)
			plot(dd[:pop],dd[:avdensity],dd[:phi],dd[:tauphi] ,layout = (1,4),size = (900,250))
		else
			println("error with ζ = $zeta, τ1 = $t1")
		end
	end
end

"""
https://github.com/floswald/LandUse.jl/issues/22
"""
function i22()
	sbs = 0.0:0.05:0.5
	cbs = 0.0:0.05:0.5
	sis = 0.1:0.1:0.99
	# ags = 0.0:0.1:1.0
	θus = 1.0:0.05:1.3
	θrs = 1.0:0.05:1.3
	esl = 0.0:10:100.0
	eps = 0.0:0.1:6.0
	t1s = range(0.6,stop = 0.98, length = 5)
	zetas = 0.05:0.05:0.5


	p = Param()
	@manipulate for θtype in Dict("Matlab θ" => 1, "Growth θ" => 2),
		            # sb in slider(sbs, label = "sbar", value = p.sbar),
		            # cb in slider(cbs, label = "cbar", value = p.cbar),
		            # sig in slider(sis, label = "σ", value = p.σ),
		            # # agg in slider(ags, label = "agg_g", value = p.θagg_g),
		            θug in slider(θus, label = "θu_g", value = p.θu_g),
		            # esl in slider(esl, label = "ϵ-slope", value = p.ϵsmax),
		            # epl in slider(eps, label = "ϵ", value = p.ϵr),
		            θrg in slider(θrs, label = "θr_g", value = p.θr_g),
		            t1 in slider(t1s, label = "τ1", value = 0.98),
		            # t1 in slider(range(0.99,0.99,length = 1)),
		            zeta in slider(zetas, label = "ζ", value = p.ζ)
		            # zeta in slider(range(0.1,0.1,length = 1))


		if θtype == 1
			x,M,p = try
				run(Region,
				             # Param(par = Dict(:ζ => zeta,:τ1 => t1,:ϵs => 0.0, :ϵsmax => esl,:ϵr => epl,
							 #                  :sbar => sb, :cbar => cb, :σ => sig)))
							 Param(par = Dict(:ζ => zeta,:τ1 => t1)))


			catch
				println("error with ζ = $zeta, τ1 = $t1")
				(0,0,0)
			end
		else
			x,M,p = try
	    		run(Region,
		             # Param(par = Dict(:ζ => zeta,:τ1 => t1,:ϵs => 0.0, :ϵsmax => esl,:ϵr => epl,
					 #                  :sbar => sb, :cbar => cb, :σ => sig,
						# 			  :θut => 0.32, :θrt => 0.32, :θu_g => θug, :θr_g => θrg)))
						Param(par = Dict(:θut => 0.32, :θrt => 0.32, :θu_g => θug, :θr_g => θrg,:ζ => zeta,:τ1 => t1)))



		    catch
				println("error with ζ = $zeta, τ1 = $t1")
				(0,0,0)

			end

        end
		# @info "model done."

		if isa(M,Vector)
			dd = ts_plots(M,p)
			plot(dd[:pop],dd[:avdensity],dd[:phi],dd[:tauphi] ,layout = (1,4),size = (900,250))
		end
    end
end


function i0()
	epsimax = 0.0:10.0
	gfac = 1.00:0.05:2.0
	sigmas = 0.1:0.1:0.99
	sbs = 0.0:0.1:0.5

	@manipulate for growthtype = Dict("orig" => 1,"u-r const" => 4),
					ugrowth in slider(gfac,value = 1.01, label = "u-growthrate"),
					rgrowth in slider(gfac,value = 1.01, label = "r-growthrate"),
					epsim in slider(epsimax, value = 1.0, label = "ϵr"),
					sigma in slider(sigmas, value = 0.9, label = "sigma"),
					sb in slider(sbs, label = "sbar")

					if growthtype == 1
						p0 = LandUse.Param(par = Dict(:ϵr => epsim,:ϵsmax => 0.0,:σ => sigma , :sbar => sb))
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
						p0 = LandUse.Param(par = Dict(:θut => [ugrowth for i in 1:14],
											 :θrt => [rgrowth for i in 1:14],
						                     :ϵr => epsim,:ϵsmax => 0.0,
											 :σ => sigma, :sbar => sb ))
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
						p0 = LandUse.Param(par = Dict(:θut => [ugrowth for i in 1:14],
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

function interact1()
	K = 2
	es = 1.0:0.5:10.0
	g2 = 1.0:0.01:1.06

    # mp = @manipulate for θr in slider(θrs, label = "θr"), θu in slider(θus, label = "θu")
    mp = @manipulate for e in slider(es, label = "ϵr", value = 4.0 ),
		                 g in slider(g2, label = "growth2", value = 1.0)
		x = LandUse.issue12_1(e; gf = [1.0,g])
		x[5]
    end
end

function interact2()
	K = 2
	es = 1.0:0.5:10.0
	g2 = 1.0:0.01:1.06

    # mp = @manipulate for θr in slider(θrs, label = "θr"), θu in slider(θus, label = "θu")
    mp = @manipulate for e in slider(es, label = "ϵr", value = 4.0 ),
		                 g in slider(g2, label = "growth2", value = 1.01)
		x = LandUse.issue12_1(e; gf = [1.0,g])
		x[5]
    end
end


function i12()
	g1 = 1.0:0.005:1.04
	esl = 0.0:20:400.0
	es = 1.0:0.5:10.0
	esm = 1.0:0.5:10.0
	zetas = 0.0:0.1:0.99
	t1s = 0.5:0.1:1.0

	mp = @manipulate for e in slider(es, label = "ϵr", value = 4.0 ),
		                 esl in slider(esl,label = "eslope",value = 0),
		                 # gr in slider(g1,label = "gr", value = 4.0),
		                 gs1 in slider(g1,label = "ug1",value = 1.019),
		                 gs2 in slider(g1,label = "ug2",value = 1.019),
		                 gs3 in slider(g1,label = "ug3",value = 1.019),
		                 # gagg in slider(g1,label = "agg"),
						 t1 in slider(t1s, label = "τ1", value = 0.9),

						 z in slider(zetas,label = "ζ", value = 0.0)

						 x = issue12(e; gu = [gs1,gs2,gs3], θu = [1.0,1.0,1.0],θagg_g = 0.0, ϵsmax = esl, zeta = z)
						 # plot(rand(10))
						 # plot(x[5],x[6],layout = (1,2),size = (900,600))
						 x[5]
					 end

end
