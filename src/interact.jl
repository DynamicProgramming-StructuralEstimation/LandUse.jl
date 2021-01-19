









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
	gfac = 1.00:0.05:2.0
	eta = 0.5:0.01:2.0
	tau = 0.1:0.1:1.0
	eta2 = 0.0:0.01:1.0
	cbars = 0.0:0.01:1.5
	ctaus = 0.0:0.01:4.0
	psis = 0.1:0.01:1.0
	p1 = Param()

	# mp = @manipulate for est in OrderedDict("from data θ" => 1, "estimate θ" => 2, "from estimation θ" => 3),
	mp = @manipulate for it in slider(1:length(p1.T), value = length(p1.T), label = "period"),
					ϕx in slider(0.05:0.05:1.0, value = 1.0, label = "ϕ1x"),
					cbar in slider(cbars, value = p1.cbar, label = "cbar"),
					sbar in slider(cbars, value = p1.sbar, label = "sbar"),
					gamma in slider(psis, value = p1.γ, label = "gamma"),
					cτ in slider(ctaus, value = p1.a, label = "cτ"),
					es in slider(0.0:0.1:10.0, value = 0.0, label = "ϵs"),
					etam in slider(eta, value = p1.ηm, label = "ηm"),
					etal in slider(eta2, value = p1.ηl, label = "ηl"),
					etaw in slider(eta2, value = p1.ηw, label = "ηw")

					p0 = LandUse.Param(par = Dict(:ηm => etam, :ηl => etal, :ηw => etaw,
											  :cbar => cbar, :sbar => sbar, :cτ => cτ, :γ => gamma, :ϵsmax => es, :ϵs => es,
											  :ϕ1x => ϕx),
									   use_estimatedθ = false)

					try
						x,M,p = LandUse.run(p0, estimateθ = false)
						pl = LandUse.ts_plots(M,p0,fixy = false)
						pc = LandUse.cs_plots(M[it], p0, it)
						plot(pl[:Lr_data],pl[:spending],pl[:pr_data],pl[:productivity],
						     pl[:n_densities], pl[:densities], pl[:mode], pl[:ctime],
							 pl[:phi] , pl[:qbar_real], pl[:r_y], pl[:r_rho],
							 pc[:ϵ] , pc[:D], pc[:q] , pc[:H],
							 layout = (4,4),size = (1200,800))
					catch e
						wdg = alert("Error!")
						print(wdg())
					end

					# plot(pl[:phi], pl[:avdensity],pl[:mode],pl[:ctime],pl[:dist_vs_time],plot(), l = (2,3))
					# plot(pl[:Lr_data],pl[:spending],pl[:qbar_real],pl[:productivity],pl[:n_densities], pl[:avdensity], layout = (2,3),link = :x)

	end
	@layout! mp vbox(hbox(:it, :ϕx),
	                 hbox(:sbar, :cbar,:gamma,:es),
	                 hbox(:cτ, :etam, :etal, :etaw),
					 observe(_))
end

function ik()
	p1 = Param()
	gfac = 1.00:0.05:2.0
	eta = 0.5:0.01:3.0
	tau = 0.1:0.1:1.0
	eta2 = 0.0:0.01:1.0
	cbars = 0.0:0.01:1.5
	ctaus = 0.0:0.01:6.0
	psis = 0.1:0.01:1.0
	mp = @manipulate for cbar in slider(cbars, value = p1.cbar, label = "cbar"),
					sbar in slider(cbars, value = p1.sbar, label = "sbar"),
					gamma in slider(psis, value = p1.γ, label = "gamma"),
					cτ in slider(ctaus, value = p1.a, label = "cτ"),
					es in slider(0.0:0.1:10.0, value = 0.0, label = "ϵs"),
					g2 in slider(1.0:0.01:1.1, value = 1.02, label = "g2"),
					etam in slider(eta, value = p1.ηm, label = "ηm"),
					etal in slider(eta2, value = p1.ηl, label = "ηl"),
					etaw in slider(eta2, value = p1.ηw, label = "ηw")

					try
						x,C,p = runk(par = Dict(:K => 2, :kshare => [0.5,0.5], :factors => [1.0,g2],
						                        :ηm => etam, :ηl => etal, :ηw => etaw,
												  :cbar => cbar, :sbar => sbar, :cτ => cτ, :γ => gamma, :ϵsmax => es))
						x = impl_plot_slopes(C)
						plot(x[1],x[2],x[3],
							 layout = (1,3),size = (1200,600))
					catch e
						wdg = alert("Error!")
						print(wdg())
					end

					# plot(pl[:phi], pl[:avdensity],pl[:mode],pl[:ctime],pl[:dist_vs_time],plot(), l = (2,3))
					# plot(pl[:Lr_data],pl[:spending],pl[:qbar_real],pl[:productivity],pl[:n_densities], pl[:avdensity], layout = (2,3),link = :x)

	end
	@layout! mp vbox(hbox(:g2),
	                 hbox(:sbar, :cbar,:gamma,:es),
	                 hbox(:cτ, :etam, :etal, :etaw),
					 observe(_))


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


function imulti3()
	g1 = 1.18:0.01:1.22
	esl = 0.0:20:400.0
	es = 1.0:0.5:10.0
	esm = 1.0:0.5:10.0
	zetas = 0.0:0.1:0.6
	t1s = 0.0:0.1:1.0
	t2s = 0.5:0.1:2.0

	K = 3
	cpar = Dict(:S => 1.0, :L => 1.0,
				:K => K,
				:kshare => [1/K for i in 1:K])


	mp = @manipulate for
		                 # esl in slider(esl,label = "eslope",value = 0.0),
		                 # gr in slider(g1,label = "gr", value = 4.0),
		                 gs1 in slider(g1,label = "u-growth 1",value = 1.19),
		                 gs2 in slider(g1,label = "u-growth 2",value = 1.2),
		                 gs3 in slider(g1,label = "u-growth 3",value = 1.21),
						 etam in slider(t2s, label = "ηm", value = 1.0),
						 etal in slider(t1s, label = "ηl", value = 0.0)
						 # z in slider(zetas,label = "ζ", value = 0.0)

						 d0 = Dict(:ηm => etam, :ηl => etal)

						 dd = Dict(1 => merge(Dict(:θut=> 1.0, :θrt=> 1.0,:θu_g => gs1,  :θr_g => 1.2), d0),
						 		   2 => merge(Dict(:θut=> 1.0, :θrt=> 1.0,:θu_g => gs2, :θr_g => 1.2),d0),
								   3 => merge(Dict(:θut=> 1.0, :θrt=> 1.0,:θu_g => gs3,  :θr_g => 1.2),d0)
								   )

   					    sols,C,cp,pp = LandUse.runk(cpar = cpar,par = dd)
   					    p1 = plot(impl_plot_ts_all(C)..., layout = (2,2))
   					    p2 = impl_plot_slopes(C)

						p1
	# 					p1
	end
end

function imulti2()
	g1 = 1.0:0.01:1.1
	esl = 0.0:20:400.0
	es = 1.0:0.5:10.0
	esm = 1.0:0.5:10.0
	zetas = 0.0:0.1:0.6
	t1s = 0.5:0.1:1.0

	K = 2
	cpar = Dict(:S => 1.0, :L => 1.0,
				:K => K,
				:kshare => [1/K for i in 1:K])


	mp = @manipulate for e in slider(es, label = "ϵr", value = 4.0 ),
		                 esl in slider(esl,label = "eslope",value = 0.0),
		                 # gr in slider(g1,label = "gr", value = 4.0),
		                 gs1 in slider(g1,label = "u-shift 1",value = 1.0),
		                 gs2 in slider(g1,label = "u-shift 2",value = 1.0),
		                 # gs3 in slider(g1,label = "u-shift 3",value = 1.0),
						 t1 in slider(t1s, label = "τ1", value = 0.9),
						 z in slider(zetas,label = "ζ", value = 0.0)

						 d0 = Dict(:ζ => z, :τ1 => t1, :ϵr => e, :ϵsmax => esl)

						 dd = Dict(1 => merge(Dict(:θut=> gs1, :θrt=> 1.0,:θu_g => 1.2,  :θr_g => 1.2), d0),
						 		   2 => merge(Dict(:θut=> gs2, :θrt=> 1.0,:θu_g => 1.2, :θr_g => 1.2),d0)
								   # 3 => merge(Dict(:θut=> gs3, :θrt=> 1.0,:θu_g => 1.2,  :θr_g => 1.2),d0)
								   )

   					    sols,C,cp,pp = LandUse.runk(cpar = cpar,par = dd)
   					    p1 = plot(impl_plot_ts_all(C)..., layout = (2,2))
   					    p2 = impl_plot_slopes(C)

						plot(p1,p2,layout = (1,2))
	# 					p1
	end
end
