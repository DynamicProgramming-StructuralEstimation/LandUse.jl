

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
	t1s = 0.1:0.05:0.99
	zetas = 0.1:0.05:0.99


	p = Param()
	@manipulate for θtype in Dict("Matlab θ" => 1, "Growth θ" => 2),
		            sb in slider(sbs, label = "sbar", value = p.sbar),
		            cb in slider(cbs, label = "cbar", value = p.cbar),
		            sig in slider(sis, label = "σ", value = p.σ),
		            # agg in slider(ags, label = "agg_g", value = p.θagg_g),
		            θug in slider(θus, label = "θu_g", value = p.θu_g),
		            esl in slider(esl, label = "ϵ-slope", value = 0.0),
		            epl in slider(eps, label = "ϵ", value = p.ϵr),
		            θrg in slider(θrs, label = "θr_g", value = p.θr_g),
		            t1 in slider(t1s, label = "τ1", value = p.τ1),
		            zeta in slider(zetas, label = "ζ", value = p.ζ)


		if θtype == 1
			x,M,p = run(Region,
			             Param(par = Dict(:ζ => zeta,:τ1 => t1,:ϵs => 0.0, :ϵsmax => esl,:ϵr => epl,
						                  :sbar => sb, :cbar => cb, :σ => sig)))
		else
	    	x,M,p = run(Region,
		             Param(par = Dict(:ζ => zeta,:τ1 => t1,:ϵs => 0.0, :ϵsmax => esl,:ϵr => epl,
					                  :sbar => sb, :cbar => cb, :σ => sig,
									  :θut => 0.32, :θrt => 0.32, :θu_g => θug, :θr_g => θrg)))
        end
		@info "model done."
		println()

	    d = dataframe(M,p.T)
	    df = @linq d |>
	         select(:year,:Ch,:Cu,:Cr,:C ) |>
	         transform(h = :Ch ./ :C,u = :Cu ./ :C,r = :Cr ./ :C) |>
	         select(:year, :h, :u , :r)
	    ds = stack(df, Not(:year))
	    pl = @df ds plot(:year,:value, group = :variable,
	               linewidth = 2, title = "Spending Shares",
				   ylims = (0.0,0.85))

		ds2 = stack(select(d,:year,:Lu, :Lr), Not(:year))
		pl2 = @df ds2 plot(:year, :value, group = :variable,
		                 linewidth = 2, title = "Population",
						 ylims = (0,1))
	    d3  = @linq d |>
		        select(:year,:θu, :θr) |>
				transform(theta_u = :θu, theta_r = :θr) |>
				select(:year,:theta_u, :theta_r)

		# ds3 = stack(d3, Not(:year))
	    # pl3 = @df ds3 plot(:year, :value, group = :variable,
		# 			      linewidth = 2, title = "Consumption",
		# 				  ylims = (0.0,3.5))

	    ds3 = stack(d3, Not(:year))
	    pl3 = @df ds3 plot(:year, :value, group = :variable,
					      linewidth = 2, title = "Productivity",
						  ylims = (0.0,3.5))
	    ds4 = stack(select(d,:year,:ϕ), Not(:year))
		pl4 = @df ds4 plot(:year, :value, group = :variable,
		                 linewidth = 2, title = "city size",
						 ylims = (0.0,0.3),
						 leg = false)
	    # ds4 = stack(select(d,:year,:pr), Not(:year))
		df4 = @linq d |>
			# transform(avgd = :Lu ./ :ϕ, tauphi = (1 .- map(x -> τ(x,x,p),:ϕ) .* :ϕ) ./ :pr)
			transform(avgd = :Lu ./ :ϕ, tauphi = (1 .- map(x -> τ(x,x,p),:ϕ) .* :ϕ) )
			# transform(avgd = :Lu ./ :ϕ, tauphi = (1 .- p.τ .* :ϕ))

	    ds4 = stack(select(df4,:year,:d0,:dq1, :dq2, :dq3, :dq4, :dr, :avgd), Not(:year))
	    ds5 = stack(select(df4,:year,:avgd), Not(:year))
	    ds6 = stack(select(df4,:year,:tauphi), Not(:year))
	    # ds4 = stack(select(df4,:year, :avgd), Not(:year))
		pl5 = @df ds4 plot(:year, :value, group = :variable,
						 linewidth = 2, title = "densities")
		pl6 = @df ds5 plot(:year, :value, group = :variable,
		                 linewidth = 2, title = "densities")
		 pl7 = @df ds6 plot(:year, :value, group = :variable,
						  linewidth = 2, title = "1 - tau(phi)",leg = false)
		# plot(pl,pl2,pl3,pl4,pl5, layout = (1,5),size = (700,250))
		# plot(pl2,pl6,pl4,pl7 ,layout = (1,4),size = (900,250))
		plot(pl2,pl6,pl4,pl7 ,layout = (1,4),size = (900,250))
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
	esl = 1.0:20:400.0
	es = 1.0:0.5:10.0
	esm = 1.0:0.5:10.0
	zetas = 0.0:0.1:3.3
	mp = @manipulate for e in slider(es, label = "ϵr", value = 4.0 ),
		                 esl in slider(esl,label = "eslope"),
		                 gr in slider(g1,label = "gr"),
		                 gs1 in slider(g1,label = "ug1"),
		                 gs2 in slider(g1,label = "ug2"),
		                 gs3 in slider(g1,label = "ug3"),
		                 gagg in slider(g1,label = "agg"),
						 z in slider(zetas,label = "ζ", value = 0.0),
						 type in Dict("TS" => 1, "phi vs Lu" => 2)

						 x = issue12(e;gr = gr, gu = [gs1,gs2,gs3], θu = [0.11,0.32,0.35],θagg_g = gagg, ϵsmax = esl, zeta = z)
						 # plot(rand(10))
						 if type == 1
						 	x[5]
						else
							x[6]
						end
					 end

end
