

# 1) plot methods for a single model instance
# visualize spatial setup
# moving cost function, construction cost function, density, cu, cr,h consumption, 

function plot_static(m::Model,p::Param)
	lvec = collect(range(0.,1.0,length=100))
	plts = Any[]
	push!(plts,plot(lvec,[τ(i,m.ϕ,p) for i in lvec], title = "tau",leg = false))
	push!(plts,plot(lvec,[w(m.Lu,i,m.ϕ,p) for i in lvec], title = "w(l)",leg = false))
	push!(plts,plot(lvec,[cu(i,p,m) for i in lvec], title = "cu(l)",leg = false))
	push!(plts,plot(lvec,[cr(i,p,m) for i in lvec], title = "cr(l)",leg = false))
	push!(plts,plot(lvec,[cost(i,m.ϕ,p) for i in lvec], title = "cost(l)",leg = false))
	plot(plts...)
end





# 2) plotting a time series of models
# for each model instance at time t, plot against t:
# fringe
# qr
# Lr, Lu
# land rent
# density in center and at fringe
# 
