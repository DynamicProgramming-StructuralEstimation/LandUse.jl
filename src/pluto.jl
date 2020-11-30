### A Pluto.jl notebook ###
# v0.12.15

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ ef8d966c-332d-11eb-3269-a32386247453
using LandUse

# ╔═╡ 374a9ed6-332f-11eb-2410-cbb208137934
using PlutoUI

# ╔═╡ 10e03f36-332e-11eb-182b-bfe90c016f06
LandUse.plot1()

# ╔═╡ 851e0720-332e-11eb-22d9-cb94753cdacd
@bind ηm Slider(0.5:0.01:2.0)

# ╔═╡ 7dd0fb98-332f-11eb-2a06-ff44e1657ee9
"""md eta is $(ηm)"""

# ╔═╡ 47e18142-332f-11eb-292d-25dd53a0418a
p0 = LandUse.Param(par = Dict(:ϵsmax => 0.0));

# ╔═╡ 5d02b118-332f-11eb-1767-71a66abbb820
begin
	x,M,p = LandUse.run(p0)
						pl= LandUse.ts_plots(M,p0,fixy = false)
						LandUse.plot(pl[:Lr_data],pl[:spending],pl[:pr_data],pl[:productivity],
						     pl[:n_densities], pl[:avdensity], pl[:mode], pl[:ctime],
							 pl[:phi] , pl[:qbar_real], pl[:r_y], pl[:r_rho],
							 layout = (4,3),link = :x,size = (700,800))
end

# ╔═╡ a002c484-3330-11eb-0c1a-130c48ad63bc


# ╔═╡ Cell order:
# ╠═ef8d966c-332d-11eb-3269-a32386247453
# ╠═10e03f36-332e-11eb-182b-bfe90c016f06
# ╠═374a9ed6-332f-11eb-2410-cbb208137934
# ╠═851e0720-332e-11eb-22d9-cb94753cdacd
# ╠═7dd0fb98-332f-11eb-2a06-ff44e1657ee9
# ╠═47e18142-332f-11eb-292d-25dd53a0418a
# ╠═5d02b118-332f-11eb-1767-71a66abbb820
# ╠═a002c484-3330-11eb-0c1a-130c48ad63bc
