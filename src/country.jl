

mutable struct Country
	R  :: Vector{Region}   # a set of regions
	pr :: Float64          # a global price of rural good
	ρr :: Float64          # a global land price
	r  :: Float64          # country-wide per capita rental income

	function Country(p::Param)
		this = new()
		this.R = [Region(p) for k in 1:p.K]
		this.pr = NaN
		this.ρr = NaN
		this.r  = NaN
		return this
	end

end