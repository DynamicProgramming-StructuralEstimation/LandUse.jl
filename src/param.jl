mutable struct Param
	γ     :: Float64 # housing weight
	ϵ     :: Float64 # housing supply elasticity
	η     :: Float64 # agglomeration forces
	ν     :: Float64 # weight of rural good consumption on consumption composite
	cbar  :: Float64 # agr cons subsistence level (0.14 for thailand)
	sbar  :: Float64 # neg urba subsistence level
	Θr    :: Vector{Float64} # set of rural sector TFP
	Θu    :: Vector{Float64} # set of urban sector TFP
	θr    :: Float64 # rural sector TFP
	θu    :: Float64 # urban sector TFP
	α     :: Float64 # labor weight on farm sector production function
	λ     :: Float64 # useless land: non-farm, non-urban land (forests, national parks...)
	τ     :: Float64 # commuting cost parameter
	χr    :: Float64 # "easiness" to convert land into housing in rural sector
	χu    :: Float64 # "easiness" to convert land into housing in urban sector
	L     :: Float64 # total population (note: total land normalized to 1)
	T     :: Vector{Float64}
	σ     :: Float64 # land-labor elasticity of substitution in farm production function
	c0    :: Float64  # construction cost function intercept
	c1    :: Float64  # construction cost function gradient
	c2    :: Float64  # construction cost function quadratic term
	Ψ     :: Float64  # urban ammenities rel to rural 
	int_nodes :: Int  # number of integration nodes

	function Param(;par=Dict())
        f = open(joinpath(dirname(@__FILE__),"params.json"))
        j = JSON.parse(f)
        close(f)
        this = new()
        for (k,v) in j
            if v["value"] isa Vector{Any}
                setfield!(this,Symbol(k),convert(Vector{Float64},v["value"]))
            else
                setfield!(this,Symbol(k),v["value"])
            end
        end

        if length(par) > 0
            # override parameters from dict par
            for (k,v) in par
                setfield!(this,k,v)
            end
        end
    	return this
	end
end