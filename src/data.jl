

function popdata(; maxrank = 20)
    d = CSV.read(joinpath(LandUse.dboutdata, "relpop.csv"), DataFrame, types = Dict(2 => String))
    d20 = subset(d, :rank => x -> x .<= maxrank)
    d20 = transform(groupby(d20, :year), :pop_data => (x -> x / sum(x)) => :popweight)

    d1 = sort(subset(d20, :year => x-> x .== 1876), :relpop_data, rev = true ) # fist year

    # classify cities
    d1[!, :group] .= 0
    d1[d1.relpop_data .== 1.0, :group] .= 1
    d1[d1.relpop_data .<  1.0, :group]  .= 2
    d1[d1.relpop_data .<  0.13, :group] .= 3
    d1[d1.relpop_data .<  0.043, :group] .= 4
    d1[d1.relpop_data .<  0.029, :group] .= 5
    d1 = transform(groupby(d1, :group), nrow => :ngroup)
    CSV.write(joinpath(LandUse.dboutdata, "relpop-classification.csv"), d1)

    # merge back into data
    d20 = leftjoin(d20,select(d1, :CODGEO, :group, :ngroup), on = :CODGEO)

    # create group means
    dm = combine(groupby(d20, [:group, :year]), :ngroup => first => :ngroup, :relpop_data => mean => :relpop_data)

    # plot
    @df subset(dm, :group => x -> x .> 1) plot(:year, :relpop_data, group = :group)
    @df subset(d20, :group => x -> x .> 1) plot(:year, :relpop_data, group = :LIBGEO)

    # manually smooth out marseille in 1936: set equal to Lyon
    # d20[(d20.LIBGEO .== "Marseille") .& (d20.year .== 1936), :relpop] .= d20[(d20.LIBGEO .== "Lyon") .& (d20.year .== 1936), :relpop]

    # create group means
    dm = combine(groupby(d20, [:group, :year]), :ngroup => first => :ngroup, :relpop_data => mean => :relpop_data, :popweight => sum => :popweight)
    transform!(dm, [:group, :ngroup] => ((a,b) -> string.(a,"(n=",b,")")) => :grouplabel)

    # plots
    pl = @df subset(dm, :group => x -> x .> 1) plot(:year, :relpop_data, group = :grouplabel, marker = :circle, title = "Population Data relative to Paris")
    pl2 = @df subset(d20, :group => x -> x .> 1) plot(:year, :relpop_data, group = :LIBGEO, marker = :circle, title = "Population Data relative to Paris")

    CSV.write(joinpath(LandUse.dboutdata, "relpop-means.csv"), dm)
    CSV.write(joinpath(LandUse.dboutdata, "relpop-full.csv"), d20)

    return (pl, pl2, d20, dm)
end

"""
find closest year in population data to model years
"""
function popdata_mapyears(p::Param)
    d = CSV.read(joinpath(LandUse.dboutdata, "relpop-full.csv"), DataFrame)
    # find closest year in data
    datayears = unique(d.year)
    modelyears = p.T
    df = DataFrame(modelyears = modelyears, datayears = zeros(Int,length(modelyears)))
    for ir in eachrow(df)
        ir[:datayears] = datayears[argmin( abs.(datayears .- ir[:modelyears]) )]
    end
    df
end

"""
load population and area data
"""
function poparea_data()
    CSV.read(joinpath(LandUse.dboutdata, "france_final.csv"), DataFrame)
end