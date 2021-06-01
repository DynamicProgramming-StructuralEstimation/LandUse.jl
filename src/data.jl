

function relpop(; maxrank = 20)
    d = CSV.read(joinpath(LandUse.dboutdata, "relpop.csv"), DataFrame, types = Dict(2 => String))
    d20 = subset(d, :rank => x -> x .< maxrank)
    d1 = sort(subset(d20, :year => x-> x .== 1876), :relpop, rev = true ) # fist year

    # classify cities
    d1[!, :group] .= 0
    d1[d1.relpop .== 1.0, :group] .= 1
    d1[d1.relpop .<  1.0, :group]  .= 2
    d1[d1.relpop .<  0.16, :group] .= 3
    d1[d1.relpop .<  0.08, :group] .= 4
    d1[d1.relpop .<  0.043, :group] .= 5
    d1 = transform(groupby(d1, :group), nrow => :ngroup)

    # merge back into data
    d20 = leftjoin(d20,select(d1, :CODGEO, :group, :ngroup), on = :CODGEO)

    # create group means
    dm = combine(groupby(d20, [:group, :year]), :ngroup => first => :ngroup, :relpop => mean => :relpop)

    # plot
    @df subset(dm, :group => x -> x .> 1) plot(:year, :relpop, group = :group)
    @df subset(d20, :group => x -> x .> 1) plot(:year, :relpop, group = :LIBGEO)

    # manually smooth out marseille in 1936: set equal to Lyon
    d20[(d20.LIBGEO .== "Marseille") .& (d20.year .== 1936), :relpop] .= d20[(d20.LIBGEO .== "Lyon") .& (d20.year .== 1936), :relpop]

     # create group means
     dm = combine(groupby(d20, [:group, :year]), :ngroup => first => :ngroup, :relpop => mean => :relpop)

     # plot
     @df subset(dm, :group => x -> x .> 1) plot(:year, :relpop, group = :group, marker = :circle)



end