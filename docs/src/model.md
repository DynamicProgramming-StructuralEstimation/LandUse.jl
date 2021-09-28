# Single City Model

This page describes the how to use the single city model, what we call a `Region`.

## Parameter and `Region` Setup

The first step is to obtain a [`Param`](@ref LandUse.Param) type, like so:

```@example 1
using LandUse
using NNlib, Flux
p = LandUse.Param()
p
```

and then we create a [`Region`](@ref LandUse.Region) :

```@example 1
m = LandUse.Region(p)
typeof(m)
```

*Running a model* refers to computing the solution to the equation system defined in the [`jm`](@ref  LandUse.jm) function via JuMP.jl at all periods:

```@docs
LandUse.jm
```

We compute a starting value ``x_0``, and then supply the period ``t`` solution vector as starting value to the period ``t+1`` problem. We obtain this model run via

```@example 1
x,M,p = LandUse.run(p)  # solution vectors, Region in each period, Param
M
```

And you can just index the `M` vector to look at a certain period in detail:

```@example 1
M[4]
```

## Change Parameter Values

The constructor to the `Param` type accepts a keyword `par` which has to be dict. You can change the default values by using the appropriate `key => value` pairs:

```@example 1
p2 = LandUse.Param(par = Dict(:α => 0.5))
p2.α
```

