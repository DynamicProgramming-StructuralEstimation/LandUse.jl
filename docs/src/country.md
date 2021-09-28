# Country

```@docs
LandUse.Country
```

Similarly to a single region, constructing a `Country` requires a parameter vector. We have to require at least 2 regions for this work. In this example we create a country with ``K=3`` regions, where the first region has 50% of all available space, and the the rest is split equally between cities 2 and 3 (key `:kshare`). The key `:factors` should be `ones(K)` by default, it allows an exogenous mulitplicative shift of productivity ``\theta_u``. Same for `:gs` which is a growth shift (by default zero). This is useful to experiment with different growth trajectories for regions according to which we would set the ``\theta_u`` of region `ik` to:

```julia
p.θu * p.factors[ik] * exp(p.gs[ik] * (p.it-1))
```

Here we create 3 regions with equal ``\theta_u``s:

```@example 2
using LandUse   # hide
using NNlib, Flux # hide
pk = LandUse.Param(par = Dict(:K => 3, 
                              :kshare => [0.5,0.25,0.25],
                              :factors => ones(3),
                              :gs => zeros(3)))
c3 = LandUse.Country(pk)      
c3
```   

## Solving for a `Country` in one period

We use the `Country` type mainly to keep country-wide variables in one place. The actual computation of the equilibrium is again performed in JuMP.jl in the [`jc`](@ref LandUse.jc) function:

```@docs
LandUse.jc
```

The main difference to the single city solution is that this imposes certain constraints on the chosen ``\theta_u`` sequence. Please consult the source code of [`jc`](@ref LandUse.jc) for more details.

## `Run`ning a Country

As before, we want to compute the country solution at a sequence of dates. This can be done with the [`runk`](@ref LandUse.runk). The default keywords run a 2-country setup:

```@example 2
x,C,p = LandUse.runk()  # sol vector, Country in each t, param
C
```

and as before, a certain index contains the country in that period:

```@example 2
C[6].R   # regions in period 6
```

