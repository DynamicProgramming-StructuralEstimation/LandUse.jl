# Single City Model

This page describes the how to use the single city model. 

## Parameter Setup

The first step is to obtain a [`Param`](@ref) type, like so:

```jldoctest
using LandUse, Flux  # you need to load Flux as well
p = LandUse.Param()
p.L

# output
1.0
```



