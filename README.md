# LandUse.jl

This is the julia code repository for our paper [Structural Change, Land Use and Urban Expansion](https://floswald.github.io/publication/landuse/). It contains all code needed to produce model outputs. Notice that there is a separate [R package (link)](https://github.com/floswald/LandUseR) which is used to perform data measurements for the same paper.

Online Documentation : [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://floswald.github.io/LandUse.jl)


Unit Tests : [![Build Status](https://github.com/floswald/LandUse.jl/workflows/CI/badge.svg)](https://github.com/floswald/LandUse.jl/actions?query=workflow%3ACI+branch%3Amaster)

## How to Use this

1. [Download julia](https://julialang.org/downloads/)
2. start julia. you see this:
    ```
                   _
       _       _ _(_)_     |  Documentation: https://docs.julialang.org
      (_)     | (_) (_)    |
       _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
      | | | | | | |/ _` |  |
      | | |_| | | | (_| |  |  Version 1.6.2 (2021-07-14)
     _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
    |__/                   |

    ```
3. Hit the `]` key to switch to package manager mode. the prompt switches to
    ```
    (v1.6) pkg>
    ```
4. Download this package by pasting this into the `(v1.6) pkg>` prompt and hitting enter. 
    ```julia
    dev https://github.com/floswald/LandUse.jl.git
    ```
5. After this is done, hit backspace or `ctrl-c` to go back to standard `julia>` prompt.
    ```julia
    julia> cd(joinpath(DEPOT_PATH[1],"dev","LandUse"))  # go to the location of LandUse
    ```
6. Go back to package mode: type `]`. then:
    ```julia
    (v1.6) pkg> activate .     # tell pkg manager to modify current directory as project
    (LandUse) pkg> instantiate    # download all dependencies
    ```
7. Done! :tada: Now try it out. Go back to command mode with `ctrl-c`. Run the standard model of a single region. Notice that the first call is slow because it has to precompile (in particular Ipopt):
    ```julia
    julia> using LandUse, Flux

    julia> @time (x,M,p) = LandUse.runm();
    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit https://github.com/coin-or/Ipopt
    ******************************************************************************
    
     48.071280 seconds (133.79 M allocations: 7.714 GiB, 4.13% gc time, 2.12% compilation time)
    
    julia> @time (x,M,p) = LandUse.runm();
      0.278983 seconds (1.07 M allocations: 166.982 MiB)
    
    julia> 

    # fast!
    ```
8. Run unit tests: (hit `]` to go back to pkg mode)
    ```julia
    (LandUse) pkg> test
     Testing Running tests...

    Test Summary: | Pass  Total
    LandUse.jl    |  326    326
         Testing LandUse tests passed 
    ```
9. Run the interactives. (Go back to command mode with `ctrl-c`)
    ```julia
    julia> LandUse.i0()  
    ```

![](docs/src/assets/single-interact.gif)
