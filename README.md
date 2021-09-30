# LandUse

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://floswald.github.io/LandUse.jl)


GitHub Actions : [![Build Status](https://github.com/floswald/LandUse.jl/workflows/CI/badge.svg)](https://github.com/floswald/LandUse.jl/actions?query=workflow%3ACI+branch%3Amaster)

## How to Use this

1. [Download julia](https://julialang.org/downloads/)
2. start julia. you see this:
    ```
    ➜  LandUse git:(eps) ✗ julia
               _
       _       _ _(_)_     |  Documentation: https://docs.julialang.org
      (_)     | (_) (_)    |
       _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
      | | | | | | |/ _` |  |
      | | |_| | | | (_| |  |  Version 1.2.0 (2019-08-20)
     _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
    |__/                   |

    ```
3. Hit the `]` key to switch to package manager mode. the prompt switches to
    ```
    (v1.2) pkg>
    ```
4. Download this package by pasting this into the `(v1.2) pkg>` prompt and hitting enter. 
    ```julia
    dev https://github.com/floswald/LandUse.jl.git
    ```
5. After this is done, hit backspace or `ctrl-c` to go back to standard `julia>` prompt.
    ```julia
    julia> cd(joinpath(DEPOT_PATH[1],"dev","LandUse"))  # go to the location of LandUse
    ```
6. Go back to package mode: type `]`. then:
    ```julia
    (v1.2) pkg> activate .     # tell pkg manager to modify current directory as project
    (LandUse) pkg> instantiate    # download all dependencies
    ```
7. Done! :tada: Now try it out. Go back to command mode with `ctrl-c`. Run the standard model of a single region:
    ```julia
    julia> using LandUse

    julia> @time (x,M,p) = LandUse.run(LandUse.Region,LandUse.Param());
    8.649693 seconds (32.37 M allocations: 1.642 GiB, 8.71% gc time)
    # first call needs to compile. slow.

    julia> @time (x,M,p) = LandUse.run(LandUse.Region,LandUse.Param())();
    0.061881 seconds (31.64 k allocations: 1.369 MiB)

    # fast!
    ```
8. Run unit tests: (hit `]` to go back to pkg mode)
    ```julia
    (LandUse) pkg> test
       Testing LandUse
     Resolving package versions...
    Test Summary: | Pass  Total
    LandUse.jl    |   44     44
       Testing LandUse tests passed
    ```
9. Run the interactives. (Go back to command mode with `ctrl-c`)
    ```julia
    # hit ? to go to help mode:
    help?> LandUse.i22
       https://github.com/floswald/LandUse.jl/issues/22

    # run it
    julia> LandUse.i22()  
    ```
