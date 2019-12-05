# LandUse

[![Build Status](https://travis-ci.com/floswald/LandUse.jl.svg?token=yCXmyQ4r4F8RyxxzHZFG&branch=master)](https://travis-ci.com/floswald/LandUse.jl)

## How to Use this

1. [Download julia](https://julialang.org/downloads/)
2. start julia. you see this:
    ```
    ➜  julia
                   _
       _       _ _(_)_     |  Documentation: https://docs.julialang.org
      (_)     | (_) (_)    |
       _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
      | | | | | | |/ _` |  |
      | | |_| | | | (_| |  |  Version 1.1.0 (2019-01-21)
     _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
    |__/                   |
    

    julia> 
    ```
3. Hit the `]` key to switch to package manager mode. the prompt switches to 
    ```
    (v1.1) pkg>
    ```
4. Download this package by pasting this into the `(v1.1) pkg>` prompt and hitting enter. (Only works for authorized users - i.e you! :smile:)
    ```julia
    dev https://github.com/floswald/LandUse.jl.git
    ```
5. After this is done, hit backspace or `ctrl-c` to go back to standard `julia>` prompt.
    ```julia
    julia> cd(joinpath(DEPOT_PATH[1],"dev","LandUse"))  # go to the location of LandUse
    ```
6. Go back to package mode: type `]`. then:
    ```julia
    (v1.1) pkg> activate .     # tell pkg manager to modify current directory as project
    (LandUse) pkg> instantiate    # download all dependencies
    ```
7. Done! :tada: Now try it out. Go back to command mode with `ctrl-c`
    ```julia
    julia> using LandUse

    julia> @time m = LandUse.main(); # first call needs to compile. slow.
        Results of Nonlinear Solver Algorithm
     * Algorithm: Trust-region with dogleg and autoscaling
     * Starting Point: [0.10726751457398807, 0.13812932090095717, 0.10726751457398807, 0.6301030070492295, 1.3352070297616616, 0.6265735501145178]
     * Zero: [0.07690277582937129, 0.02749276302115747, 0.07733042441525058, 0.665630178180822, 1.2160416979906072, 0.9077109356334955]
     * Inf-norm of residuals: 0.000000
     * Iterations: 6
     * Convergence: true
       * |x - x'| < 0.0e+00: false
       * |f(x)| < 1.0e-08: true
     * Function Calls (f): 7
     * Jacobian Calls (df/dx): 7
    qr = 0.07690277582937129
    ϕ = 0.02749276302115747
    r = 0.07733042441525058
    Lr = 0.665630178180822
    pr = 1.2160416979906072
    Sr = 0.9077109356334955
    12.316737 seconds (30.61 M allocations: 1.556 GiB, 6.89% gc time)

    julia> @time m = LandUse.main();  # much faster!
    Results of Nonlinear Solver Algorithm
     * Algorithm: Trust-region with dogleg and autoscaling
     * Starting Point: [0.10726751457398807, 0.13812932090095717, 0.10726751457398807, 0.6301030070492295, 1.3352070297616616, 0.6265735501145178]
     * Zero: [0.07690277582937129, 0.02749276302115747, 0.07733042441525058, 0.665630178180822, 1.2160416979906072, 0.9077109356334955]
     * Inf-norm of residuals: 0.000000
     * Iterations: 6
     * Convergence: true
       * |x - x'| < 0.0e+00: false
       * |f(x)| < 1.0e-08: true
     * Function Calls (f): 7
     * Jacobian Calls (df/dx): 7
    qr = 0.07690277582937129
    ϕ = 0.02749276302115747
    r = 0.07733042441525058
    Lr = 0.665630178180822
    pr = 1.2160416979906072
    Sr = 0.9077109356334955
      0.005716 seconds (2.96 k allocations: 555.266 KiB)
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

