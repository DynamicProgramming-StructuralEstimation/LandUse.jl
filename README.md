# LandUse

[![Build Status](https://travis-ci.com/floswald/LandUse.jl.svg?token=yCXmyQ4r4F8RyxxzHZFG&branch=master)](https://travis-ci.com/floswald/LandUse.jl)

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
4. Download this package by pasting this into the `(v1.2) pkg>` prompt and hitting enter. (Only works for authorized users - i.e you! :smile:)
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
7. Done! :tada: Now try it out. Go back to command mode with `ctrl-c`
    ```julia
    julia> using LandUse

    julia> @time (x,M,p) = LandUse.run();
    8.649693 seconds (32.37 M allocations: 1.642 GiB, 8.71% gc time) 
    # first call needs to compile. slow.

    julia> @time (x,M,p) = LandUse.run();
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
9. Compare to matlab performance.  
    a. Matlab code is in `matlab/`
    ```
    starting values or each period:

    ans =
    
        0.0593    0.0098    0.0622    0.5545    1.0918    0.9239
        0.0587    0.0105    0.0621    0.5334    1.0764    0.9202
        0.0582    0.0112    0.0619    0.5137    1.0617    0.9168
        0.0571    0.0125    0.0618    0.4781    1.0341    0.9107
        0.0561    0.0137    0.0617    0.4467    1.0087    0.9053
        0.0547    0.0154    0.0617    0.4060    0.9743    0.8982
        0.0519    0.0185    0.0624    0.3333    0.9075    0.8855
        0.0460    0.0250    0.0680    0.2091    0.7706    0.8633
        0.0389    0.0332    0.0947    0.0968    0.5958    0.8410
        0.0367    0.0378    0.1406    0.0556    0.4986    0.8310
        0.0370    0.0397    0.1760    0.0436    0.4618    0.8275
        0.0395    0.0419    0.2532    0.0317    0.4180    0.8236
        0.0406    0.0424    0.2792    0.0295    0.4087    0.8228
        0.0418    0.0428    0.3052    0.0277    0.4009    0.8221
    
    period: 1
    starting epsilon search
    epsilon search step 2
    epsilon search step 3
    epsilon search step 4
    epsilon search step 5
    epsilon search step 6
    epsilon search step 7
    epsilon search step 8
    epsilon search step 9
    epsilon search step 10
    period: 2
    period: 3
    period: 4
    period: 5
    period: 6
    period: 7
    period: 8
    period: 9
    period: 10
    period: 11
    period: 12
    period: 13
    period: 14
    Elapsed time is 23.291344 seconds.
    ``` 
    
    b. Julia starting values and performace:
    ```
    julia> vcat(LandUse.get_starts()'...)
    14×6 Array{Float64,2}:
     0.0593117  0.00981358  0.0622341  0.554457   1.0918    0.923865
     0.0587322  0.0105377   0.0620709  0.533401   1.0764    0.920236
     0.05817    0.0112321   0.0619381  0.513739   1.06167   0.916847
     0.0570968  0.0125367   0.0617568  0.478103   1.03407   0.910699
     0.0560892  0.0137382   0.0616772  0.44669    1.00871   0.905272
     0.0546927  0.0153697   0.0617228  0.406043   0.974296  0.898232
     0.0518981  0.018535    0.0624406  0.333337   0.90746   0.885548
     0.0459825  0.0249714   0.068008   0.209129   0.770629  0.863256
     0.0389235  0.0331633   0.094704   0.0968216  0.595776  0.840951
     0.036694   0.0378419   0.140647   0.0556009  0.498566  0.830961
     0.0369503  0.0396774   0.17599    0.0436126  0.461788  0.827514
     0.0394578  0.0419417   0.253196   0.0316943  0.41795   0.823577
     0.0405955  0.0424347   0.279154   0.0294833  0.408668  0.822762
     0.0418213  0.0428495   0.305176   0.0277212  0.400933  0.822086

     julia> @time (x,M,p) = LandUse.run();
    0.067674 seconds (31.65 k allocations: 1.372 MiB)
    ```

