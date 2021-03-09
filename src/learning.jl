function make_grid_df(p0::Param,n)
    so = search_over()
    arrs = []
    for (k,v) in so
        gup   = range(getfield(p0,k), stop = v[2], length = div(n, 2))  #
        gdown = range(getfield(p0,k), stop = v[1], length = div(n, 2))[2:end]
        push!(arrs, [gup; gdown])
    end
    d = DataFrame(gridmake(arrs...), collect(keys(so)) )
    d.id = 1:nrow(d)
    d

end

function grid_arr(p::Param,n)
    x = make_grid_df(p,n)
    Array(x[!,1:(ncol(x)-1)])
end

function eval_grid(g::Matrix)
    n = size(g)[1]
    r = zeros(n)
    @showprogress for i in 1:n
        r[i] = objective(g[i,:])
    end
    r
end

function worker_receiver(ids::Vector{Int},n::Int)
    # ids are indices to work on for this worker
    jobs = grid_arr(Param(),n)
    out = zeros(length(ids))
    for ix in 1:length(ids)
        out[ix] = objective(jobs[ix,:])
    end
    return out
end




function peval_grid(g::SharedArray)
    n = size(g)[1]
    out = @showprogress pmap(1:n) do x
        objective(g[x,:])
    end
    out
end

# function parallel_setup(nworkers)

#     # create workers
#     w = addprocs(nworkers, exeflags = "--project=.")

#     # load code on all workers
#     @everywhere using LandUse

#     return w

# end

function parallel_grid(npoints,nworkers)
    
    # create shared Array
    sg = SharedArray{Float64,2}(grid_arr(Param(),npoints))

    # run in parallel
    r = peval_grid(sg)

    # add to grid, save and return
    sg = hcat(sg,r)
    writedlm(joinpath(@__DIR__,"..","out","par_grid.txt"), sg)
    sg

end

# also need to permutedims and exchange which is the first column!
# only first col will work, because we start from know solution and extrapolte
function shared_grid(x::DataFrame) 
    a = Array(x[!,1:(ncol(x)-1)])
    a = hcat(a, fill(NaN,size(a)[1]))
    SharedArray{Float64,2}(a)
end

# get all reachable starting points from x0
# given our current starting value in startval(p) - which is constant - 
# compute for a grid of params all solutions that converge and record their first period solution
# from that dataset - param input and implied solution

"""
checks whether a certain row of the grid together with a starting value x0 
produces a valid solution in all periods.
"""
function isfeasible_grid_row(r::Array,x0::NamedTuple)
    p1 = Param(par = x2dict(r[1:(end-1)]))
    try
        x1,M1,p2 = run(x0,p1)
        return 1
    catch
        return 0
    end

end



function gridsearchPOC(n)

    # 1. starting value we know works is
    p0 = Param()
    x0,M0,p0 = run(p0)   # x0[1] is first period of a working solution 

    solsup = []
    solsdown = []

    # 2. make grids
    so = search_over()
    for (k,v) in so
        gup   = range(getfield(p0,k), stop = v[2], length = n)
        gdown = range(getfield(p0,k), stop = v[1], length = n)
        # go up 
        push!(solsup, x0[1])
        for (ip,p) in enumerate(gup[2:end])
            println(p)
            p1 = Param(par = Dict(k => p))
            x1,M1,p2 = run(solsup[ip],p1)
            push!(solsup, x1[1])
        end
        push!(solsdown, x0[1])

        for (ip,p) in enumerate(gdown[2:end])
            println(p)
            p1 = Param(par = Dict(k => p))
            x1,M1,p2 = run(solsdown[ip],p1)
            push!(solsdown, x1[1])
        end
    end
    [solsup ; solsdown]
end
