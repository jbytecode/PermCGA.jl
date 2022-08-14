module PermCGA

using StatsBase
using Memoize


export permcga

struct PermutationProbs
    mat::Matrix{Float64}
    n::Int
end

abstract type AbstractMonitor end 

struct DefaultMonitor <: AbstractMonitor
    probs               :: PermutationProbs 
    maxiteration        :: Int 
    steps               :: Int
    currentiteration    :: Int
    bestcost            :: Float64
    bestperm            :: Vector        
end 


function PermutationProbs(mat::Matrix{Float64})
    n, _ = size(mat)
    return PermutationProbs(mat, n)
end

function defaultmonitorfunction(monitor::T) where {T <: AbstractMonitor}
    if monitor.currentiteration % monitor.steps == 0
    @info """
        Iteration $(monitor.currentiteration) of $(monitor.maxiteration), Cost: $(monitor.bestcost)
    """
    end 
end

initialscorematrix(n::Int)::PermutationProbs = PermutationProbs(ones(Float64, (n, n)))

function sampleperm(pm::PermutationProbs)::Vector{Int}
    n = pm.n
    cpm = copy(pm.mat)
    items = collect(1:n)
    result = Int[]
    for i = 1:n
        lucky = sample(items, Weights(cpm[:, i]))
        push!(result, lucky)
        cpm[lucky, :] .= 0.0
    end
    return result
end

function extract(pm::PermutationProbs)::Vector{Int}
    map(x -> findmax(pm.mat[:, x])[2], 1:pm.n)
end

function isvalid(v::Vector{Int})::Bool
    return (v |> unique |> length) == length(v)
end

function permcga(
    costfn::Function, 
    n::Int, iterations::Int; 
    monitorfunction::Union{Function, Nothing} = defaultmonitorfunction)

    @memoize function memoizedcost(x)
        return costfn(x)
    end

    sm = initialscorematrix(n)
    mutation = 0.0
    fiterations = Float64(iterations)
    
    bestsolution = sampleperm(sm)
    bestcost = costfn(bestsolution)

    for i = 1:iterations

        mutation = n * (i / fiterations)

        cand1 = sampleperm(sm)
        cand2 = sampleperm(sm)
        cost1 = memoizedcost(cand1)
        cost2 = memoizedcost(cand2)
        winner = cand1
        loser = cand2
        currentmincost = cost1
        if cost2 < cost1
            winner = cand2
            loser = cand1
            currentmincost = cost2
        end

        for p = 1:n
            if winner[p] != loser[p]
                sm.mat[winner[p], p] += mutation
            end
        end

        if currentmincost < bestcost
            bestcost = currentmincost
            bestsolution = winner
        end

        if(!isnothing(monitorfunction))
            monitorfunction(DefaultMonitor(sm, iterations, 100, i, bestcost, bestsolution))
        end
        
    end

    # Sample section 
    for i = 1:iterations
        solution = sampleperm(sm)
        cost = costfn(solution)

        if isvalid(solution)
            if cost < bestcost
                bestcost = cost
                bestsolution = solution
                if(!isnothing(monitorfunction))
                     monitorfunction(DefaultMonitor(sm, iterations, 100, iterations, bestcost, bestsolution))
                end
            end
        end

    end

    resultdict = Dict(
        :solution => bestsolution,
        :cost => bestcost,
        :valid => isvalid(bestsolution),
        :scorematrix => sm.mat,
    )

    return resultdict
end


end # module PermCGA
