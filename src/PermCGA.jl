module PermCGA

using StatsBase
using Memoize

include("ocx.jl")

export permcga
export isvalid

struct PermutationProbs
    mat::Matrix{Float64}
    n::Int
end

abstract type AbstractMonitor end

struct DefaultMonitor <: AbstractMonitor
    probs::PermutationProbs
    maxiteration::Int
    currentiteration::Int
    bestcost::Float64
    bestperm::Vector
end


function PermutationProbs(mat::Matrix{Float64})
    n, _ = size(mat)
    return PermutationProbs(mat, n)
end

function defaultmonitorfunction(monitor::T) where {T<:AbstractMonitor}
    @info """
        Iteration $(monitor.currentiteration) of $(monitor.maxiteration), Cost: $(monitor.bestcost)
    """
end

initialscorematrix(n::Int)::PermutationProbs = PermutationProbs(ones(Float64, (n, n)))

function sampleperm(pm::PermutationProbs)::Vector{Int}
    n = pm.n
    cpm = copy(pm.mat)
    items = collect(1:n)
    result = Int[]
    lucky::Int64 = 0
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

function mutate!(sm::PermutationProbs, winner, loser, mutation)
    for p = 1:sm.n
        if winner[p] != loser[p]
            sm.mat[winner[p], p] += mutation
        else
            if sm.mat[loser[p], p] >= mutation
                sm.mat[loser[p], p] -= mutation 
            end
        end
    end
end

function flip(aperm::Vector{Int})::Vector{Int}
    result = copy(aperm)
    indices = StatsBase.sample(1:length(result), 2, replace = false)
    a, b = result[indices]
    c = result[a]
    result[a] = result[b]
    result[b] = c
    return result
end


function permcga(
    costfn::Function,
    n::Int,
    iterations::Int;
    monitorfunction::Union{Function,Nothing} = defaultmonitorfunction
)
   
    @memoize memoizedcostfn(x) = costfn(x)
    @warn "Memoized function will be used."
    

    sm            :: PermutationProbs = initialscorematrix(n)
    mutation      :: Float64 = 0.0
    fiterations    :: Float64 = Float64(iterations)
    floatn         :: Float64 = Float64(n)
    bestsolution  :: Vector{Int} = sampleperm(sm)
    bestcost      :: Float64 = memoizedcostfn(bestsolution)

    for i = 1:iterations

        mutation = floatn * (1.0 - i / fiterations)

        cand1 = sampleperm(sm) 
        cand2 = sampleperm(sm) 
        cost1 = memoizedcostfn(cand1)
        cost2 = memoizedcostfn(cand2)
        
        winner = cand1
        loser = cand2 
        currentmincost = cost1

        if (cost2 < cost1) 
            winner = cand2
            loser = cand1
            currentmincost = cost2
        end 
<

        mutate!(sm, winner, loser, mutation)
        

        if currentmincost < bestcost
            bestcost = currentmincost
            bestsolution = winner
            if (!isnothing(monitorfunction))
                monitorfunction(DefaultMonitor(sm, iterations, i, bestcost, bestsolution))
            end
        end

    end # end of iterations 


    resultdict = Dict(
        :solution => bestsolution,
        :cost => bestcost,
        :valid => isvalid(bestsolution),
        :scorematrix => sm.mat,
    )

    return resultdict
end


end # module PermCGA
