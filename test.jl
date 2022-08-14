using PermCGA 
using Plots

function createpoints()
    pts = Array{Float64, 2}(undef, (0, 2))
    for i in 1:10
       pts = vcat(pts, [i, 10]')
    end 
    for i in 9:(-1):1
       pts = vcat(pts, [10, i]')
    end
    for i in 9:(-1):1
       pts = vcat(pts, [i, 1]')
    end
    for i in 2:9
       pts = vcat(pts, [1, i]')
    end
    return pts 
end 

euclidean(u, v) = (u .- v) .|> a -> a * a |> sum |> sqrt 

pts = createpoints()
n, p = size(pts)
distmat = Array{Float64, 2}(undef, (n , n ))
for i in 1:n 
    for j in 1:n
        distmat[i, j] = euclidean(pts[i], pts[j])
    end
end 

function costfn(permutation)
    n = length(permutation)
    totaldist = 0.0
    for i in 1:(n-1)
        totaldist += distmat[permutation[i], permutation[i + 1]]
    end 
    totaldist += distmat[permutation[n], permutation[1]]
    return totaldist 
end



result = permcga(costfn, n, 500000)
sol = result[:solution]
push!(sol, sol[1])
Plots.scatter(pts[:,1], pts[:,2])
x = pts[sol, 1]
y = pts[sol, 2]
Plots.plot!(x, y)

