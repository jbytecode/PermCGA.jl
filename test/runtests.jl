using PermCGA
using Test


euclidean(u, v) = (u .- v) |> (a -> a .* a) |> sum |> sqrt

@testset "Euclidean distance" begin
    @test euclidean([0.0, 0.0], [0.0, 5.0]) == 5.0
end

@testset "isvalid" begin
    @test PermCGA.isvalid([1, 2, 3, 4, 5])
    @test PermCGA.isvalid([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    @test PermCGA.isvalid([10, 9, 8, 7, 6, 1, 2, 3, 4, 5])
    @test PermCGA.isvalid([1, 2, 3, 4, 5, 5]) == false
end

@testset "ocx" begin
    p1 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    p2 = [8, 5, 7, 1, 2, 4, 9, 3, 6]
    
    result = PermCGA.ocx(p1, p2, c1 = 3, c2 = 6)
    
    @test length(result) == 2
    @test result[1] == [4, 5, 7, 1, 2, 6, 8, 9, 3]
    @test result[2] == [1, 2, 3, 4, 5, 9, 6, 8, 7]
end

@testset "Simple test with n = 5" begin

    thesolution = [1, 2, 3, 4, 5]

    function cost(x)::Float64
        result = thesolution
        costval = findall(a -> a != 0, abs.(result .- x)) |> length
        return Float64(costval)
    end

    result = permcga(cost, 5, 500)

    #result[:scorematrix] |> display

    @test result[:valid] == true
    @test result[:solution] == thesolution
    @test result[:cost] == 0.0

end

@testset "Simple test with n = 10" begin

    thesolution = [1, 3, 5, 7, 9, 2, 4, 6, 8, 10]

    function cost(x)::Float64
        costval = findall(a -> a != 0, abs.(thesolution .- x)) |> length
        return Float64(costval)
    end

    result = permcga(cost, 10, 2000)

    #result[:scorematrix] |> display

    @test result[:valid] == true
    @test result[:solution] == thesolution
    @test result[:cost] == 0.0

end



@testset "Traveling Salesman with 5x5" begin

    datapoints = [
        1 1
        10 1
        10 20
        1 20
    ]




    distances = zeros(4, 4)
    for i = 1:4
        for j = 1:4
            distances[i, j] = euclidean(datapoints[i, :], datapoints[j, :])
        end
    end

    function cost(x)
        y = copy(x)
        push!(y, first(x))
        mysum = 0.0
        for i = 1:3
            mysum = mysum + distances[y[i], y[i+1]]
        end
        return mysum
    end

    result = permcga(cost, 4, 500)


    @test result[:cost] == 37.0
    @test result[:valid]

    result |> display

end



@testset "TSP with 16 nodes" begin

    function createpoints()
        pts = Array{Float64,2}(undef, (0, 2))
        for i = 1:5
            pts = vcat(pts, [i, 5]')
        end
        for i = 4:(-1):1
            pts = vcat(pts, [5, i]')
        end
        for i = 4:(-1):1
            pts = vcat(pts, [i, 1]')
        end
        for i = 2:4
            pts = vcat(pts, [1, i]')
        end
        return pts
    end



    pts = createpoints()
    n, p = size(pts)
    distmat = Array{Float64,2}(undef, (n, n))
    for i = 1:n
        for j = 1:n
            distmat[i, j] = euclidean(pts[i, :], pts[j, :])
        end
    end

    function costfn(permutation)
        n = length(permutation)
        totaldist = 0.0
        for i = 1:(n-1)
            totaldist += distmat[permutation[i], permutation[i+1]]
        end
        totaldist += distmat[permutation[n], permutation[1]]
        return totaldist
    end

    result = permcga(costfn, n, 50000)

    @test result isa Dict
    @test result[:cost] == 16.0
    @test result[:solution] |> length == n
end
