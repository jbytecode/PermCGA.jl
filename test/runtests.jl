using PermCGA
using Test

@testset "isvalid" begin
    @test PermCGA.isvalid([1,2,3,4,5])
    @test PermCGA.isvalid([1,2,3,4,5,6,7,8,9,10])
    @test PermCGA.isvalid([10,9,8,7,6,1,2,3,4,5])
    @test PermCGA.isvalid([1,2,3,4,5,5]) == false
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
