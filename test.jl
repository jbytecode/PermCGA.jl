using PermCGA 

thesolution = collect(1:5)

    function cost(x)::Float64
        result = thesolution
        costval = findall(a -> a != 0, abs.(result .- x)) |> length
        return Float64(costval)
    end

    result = permcga(cost, length(thesolution), 10000)
