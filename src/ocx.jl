# Source:
# https://www.redalyc.org/pdf/2652/265219618002.pdf
function ocx(p1, p2; c1=0, c2=0)::Tuple{Array{Int64, 1}, Array{Int64, 1}}
    n = length(p1)
    #c1 = 3
    #c2 = 6
    if c1 == 0 || c2 == 0
        c1, c2 = sort(sample(2:(n-1), 2, replace = false))
    end

    part11 = p1[1:(c1-1)]
    part12 = p1[c1:(c2-1)]
    part13 = p1[c2:end]

    part21 = p2[1:(c1-1)]
    part22 = p2[c1:(c2-1)]
    part23 = p2[c2:end]

    o1 = vcat(part13, part11, part12)
    o1 = filter(x -> !(x in part22), o1)
    o1 = vcat(o1[end-length(part11)+1:end], part22, o1[1:length(part13)])

    o2 = vcat(part23, part21, part22)
    o2 = filter(x -> !(x in part12), o2)
    o2 = vcat(o2[end-length(part21)+1:end], part12, o2[1:length(part23)])

    return (o1, o2)
end