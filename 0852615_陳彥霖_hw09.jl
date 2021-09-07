import Statistics.mean
import Random.shuffle!

function setd(a::Vector, b::Vector)
    i = 1
    while i <= length(a)
        if a[i][1] < b[1] < a[i][2]
            push!(a, [a[i][1], b[1]])
            push!(a, [b[1], a[i][2]])
        elseif a[i][1] < b[2] < a[i][2]
            push!(a, [a[i][1], b[2]])
            push!(a, [ b[2], a[i][2]])
        end
        i += 1
    end
    filter!(x -> !(b[1] <= x[1] < x[2] <= b[2]), a)
    filter!(x -> !(x[1] < b[1] < x[2]), a)
    filter!(x -> !(x[1] < b[2] < x[2]), a)
    return a
end

function iunif(x::Vector)
    d = [a[2]-a[1] for a in x]
    v = rand(1)[1] * sum(d)
    i = findfirst(v .<= cumsum(d))
    return x[i][2] - cumsum(d)[i] + v
end

function p1(nsim = 1, k = 2e2, d = 0.1)
    output = Array{Array{Float64}}(undef, nsim)
    for i in 1:nsim
        y = shuffle!(collect(range(0, 1, length = 9)))
        for j in 1:k
            s = [[0., 1]]
            for r in setdiff(1:9, (j%9 + 1))
                s = setd(s, [y[r]-d, y[r]+d])
            end
            y[Int(j%9 + 1)] = iunif(s)
        end
        output[i] = y
    end
    return output
end

p1()

using Plots
histogram(vcat(p1(10000)...))


function p2a(nsim = 1, k = 2e2)
    output = Array{Array{Float64}}(undef, nsim)
    for i in 1:nsim
        x = [1.5, 2.5, 3.5]
        for j in 1:k
            if (j%3 + 1) == 1
                x[1] = -log(rand(1)[1]) + max(0, 15 - 2x[2] - 3x[3])
            elseif (j%3 + 1) == 2
                x[2] = -log(rand(1)[1]) + max(0, (15 - x[1] - 3x[3])/2)
            else
                x[3] = -log(rand(1)[1]) + max(0, (15 - x[1] - 2x[2])/3)
            end
        end
        output[i] = x
    end
    return output
end

mean([x[1]+2x[2]+3x[3] for x in p2a(100000)])

function p2b(nsim = 1, k = 2e2)
    output = Array{Array{Float64}}(undef, nsim)
    for i in 1:nsim
        x = [0.05, 0.1, 0.2]
        for j in 1:k
            if (j%3 + 1) == 1
                x[1] = -log(1-rand(1)[1]*(1-exp(-(1-2x[2]-3x[3]))))
            elseif (j%3 + 1) == 2
                x[2] = -log(1-rand(1)[1]*(1-exp(-((1-x[1]-3x[3])/2))))
            else
                x[3] = -log(1-rand(1)[1]*(1-exp(-((1-x[1]-2x[2])/3))))
            end
        end
        output[i] = x
    end
    return output
end

mean([x[1]+2x[2]+3x[3] for x in p2b(100000)])
