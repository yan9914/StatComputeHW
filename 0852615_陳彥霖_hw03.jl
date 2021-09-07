function Box_Muller(n::Int64)
    k = iseven(n) ? Int64(n/2) : Int64((n+1)/2)
    output = Array{Float64}(undef, 2k)
    r = -2log.(rand(k))
    θ = rand(k)*.2π
    @inbounds @simd for i in 1:k
        output[2i-1] = sqrt(r[i])*cos(θ[i])
        output[2i] = sqrt(r[i])*sin(θ[i])
    end
    return iseven(n) ? output : output[1:length(output)-1]
end

x = Box_Muller(10000)

using Gadfly, DataFrames, Distributions

plot(DataFrame(x = x),
    layer(x -> pdf(Normal(0, 1), x), min(x...), max(x...), Geom.line, Theme(default_color="red")),
    layer(x = :x, Geom.histogram(bincount = 100, density = true))
)

using Optim
function NHPP(FUN::Function, start::Float64, endT::Float64)
    t = start
    λ = -minimum(optimize(x -> -FUN(x), start, endT))
    S = Float64[]
    while true
        tmp = t - 1/λ * log(rand())
        if rand() < FUN(tmp)/λ
            t = tmp
            if t > endT
                return S
            else
                push!(S, t)
            end
        end
    end
end

using BenchmarkTools
NHPP(x -> 3 + 4/(x+1), 0.0, 10.0)
@code_warntype NHPP(x -> 3 + 4/(x+1), 0.0, 10.0)
@btime NHPP(x -> 3 + 4/(x+1), 0.0, 10.0)
