import SpecialFunctions.gamma
using StatsFuns
using Statistics
using Plots
using JuMP
using Ipopt
using PoissonRandom
using Distributions

# 1

function p1b(nsim = 10000)
    function find_y()
        y = Array{Float64}(undef, 3)
        y[3] = -1
        while y[3] < 0
            y[1:2] = rand(2)
            y[3] = 1 - sum(y[1:2])
        end
        return y
    end
    f(x) = gamma(sum([2, 3, 4])) / prod(gamma.([2, 3, 4])) * x[1] * x[2]^2 * x[3]^3
    output = Array{Float64}(undef, nsim, 3)
    N = Array{Float64}(undef, nsim)

    # Step1
    m = Model(Ipopt.Optimizer)
    @variable(m, 0 <= x[1:3] <= 1)
    @NLobjective(m, Max, 40320 / 12 * x[1] * x[2]^2 * x[3]^3)
    @constraint(m, x[1] + x[2] + x[3] == 1)
    optimize!(m)
    c = JuMP.objective_value(m)

    i = 1
    j = 1
    while i <= nsim
        y = find_y()        # Step2 + Step3
        U₂ = rand(1)[1]     # Step4
        if U₂ <= f(y)/c     # Step5
            output[i,:] = y
            N[i] = j
            i += 1
            j = 1
        else
            j += 1
        end
    end
    return (output = output, N = N)
end

a = p1b();

gr()
histogram(a.output[:,1])
histogram(a.output[:,2])
histogram(a.output[:,3])

mean(a.output[:,1]), mean(a.output[:,2]), mean(a.output[:,3])
cov(a.output)



mean(a.N), var(a.N)

function p1d(nsim = 10000)
    output = Array{Float64}(undef, nsim, 3)
    for i in 1:nsim
        tmp = -log.(rand(9))
        y = [sum(tmp[1:2]), sum(tmp[3:5]), sum(tmp[6:9])]
        output[i,:] = y / sum(y)
    end
    return output
end

a = p1d()

histogram(a[:,1])
histogram(a[:,2])
histogram(a[:,3])

mean(a[:,1]), mean(a[:,2]), mean(a[:,3])


# 2

p2a(nsim = 1e5) = sum.(randn.([pois_rand(5) for _ in 1:nsim])) .> 10
a = p2a()
mean(a)
var(a)


function p2b(nsim = 1e5)
    N = [pois_rand(5) for _ in 1:nsim]
    x = sum.(randn.(N)) .> 10
    return x .- cov(Int64[x...], N)/var(N) * (N .- 5)
end

a = p2b()
mean(a)
var(a)

function p2c(nsim = 1e5)
    N = [pois_rand(5) for _ in 1:nsim]
    return normccdf.(0, sqrt.(N), 10)
end

a = p2c()
mean(a)
var(a)

function p2d(nsim = 1e5)
    N = [pois_rand(10) for _ in 1:nsim]
    x = sum.(randn.(N)) .> 10
    return x .* poispdf.(5, N) ./ poispdf.(10, N)
end

a = p2d()
mean(a)
var(a)

function p2e(nsim = 100000)
    output = Array{Float64}(undef, nsim)
    N = [pois_rand(5) for _ in 1:nsim]
    x = map(x -> x .+ 2, randn.(N))
    for i in 1:nsim
        output[i] = (sum(x[i]) > 10) .* prod(normpdf.(x[i])) ./ prod(normpdf.(2, 1, x[i]))
    end
    return output
end

a = p2e()
mean(a)
var(a)

function p2f(nsim = 100000)
    output = Array{Float64}(undef, nsim)
    N = [pois_rand(10) for _ in 1:nsim]
    x = map(x -> x .+ 2, randn.(N))
    for i in 1:nsim
        r₁ = prod(normpdf.(x[i])) / prod(normpdf.(2, 1, x[i]))
        r₂ = poispdf(5, N[i]) / poispdf(10, N[i])
        output[i] = (sum(x[i]) > 10) .* r₁ .* r₂
    end
    return output
end

a = p2f()
mean(a)
var(a)


# 3

data = [0.0839, 0.0205, 0.3045, 0.7816, 0.0003,
        0.0095, 0.4612, 0.9996, 0.9786, 0.7580,
        0.0002, 0.7310, 0.0777, 0.4483, 0.4449,
        0.7943, 0.1447, 0.0431, 0.8621, 0.3273]

function p3a(data)
    x = [deleteat!(copy(data), i) for i in 1:length(data)]
    p = histogram(vcat(x...))
    m = median.(x)
    sd = sqrt(sum((m .- mean(m)).^2) * (length(data)-1) / length(data))
    CI = (mean(data) - 1.96 * sd,
          mean(data) + 1.96 * sd)
    return (CI = CI, plot = p)
end

a = p3a(data);
a.plot
a.CI

function p3b(data)
    x = [data[rand(1:length(data), length(data))] for _ in 1:1e4]
    p = histogram(vcat(x...))
    m = median.(x)
    sd = std(m)
    θ₀ = median(data)
    z = (m .- θ₀) ./ std.(x)
    t = sort(z)[[9750, 250]]
    CI = Tuple(θ₀ .- t * std(data))
    return (CI = CI, plot = p)
end

a = p3b(data);
a.plot
a.CI

function p3c(data)
    f = fit(Beta, data)
    x = [betainvcdf.(f.α, f.β, rand(length(data))) for _ in 1:1e4]
    p = histogram(vcat(x...))
    m = sort(median.(x))
    CI = (m[250], m[9750])
    return (CI = CI, plot = p)
end

a = p3c(data);
a.plot
a.CI
