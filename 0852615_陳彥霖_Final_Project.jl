using StatsBase
using StatsFuns
using Plots
using Statistics
using Random
using HypothesisTests

Random.seed!(1)

# problem 1

function FromBin(x ; n = 8, nBootstrap = 1e4)
    function max_dif(x, p)
        f = ecdf(x)
        sx = sort(x)
        a = abs.(binomcdf.(n, p, sx) .- f.(sx))
        return max(a...)
    end
    function SimBin(nsim, n, p)
        u = rand(nsim)
        tmp = map(x -> x .>= binomcdf.(n, p, -1:n), u)
        return findlast.(tmp) .- 1
    end
    l = length(x)
    p = sum(x) / (n * l)
    d = max_dif(x, p)
    Bdata = [SimBin(length(x), n, p) for _ in 1:nBootstrap]
    ds = max_dif.(Bdata, map(x -> sum(x) / (n * l), Bdata))
    p_value = mean(ds .> d)
    return (KS_statistic = d, p_value = p_value)
end

Ans = FromBin([6, 7, 3, 4, 7, 3, 7, 2, 6, 3, 7, 8, 2, 1, 3, 5, 8, 7]);
Ans

# problem 2

function HM_MCMC(d, lower, upper ; nsim = 1e4, iter = 1e3, init = 0)
    output = Array{Float64}(undef, Int(nsim))
    x = init
    l = upper - lower
    i = 0
    while i < iter * nsim
        y = rand() * 2l/3 + x - l/3
        if (y < 0) || (2 < y) ; continue ; end
        α = (prod(normpdf.(d.-y))*normpdf(y)) / (prod(normpdf.(d.-x))*normpdf(x))
        if rand() <= min(α, 1)
            x = y
        end
        i += 1
        if isinteger(i/iter)
            output[Int(i/iter)] = x
        end
    end
    return output
end

d2 = [0.8605, 1.2175, -0.9772, -0.0378, 2.9478,
        -0.2710, 0.0380, 0.6765, 2.7070, 0.5617,
        1.1110, 2.4136, 1.0066, 2.3637, 1.4502,
        0.2516, 0.3485, 1.6041, 1.2023, 1.6049]


Ans = HM_MCMC(d2, 0, 2, nsim = 1e4)
histogram(Ans)
mean(Ans)
std(Ans)

plot(bar(autocor(Ans, 1:20)), legend = nothing)
plot!([-1.96 * sqrt(1/1e4), 1.96 * sqrt(1/1e4)], seriestype = "hline", legned = nothing)
xaxis!("lag")
yaxis!("ACF")

print(LjungBoxTest(Ans, 20))


# Problem 3

function Gibbs_MCMC(α, β, λ ; nsim = 1e4, iter = 3e3, init = [3, 1/2, 5])
    output = Array{Float64}(undef, Int(nsim), 3)
    x = copy(init)
    i = 0
    while i < nsim * iter
        j = sample([1, 2, 3])
        if j == 1
            global x[1] = binominvcdf(x[3], x[2], rand())
        elseif j == 2
            global x[2] = betainvcdf(α + x[1], β + x[3] - x[1], rand())
        else
            global x[3] = poisinvcdf(λ*(1 - x[2]), rand()) + x[1]
        end
        i += 1
        if isinteger(i/iter)
            output[Int(i/iter), :] = [x...]
        end
    end
    return output
end

Ans = Gibbs_MCMC(5, 5, 7)
mean(Ans[:,1])
mean(Ans[:,2])
mean(Ans[:,3])


# Problem 4

u = Array{Float64}(undef, 11, 10)
for i in 1:11, j in setdiff(1:10, i)
    u[i,j] = rand()
end
u

function Simulated_Annealing2(m ; iter = 60, init = collect(0:10), C = 1)
    function reward(v)
        r = m[11, v[2]]
        r += sum([m[v[i], v[i+1]] for i in 2:10])
        return r
    end
    function find_neighbor(v)
        t = copy(v)
        (i, j) = rand(2:11, 2)
        t[i], t[j] = t[j], t[i]
        return t
    end
    x = copy(init)
    result = copy(x)
    a = reward(result)
    k = 1
    while k <= iter
        y = find_neighbor(x)
        λ = C * log(k + 1)
        if rand() <= exp(λ*(reward(y) - reward(x)))
            x = copy(y)
            if reward(y) >= a
                a = reward(y)
                result = copy(y)
            end
            k += 1
        end
    end
    return (path = result, reward = a)
end

a4 = Simulated_Annealing2(u, iter = 10000, C = 1);
a4


u = reshape(rand(11*10*Int(1e6)), (11,10,Int(1e6)))
@time p4 = [Simulated_Annealing2(u[:,:,i], iter = 10000, C = 1)[2] for i in 1:Int(1e6)]
mean(p4)
var(p4)

# Problem 5

function EM_algorithm(v, ξ₀, λ₀ ; ϵ = 1e-16)
    ξ = copy(ξ₀)
    λ = copy(λ₀)
    N = sum(v)
    t = 0
    d = Inf
    table = [Text.(["t" "ξ" "λ" "n_A" "n_B"]);Array{Any}(nothing, 6, 5)]
    while d > ϵ
        n_A = v[1]*ξ / (ξ+(1-ξ)*exp(-λ))
        n_B = v[1] - n_A
        if t <= 5
            table[t+2,:] = Any[Int8(t) ξ λ n_A n_B]
        end
        t += 1
        D = [ξ, λ] - [n_A/N, sum(collect(1:6).*v[2:end])/(N-n_A)]
        d = sqrt(sum(D.^2))
        (ξ, λ) = [ξ, λ] - D
    end
    table
end

data = [3062, 587, 284, 103, 33, 4, 2]
EM_algorithm(data, 0.75, 0.4)








EM_algorithm(data, 0.5, 0.6)
