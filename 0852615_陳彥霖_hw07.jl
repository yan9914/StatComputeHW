import Statistics: mean, var
function Run(arrival_t, service_t, n)
    depart_t = Inf
    N = 0
    i = 1
    new = 1
    spend = Array{Float64}(undef, n)
    while i in 1:n
        if min(depart_t, arrival_t[new]) == depart_t
            spend[i] = depart_t - arrival_t[i]
            N -= 1
            i += 1
            if N == 0
                depart_t = Inf
            else
                depart_t += service_t[i]
            end
        else
            if N == 0
                depart_t = arrival_t[new] + service_t[i]
            end
            N += 1
            new += 1
        end
    end
    return sum(spend)
end

# (a)
function f1(n::Int, arrival_rate, service_rate)
    arrival_t = cumsum(-log.(rand(n)) / arrival_rate)
    push!(arrival_t, Inf)
    service_t = -log.(rand(n)) / service_rate
    return Run(arrival_t, service_t, n)
end

using Plots
a1 = [f1(10, 2, 1) for _ in 1:100000]

mean(a1)
variance_1 = var(a1)
histogram(a1)


# (b)
function f2(n::Int, arrival_rate, service_rate)
    u = rand(n)
    arrival_t1 = cumsum(-log.(u) / arrival_rate)
    arrival_t2 = cumsum(-log.(1 .- u) / arrival_rate)
    push!(arrival_t1, Inf)
    push!(arrival_t2, Inf)
    u = rand(n)
    service_t1 = -log.(u) / service_rate
    service_t2 = -log.(1 .- u) / service_rate
    return (Run(arrival_t1, service_t1, n) + Run(arrival_t2, service_t2, n))/2
end

a2 = [f2(10, 2, 1) for _ in 1:100000]
mean(a2)
var(a2) / variance_1


# (c)
import LinearAlgebra.dot

function f3(n::Int, arrival_rate, service_rate)
    arrival_t = cumsum(-log.(rand(n)) / arrival_rate)
    push!(arrival_t, Inf)
    service_t = -log.(rand(n)) / service_rate
    return [Run(arrival_t, service_t, n), sum(service_t)]
end

nsim = 100000
x = Array{Float64}(undef, nsim)
y = Array{Float64}(undef, nsim)

for i in 1:nsim, (a, b) in tuple(f3(10, 2, 1))
    x[i], y[i] = a, b
end
c = -dot(x .- mean(x), y .- mean(y)) / dot(y .- mean(y), y .- mean(y))
q3 = x + c * (y .- 10)
mean(q3)
var(q3) / variance_1

# (d)
function f4(n::Int, arrival_rate, service_rate)
    arrival_t = cumsum(-log.(rand(n)) / arrival_rate)
    push!(arrival_t, Inf)
    service_t = -log.(rand(n)) / service_rate
    return [Run(arrival_t, service_t, n), sum(service_t) - sum(diff(arrival_t[1:length(arrival_t)-1]))]
end

nsim = 100000
x = Array{Float64}(undef, nsim)
y = Array{Float64}(undef, nsim)

for i in 1:nsim, (a, b) in tuple(f4(10, 2, 1))
    x[i], y[i] = a, b
end
c = -dot(x .- mean(x), y .- mean(y)) / dot(y .- mean(y), y .- mean(y))
q4 = x + c * (y .- (10 - 9/2))
mean(q4)
var(q4) / variance_1

# (e)
function Run2(arrival_t, service_t, n)
    depart_t = Inf
    N = 0
    i = 1
    new = 1
    output = Array{Float64}(undef, n)
    while i in 1:n
        if min(depart_t, arrival_t[new]) == depart_t
            i += 1
            N -= 1
            if N == 0
                depart_t = Inf
            else
                depart_t += service_t[i]
            end
        else
            if N == 0
                depart_t = arrival_t[new] + service_t[i]
            end
            N += 1
            output[new] = N
            new += 1
        end
    end
    return sum(output)
end

function f5(n::Int, arrival_rate, service_rate)
    arrival_t = cumsum(-log.(rand(n)) / arrival_rate)
    push!(arrival_t, Inf)
    service_t = -log.(rand(n)) / service_rate
    return Run2(arrival_t, service_t, n) / service_rate
end

a5 = [f5(10, 2, 1) for _ in 1:100000]
mean(a5)
var(a5) / variance_1
