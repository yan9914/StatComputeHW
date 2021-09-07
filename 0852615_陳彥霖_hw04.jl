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

using LinearAlgebra

C = [1 0.5 0.5 ; 0.5 1 0.5 ; 0.5 0.5 1]

A = cholesky(C)

X = Array{Float64}(undef, 10000, 3)
for i in 1:10000
    X[i,:] = A.L * Box_Muller(3)
end
X

using Statistics
[mean(X[:,i]) for i in 1:3]




[var(X[:,i]) for i in 1:3]




[cor(X[:,i], X[:,j]) for (i,j) = Iterators.product(1:3,1:3)]






import StatsFuns.normcdf

Y = -log.(1 .- normcdf.(0, 1, X))

using Plots
using StatsBase
f1 = ecdf(Y[:,1])
f2 = ecdf(Y[:,2])

gr()
Plots.GRBackend()

scatter(f1.(Y[:,1]), f2.(Y[:,2]), markerstrokewidth = 0, markersize = [1.5], legend = false)
xlabel!("y_1"); ylabel!("y_2")
title!("Rank plot")
