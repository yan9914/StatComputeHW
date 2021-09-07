using Distributions
using StatsBase
using StatsFuns
import Statistics.mean

function FromExpDist(x, nBootstrap)
    # 估計參數
    θ = fit_mle(Exponential, x).θ
    # 找 KS 統計量
    function max_dif(x, θ)
        f = ecdf(x)
        sx = sort(x)
        a = abs.(gammacdf.(1, θ, sx) .- f.([0, sx[1:length(sx)-1]...]))
        return max(a..., abs.(gammacdf.(1, θ, sx) .- f.(sx))...)
    end
    d = max_dif(x, θ)
    # Bootstrap 生成資料, 並計算統計量
    Bdata = [-θ*log.(rand(length(x))) for _ in 1:nBootstrap]
    ds =  max_dif.(Bdata, map(x -> fit_mle(Exponential, x).θ, Bdata))
    # 計算模擬的 p-value
    p_value = mean(ds .> d)
    return p_value
end

FromExpDist([122, 133, 106, 128, 135, 126], 1e6)
