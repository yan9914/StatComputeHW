function mtn(data::Vector ; init = 1/2, ϵ = 1e-16)
    if length(data) != 3
        throw(DomainError())
    end
    θ₀ = init
    d = Inf
    while d > ϵ
        global x₃ = round(data[3] * (θ₀/4) / (1/2 + θ₀/4))
        d = abs(θ₀ - (data[2] + x₃) / (data[1] + data[2] + x₃))
        θ₀ = (data[2] + x₃) / (data[1] + data[2] + x₃)
    end
    return θ₀, x₃
end

mtn([38, 34, 125], init = 0.2)
mtn([38, 34, 125], init = 0.5)
mtn([38, 34, 125], init = 0.8)
