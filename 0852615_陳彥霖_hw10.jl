f(x) = 0.2 + x[1]^2 + x[2]^2 - 0.1cos(6pi*x[1]) - 0.1cos(6pi*x[2])

function mini(fun, k ; C = 1)
    f(x) = -fun(x)
    x = [100.,100]
    n = 1
    obj = Array{Array}(undef, k)
    value = Array{Float64}(undef, k)
    while n <= k
        λ = C * log(n + 1)
        i = rand([1, 2])
        tmp = randn(1)[1]
        y = copy(x)
        y[i] = tmp
        if f(y) >= f(x)
            x = y
            obj[n] = y
            value[n] = f(y)
            n += 1
        else
            if rand(1)[1] <= exp(λ*(f(y)-f(x)))
                x = y
                obj[n] = y
                value[n] = f(y)
                n += 1
            end
        end
    end
    i = argmax(value)
    return (obj = obj[i], value = -value[i])
end

a = mini(f, 100000, C = 10)
a.obj


a.value
