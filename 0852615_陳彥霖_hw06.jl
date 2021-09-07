function f(arrival_rate, service_rate, endTime)
    t = 0   # 時間軸
    n = 0   # 系統總人數
    depart = Inf    # 正在被服務的人的離開時間
    arrival = []    # 現在系統中的客人各自的到達時間
    spend = []      # 各個客人花在系統的時間
    next = t - log(rand(1)[1]) / arrival_rate   # 下一位潛在客戶到達的時間

    # 當還沒超過結束時間
    while min(depart, next, endTime) == depart || min(depart, next, endTime) == next

        # 判斷下次事件是否為: 客人被服務完而離開
        if min(depart, next, endTime) == depart
            t = depart
            n -= 1
            push!(spend, t - arrival[1])
            deleteat!(arrival, 1)
            if n == 0
                depart = Inf
            else
                depart = t - log(rand(1)[1]) / service_rate
            end

        # 下個事件為: 客人準備進門
        else
            t = next
            # 系統內少於或等於三人才加入
            if n <= 3
                if n == 0
                    depart = t - log(rand(1)[1]) / service_rate
                end
                n += 1
                push!(arrival, t)
            end
            next = t - log(rand(1)[1]) / arrival_rate
        end
    end

    # 當時間超過結束時間, 不再有人進門, 把剩下的客人服務完即可
    while depart != Inf
        t = depart
        n -= 1
        push!(spend, t - arrival[1])
        deleteat!(arrival, 1)
        if n == 0
            depart = Inf
        else
            depart = t - log(rand(1)[1]) / service_rate
        end
    end
    return spend
end

import Statistics: mean, var

sample = f(4, 4.2, 8)
mean(sample)


using Plots
histogram([mean(f(4, 4.2, 8)) for _ in 1:10000])

function BS_mse(x, B)
    a = mean.([x[rand(1:length(x), length(x))] for _ in 1:B])
    bias = mean(a) - mean(sample)
    variance = var(a)
    return bias^2 + variance
end

BS_mse(sample, 100000)
