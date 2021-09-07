using Plots

function ep(;population::Int=4309, radius::Float64=1.5/1000, prob::Float64=0.01,
            step::Int=5000, span::Float64=0.7/1000, remove::Int=step*6)
    t = 0
    pS, pI, pR = [1-1/population], [1/population], [0.0]
    x = rand(population)
    y = rand(population)
    angle = rand(population)*2π
    condition = repeat([1], population)
    condition[rand(1:population)] = 2
    rm_time = repeat([Inf], population)
    rm_time[condition.==2] .= rand((remove-step):(remove+step))[1]
    while pI[length(pI)] != 0
        t += 1
        # 更新康復或死亡的人口狀態
        if any(rm_time .== t)
            condition[(1:population)[rm_time .== t]] .= 3
        end
        # 人群行走
        for i in (1:population)[condition .!= 3]
            angle[i] += rand(1)[1]*π/2 - π/4
            while !(0 < x[i]+span*cos(angle[i]) < 1 || 0 < y[i]+span*sin(angle[i]) < 1)
                angle[i] = rand(1)[1]*2π
            end
            x[i] = x[i] + span*cos(angle[i])
            y[i] = y[i] + span*sin(angle[i])
        end
        # 染病
        for i in (1:population)[condition .== 2]
            for j in (1:population)[condition .== 1]
                if (x[i]-x[j])^2 + (y[i]-y[j])^2 < radius^2
                    if rand(1)[1] < prob
                        condition[j] = 2
                        rm_time[j] = t + rand((remove-step):(remove+step))[1]
                    end
                end
            end
        end
        append!(pS, sum(condition.==1)/population)
        append!(pI, sum(condition.==2)/population)
        append!(pR, sum(condition.==3)/population)
    end
    gr()
    Plots.GRBackend()
    plot((0:t)/step, [repeat([1], t+1) 1 .- collect(pR) collect(pI)], fill = 0,
        color = permutedims(["Gray", "DarkSlateGray", "Firebrick"]), legend = :outertopright,
        label = ["Removed" "Susceptible" "Infectious"])
end

ep(radius = 1.5/1000)

ep(radius = 3.0/1000)
