include("timeseries.jl")

function swing_filter_number(ts::Timeseries, remainingPoints::Integer)
    minEps = 0.001
    maxEps = max(valuesof(ts)...) - min(valuesof(ts)...)
    currentEps = (minEps + maxEps) / 2
    for i in 1:20
        currentEps = (minEps + maxEps) / 2
        compressed = swing_filter(ts, currentEps)
        if length(compressed[1]) == remainingPoints
            return compressed
        end
        if length(compressed[1]) > remainingPoints
            minEps = currentEps
        else
            maxEps = currentEps
        end
    end
    return swing_filter(ts, currentEps)
end

function swing_filter_number_better(ts::Timeseries, remainingPoints::Integer)
    minEps = 0.0
    maxEps = max(valuesof(ts)...) - min(valuesof(ts)...)
    currentEps = (minEps + maxEps) / 2
    for i in 1:20
        currentEps = (minEps + maxEps) / 2
        compressed = swing_filter(ts, currentEps)
        if length(compressed[1]) == remainingPoints
            return compressed
        end
        if length(compressed[1]) > remainingPoints
            minEps = currentEps
        else
            maxEps = currentEps
        end
    end
    return swing_filter(ts, currentEps)
end

function epsilon_sampling(ts::Timeseries, ε::Float64)
    compressed::Timeseries = []
    geno = falses(length(ts)-2)
    i = 1
    while i < length(ts)
        push!(compressed, ts[i])
        if i > 1
            geno[i-1] = true
        end
        from = ts[i]
        i += 1
        minGradient = -Inf
        maxGradient = Inf
        while i <= length(ts)
            new_point = ts[i]
            newMinGradient = slope(from, addvalue(new_point, -ε))
            newMaxGradient = slope(from, addvalue(new_point, ε))
            minGradient = max(minGradient, newMinGradient)
            maxGradient = min(maxGradient, newMaxGradient)
            if minGradient ≤ slope(from, new_point) ≤ maxGradient
                i += 1
            else
                i -= 1
                break
            end
        end
    end
    if timeof(compressed[end]) < timeof(ts[end])
        push!(compressed, ts[end])
    end
    return compressed, geno
end

function bestline_fixpoint(ts::Timeseries, fixed::Point)
    squaredXDifference = sum(p -> (timeof(p) - timeof(fixed))^2, ts)
    multipliedXYDifference = sum(p -> (timeof(p) - timeof(fixed)) * (valueof(p) - valueof(fixed)), ts)
    m = multipliedXYDifference / squaredXDifference
    b = valueof(fixed) - m*timeof(fixed)
    return m, b
end

function swing_filter(ts::Timeseries, ε::Float64)
    compressed::Timeseries = []
    geno = falses(length(ts)-2)
    from = ts[1]
    push!(compressed, ts[1])
    indexLastSent = 1
    minGradient = -Inf
    maxGradient = Inf
    for i in 2:length(ts)
        new_point = ts[i]
        newMinGradient = max(slope(from, addvalue(new_point, -ε)), minGradient)
        newMaxGradient = min(slope(from, addvalue(new_point, ε)), maxGradient)
        if newMaxGradient < newMinGradient
            m, b = bestline_fixpoint(ts[indexLastSent+1:i-1], from)
            m = max(min(m, maxGradient), minGradient)
            from = (timeof(ts[i-1]), valueof(from) + m * (timeof(ts[i-1]) - timeof(from)))
            push!(compressed, from)
            geno[i-2] = true
            indexLastSent = i-1
            minGradient = slope(from, addvalue(new_point, -ε))
            maxGradient = slope(from, addvalue(new_point, ε))
        else
            minGradient = newMinGradient
            maxGradient = newMaxGradient
        end
    end
    if timeof(compressed[end]) < timeof(ts[end])
        push!(compressed, ts[end])
    end
    return compressed, geno
end