using Statistics

include("timeseries.jl")
include("util.jl")

@views function optimalValues(original::Timeseries, times::Array{Float64, 1})  
    @assert timeof(original[1]) == times[1] && timeof(original[end]) == times[end]
    if isempty(times)
        return Timeseries([])
    end
    optimalValuesResults::Array{Array{Float64, 1}, 1} = []
    originalStartIndex::Integer = 1
    originalEndIndex::Integer = 1
    timesStartIndex::Integer = 1
    for i in 1:length(times)-1
        #originalEndIndex = findwhile(p -> timeof(p) <= times[i], original, originalEndIndex)
        while timeof(original[originalEndIndex]) <= times[i]
            originalEndIndex += 1
        end
        originalEndIndex -= 1
        if times[i+1] <= timeof(original, originalEndIndex + 1)
            push!(optimalValuesResults, optimalYsBlock(original[originalStartIndex:originalEndIndex], times[timesStartIndex:i]))
            originalStartIndex = originalEndIndex + 1
            originalEndIndex = originalStartIndex
            timesStartIndex = i + 1
        end
    end
    push!(optimalValuesResults, optimalYsBlock(original[originalStartIndex:end], times[timesStartIndex:end]))
    optimalVs = vcat(optimalValuesResults...)
    return Timeseries(collect(zip(times, optimalVs)))
end

function optimalYsBlock(original::AbstractTimeseries, times::AbstractArray{Float64, 1})
    if isempty(original)
        return map(x -> NaN, times)
    end
    if length(times) == 1
        return [mean(valuesof(original))]
    end
    equations::Matrix{Float64} = Matrix{Float64}(undef, 4, length(times))        #left, own, right, bias
    fill!(equations, 0.0)
    originalStartIndex = 1
    originalEndIndex = 1
    for i in 1:length(times) - 1
        #originalEndIndex = findwhile(p -> timeof(p) < times[i+1], original, originalEndIndex)
        while timeof(original[originalEndIndex]) < times[i+1]
            originalEndIndex += 1
        end
        originalEndIndex -= 1
        for j in originalStartIndex:originalEndIndex
            p = original[j]
            normalized_dist = (timeof(p) - times[i]) / (times[i+1] - times[i])
            equations[2, i] += (1-normalized_dist)^2
            equations[2, i+1] += normalized_dist^2
            equations[1, i+1] += (1-normalized_dist)*normalized_dist
            equations[3, i] += (1-normalized_dist)*normalized_dist
            equations[4, i] -= valueof(p) * (1-normalized_dist)
            equations[4, i+1] -= valueof(p) * normalized_dist
        end
        originalStartIndex = originalEndIndex + 1
    end
    if timeof(original[end]) == times[end]
        equations[2, end] += 1
        equations[4, end] += -valueof(original[end])
    end
    #equations[1, 1] = -1.0
    #equations[2, 1] = 0.0
    equations[3, 1] = - equations[3,1] / equations[2,1]
    equations[4, 1] = - equations[4,1] / equations[2,1]
    for i in 2:length(times)
        interown = equations[2,i] + equations[1,i] * equations[3,i-1]
        interbias = equations[4,i] + equations[1,i] * equations[4,i-1]
        #equations[1,i] = -1.0
        #equations[2,i] = 0.0
        equations[3,i] = - equations[3,i] / interown
        equations[4,i] = - interbias / interown
    end
    results::Array{Float64, 1} = Array{Float64, 1}(undef, length(times))
    results[end] = equations[4, end]
    for i in length(times)-1:-1:1
        results[i] = equations[3,i] * results[i+1] + equations[4,i]
    end
    return results
end