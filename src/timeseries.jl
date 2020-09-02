include("util.jl")

const Point = Tuple{Float64, Float64}   #why not a point? because (1, 2) is a tuple
const Timeseries = Array{Point, 1}
const AbstractTimeseries = AbstractArray{Point, 1}

#The points of a timeseries have to be ordered according to their time

function timeof(p::Point)
    return p[1]
end

function valueof(p::Point)
    return p[2]
end

function addvalue(p::Point, x::Number)
    return (p[1], p[2]+x)
end

function timeof(ts::AbstractTimeseries, index::Int64)
return timeof(ts[index])
end

function valueof(ts::AbstractTimeseries, index::Int64)
return valueof(ts[index])
end

function timesof(ts::AbstractTimeseries)
    return map(p::Point -> timeof(p), ts)
end

function valuesof(ts::AbstractTimeseries)
    return map(p::Point -> valueof(p), ts)
end

function plotTs(ts::AbstractTimeseries; args...)
    plot(timesof(ts), valuesof(ts); args...)
end

function scatterTs(ts::AbstractTimeseries; args...)
    scatter(timesof(ts), valuesof(ts); args...)
end

function slope(a::Point, b::Point)
    return (valueof(b) - valueof(a)) / (timeof(b) - timeof(a))
end

function mse(approximation::Timeseries, original::Timeseries)
    if isempty(approximation)
        return mse([(0.0, 0.0)], original)
    end
    @assert(timeof(approximation[end]) <= timeof(original[end]), "method not suitable for this special case")
    error = 0.0
    startPoint = approximation[1]
    toIndex = findwhile(p -> timeof(p) <= timeof(startPoint), original)
    error += sum_float(p -> (valueof(p) - valueof(startPoint))^2, view(original, 1:toIndex))
    for i in 1:length(approximation)-1
        startPoint = approximation[i]
        endPoint = approximation[i+1]
        m = (valueof(endPoint) - valueof(startPoint)) / (timeof(endPoint) - timeof(startPoint))
        fromIndex = toIndex + 1
        toIndex = findwhile(p -> timeof(p) <= timeof(endPoint), original, fromIndex)
        error += sum_float(p -> (valueof(p) - (m * (timeof(p) - timeof(startPoint)) + valueof(startPoint)))^2, view(original, fromIndex:toIndex))
    end
    endPoint = approximation[end]
    fromIndex = toIndex + 1
    error += sum_float(p -> (valueof(p) - valueof(endPoint))^2, view(original, fromIndex:length(original)))
    return error/length(original)
end