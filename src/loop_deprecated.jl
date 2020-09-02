include("timeseries.jl")
include("util.jl")

function bestvalueanderror(leftNeighbor::Point, xPosition::Float64, rightNeighbor::Point, original::Timeseries)
    changeIndex = findwhile(p -> timeof(p) <= xPosition, original)
    firstSegment = original[1:changeIndex]
    secondSegment = original[changeIndex+1:end]
    dividend = 0.0
    divisor = 0.0
    x0 = timeof(leftNeighbor)
    y0 = valueof(leftNeighbor)
    xk = xPosition
    xn = timeof(rightNeighbor)
    yn = valueof(rightNeighbor)
    dividend += sum_float(p -> valueof(p) + (y0+valueof(p))*(timeof(p)-xk)/(xk-x0) + y0*(timeof(p)-xk)^2/(xk-x0)^2, firstSegment)
    dividend += sum_float(p -> valueof(p) + (yn+valueof(p))*(timeof(p)-xk)/(xk-xn) + yn*(timeof(p)-xk)^2/(xk-xn)^2, secondSegment)
    divisor += sum_float(p -> ((timeof(p)-xk)/(xk-x0)+1)^2, firstSegment)
    divisor += sum_float(p -> ((timeof(p)-xk)/(xk-xn)+1)^2, secondSegment)
    y = dividend/divisor
    errorOnInterval = 0.0
    errorOnInterval += sum_float(p -> (valueof(p) - y - (y-y0)*(timeof(p)-xk)/(xk-x0))^2, firstSegment)
    errorOnInterval += sum_float(p -> (valueof(p) - y - (y-yn)*(timeof(p)-xk)/(xk-xn))^2, secondSegment)
    return y, errorOnInterval
end

function bestErrorPoint(original::Timeseries, leftNeighbor::Point, rightNeighbor::Point)
    bestIndex = 0
    bestError = Inf
    bestY = 0.0
    for i in 1:length(original)
        y, error = bestvalueanderror(leftNeighbor, timeof(original[i]), rightNeighbor, original)
        if (error < bestError)
            bestError = error
            bestIndex = i
            bestY = y
        end
    end
    return bestIndex, bestY
end

"""Computes slope m and y-intercept b of the least square regression line with the additional constraint that it has to go through fixed"""
function bestline_fixpoint(ts::Timeseries, fixed::Point)
    #returns tuple m,b of line mx+b
    squaredXDifference = sum(p -> (timeof(p) - timeof(fixed))^2, ts)
    multipliedXYDifference = sum(p -> (timeof(p) - timeof(fixed)) * (valueof(p) - valueof(fixed)), ts)
    m = multipliedXYDifference / squaredXDifference
    b = valueof(fixed) - m*timeof(fixed)
    return m, b
end 

"""slow version of loop, doing an exhaustive search"""
function loop_deprecated(original::Timeseries , compressed::Timeseries)
    looped = copy(compressed)               #input should not change
    leftIndex=1
    rightIndex=1
    for i in 2:length(looped)-1
        leftIndex = findnext(p -> timeof(looped[i-1]) < timeof(p), original, leftIndex)
        rightIndex = findwhile(p -> timeof(p) < timeof(looped[i+1]), original, rightIndex)
        bestIndex, bestY = bestErrorPoint(original[leftIndex:rightIndex], looped[i-1], looped[i+1])
        bestIndex += leftIndex - 1
        looped[i] = (timeof(original[bestIndex]), bestY)
    end
    rightIndex = findwhile(p -> timeof(p) < timeof(looped[2]), original, 1)
    m, b = bestline_fixpoint(original[1:rightIndex], looped[2])
    looped[1] = (timeof(looped[1]), m * timeof(looped[1]) + b)
    leftIndex = findlast(p -> timeof(p) <= timeof(looped[end-1]), original) + 1
    m, b = bestline_fixpoint(original[leftIndex:end], looped[end-1])
    looped[end] = (timeof(looped[end]), m * timeof(looped[end]) + b)
    @assert length(looped) == length(compressed)
    return looped
end