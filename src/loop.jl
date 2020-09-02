include("timeseries.jl")

"""A helper function for the loop heuristic."""
function bestindexandvalue(leftNeighbor::Point, rightNeighbor::Point, original::AbstractTimeseries, se_left::Matrix{Float64}, se_right::Matrix{Float64}, candidates::UnitRange{Int64}=1:length(original))
    @assert !isempty(original)
    mse_polynomials_by_slope(original, leftNeighbor, se_left)
    rightside_mse_polynomials_by_slope(original, rightNeighbor, se_right)
    best_index = -1
    best_value = NaN
    best_se = Inf
    for i in candidates
        dl = timeof(original[i]) - timeof(leftNeighbor)
        dr = timeof(rightNeighbor) - timeof(original[i])
        vl = valueof(leftNeighbor)
        vr = valueof(rightNeighbor)
        #l = se_left[:, i]
        l1 = se_left[1, i]
        l2 = se_left[2, i]
        l3 = se_left[3, i]
        #r = se_right[:, i]
        r1 = se_right[1, i]
        r2 = se_right[2, i]
        r3 = se_right[3, i]
        y = (l1*2vl/dl^2 - l2/dl + r1*2vr/dr^2 + r2/dr) / (2l1/dl^2 + 2r1/dr^2)
        sl = (y - vl)/dl
        sr = (vr - y)/dr
        squared_error = l1 * sl^2 + l2 * sl + l3 + r1 * sr^2 + r2 * sr + r3
        if squared_error < best_se
            best_index = i
            best_se = squared_error
            best_value = y
        end
    end
    return best_index, best_value
end

function rightside_mse_polynomials_by_slope(original::AbstractTimeseries, point::Point, polynomials::Matrix{Float64})
    polynomials[1,length(original)] = 0.0
    polynomials[2,length(original)] = 0.0
    polynomials[3,length(original)] = 0.0
    for i in 1:length(original)-1
        dx = timeof(original[i+1]) - timeof(point)
        dy = valueof(original[i+1]) - valueof(point)
        polynomials[1,i] = dx^2  
        polynomials[2,i] = -2.0dy*dx   
        polynomials[3,i] = dy^2                           
    end
    for i in length(original)-1:-1:1
        polynomials[1,i] += polynomials[1,i+1]
        polynomials[2,i] += polynomials[2,i+1]
        polynomials[3,i] += polynomials[3,i+1]
    end
    return polynomials
end

function mse_polynomials_by_slope(original::AbstractTimeseries, point::Point, polynomials::Matrix{Float64})
    for i in 1:length(original)
        dx = timeof(original[i]) - timeof(point)
        dy = valueof(original[i]) - valueof(point)
        polynomials[1,i] = dx^2  
        polynomials[2,i] = -2.0dy*dx   
        polynomials[3,i] = dy^2                           
    end
    for i in 2:length(original)
        polynomials[1,i] += polynomials[1,i-1]
        polynomials[2,i] += polynomials[2,i-1]
        polynomials[3,i] += polynomials[3,i-1]
    end
    return polynomials
end

function bestline_fixpoint(ts::AbstractTimeseries, fixed::Point)
    squaredXDifference = sum(p -> (timeof(p) - timeof(fixed))^2, ts)
    multipliedXYDifference = sum(p -> (timeof(p) - timeof(fixed)) * (valueof(p) - valueof(fixed)), ts)
    m = multipliedXYDifference / squaredXDifference
    b = valueof(fixed) - m*timeof(fixed)
    return m, b
end

"""This code relies on the assumption that times(compressed) is a subset of times(original)."""
function loop(compressed::Timeseries, original::Timeseries, se_left::Matrix{Float64}=Matrix{Float64}(undef, 3, length(original)), se_right::Matrix{Float64}=Matrix{Float64}(undef, 3, length(original)))
    looped::Timeseries = copy(compressed)
    return loop!(looped, original, se_left, se_right)
end

"""Inplace version of loop"""
@views function loop!(compressed::Timeseries, original::Timeseries, se_left::Matrix{Float64}, se_right::Matrix{Float64})
    geno = falses(length(original) - 2)
    leftIndex = 1
    rightIndex = 1
    middleIndex = findwhile(p -> timeof(p) <= timeof(compressed, 2), original, 1)
    for i in 2:length(compressed)- 1
        rightIndex = findwhile(p -> timeof(p) <= timeof(compressed, i+1), original, middleIndex+1)
        bestIndex, bestY = bestindexandvalue(compressed[i-1], compressed[i+1], original[leftIndex+1:rightIndex-1], se_left, se_right)
        bestIndex += leftIndex
        compressed[i] = (timeof(original, bestIndex), bestY)
        geno[bestIndex - 1] = true
        leftIndex = bestIndex
        middleIndex = rightIndex
    end
    rightIndex = findwhile(p -> timeof(p) <= timeof(compressed, 2), original)
    m, b = bestline_fixpoint(original[1:rightIndex-1], compressed[2])
    compressed[1] = (timeof(compressed, 1), m*timeof(compressed, 1) + b)
    m, b = bestline_fixpoint(original[leftIndex+1:end], compressed[end-1])
    compressed[end] = (timeof(compressed, length(compressed)), m*timeof(compressed, length(compressed)) + b)
    return compressed, geno
end