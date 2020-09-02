using DataStructures
include("timeseries.jl")

module _Deletable
    mutable struct DeletablePoint
        index::Integer
        prev::Union{DeletablePoint, Nothing}
        next::Union{DeletablePoint, Nothing}
    end
end

function visvalingam_area_error_helper(left::Point, middle::Point, right::Point)
    dy = abs(valueof(middle) - (valueof(left) + (timeof(middle) - timeof(left)) * slope(left, right)))
    return dy * (timeof(right) - timeof(left)) / 2.0
end

function visvalingam_approximated_se_helper(left::Point, middle::Point, right::Point)
    slope_lm = slope(left, middle)
    slope_mr = slope(middle, right)
    slope_lr = slope(left, right)
    dl = timeof(middle) - timeof(left)
    dr = timeof(right) - timeof(middle)
    left_error = (slope_lm - slope_lr)^2 / 3.0 * dl^3
    right_error = (slope_mr - slope_lr)^2 / 3.0 * dr^3
    return left_error + right_error
end

function visvalingam_area(original::Timeseries, remainingPoints::Integer)
    return visvalingam(original::Timeseries, remainingPoints::Integer, visvalingam_area_error_helper)
end

function visvalingam_squarederror(original::Timeseries, remainingPoints::Integer)
    return visvalingam(original::Timeseries, remainingPoints::Integer, visvalingam_approximated_se_helper)
end

function visvalingam(original::Timeseries, remainingPoints::Integer, errorFunction)
    function priority(deletablePoint)
        return errorFunction(original[deletablePoint.prev.index], original[deletablePoint.index], original[deletablePoint.next.index])
    end

    @assert remainingPoints >= 2
    start = _Deletable.DeletablePoint(1, nothing, nothing)
    prev = start
    pq = PriorityQueue{_Deletable.DeletablePoint, Float64}()
    for i in 2:length(original)
        current = _Deletable.DeletablePoint(i, prev, nothing)
        prev.next = current
        prev = current
    end
    current = start.next
    while current.next != nothing
        enqueue!(pq, current, priority(current))
        current = current.next
    end
    for i in 1:length(original) - remainingPoints
        removed = dequeue!(pq)
        removed.next.prev = removed.prev
        removed.prev.next = removed.next
        if removed.prev.prev != nothing pq[removed.prev] = priority(removed.prev) end
        if removed.next.next != nothing pq[removed.next] = priority(removed.next) end
    end
    result::Timeseries = []
    geno = falses(length(original) - 2)
    current = start
    while current != nothing
        push!(result, original[current.index])
        if 1 < current.index < length(original)
            geno[current.index - 1] = true
        end
        current = current.next
    end
    return result, geno
end