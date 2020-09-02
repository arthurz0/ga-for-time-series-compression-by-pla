using DataStructures
include("timeseries.jl")

function ramer_douglas_peucker(ts::Timeseries, numberOfPointsToRemember::Int)
    #potentially unsafe if you suffer from bad luck (PriorityQueue can't handle duplicate keys)
    function getMaxError(a::Int64, b::Int64)
        startPoint = ts[a]
        endPoint = ts[b]
        gradient = (valueof(endPoint) - valueof(startPoint)) / (timeof(endPoint) - timeof(startPoint))
        index = 0
        maxDist = 0
        for i in a:b
            p = ts[i]
            dist = abs(valueof(p) - gradient * (timeof(p) - timeof(startPoint)) - valueof(startPoint))
            if dist > maxDist
                maxDist = dist
                index = i
            end
        end
        return index, maxDist
    end

    result = Timeseries()
    geno = falses(length(ts)-2)
    push!(result, ts[1])
    push!(result, ts[end])
    numberOfPointsToRemember = numberOfPointsToRemember - 2
    pq = PriorityQueue()
    i, d = getMaxError(1, length(ts))
    if i > 0
        enqueue!(pq, [1, i, length(ts)] => -d)
    end
    while numberOfPointsToRemember > 0 && length(pq) > 0
        x = dequeue!(pq)
        push!(result, ts[x[2]])
        geno[x[2]-1] = true
        numberOfPointsToRemember = numberOfPointsToRemember - 1
        i, d = getMaxError(x[1], x[2])
        if i > 0
            enqueue!(pq, [x[1], i, x[2]] => -d)
        end
        i, d = getMaxError(x[2], x[3])    
        if i > 0
            enqueue!(pq, [x[2], i, x[3]] => -d)
        end
    end
    result = sort(result)
    return result, geno
end
