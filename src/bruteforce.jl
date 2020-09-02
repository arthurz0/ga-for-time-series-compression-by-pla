include("timeseries.jl")
include("optimalYs.jl")

function nextsolution!(solution, maxVal)
    hasChanged = false
    i = 1
    for i in 1:length(solution)
        if i< length(solution) && solution[i] + 1 < solution[i+1] || i == length(solution) && solution[i] < maxVal
            solution[i] += 1
            for j in 1:i-1
                solution[j] = j+1
            end
            hasChanged = true
            break
        end
    end
    if hasChanged
        return solution
    else
        return nothing
    end
end

function bruteforce(ts, remainingPoints)
    bestError = Inf
    bestSolution = nothing
    bestCompressed = nothing
    bestGeno = nothing
    n = length(ts)-1

    solution = [i+1 for i in 1:remainingPoints-2]
    numbersOfTraversedSolutions = 0
    while solution != nothing
        numbersOfTraversedSolutions += 1
        times = vcat([timeof(ts[1])], map(x -> timeof(ts[x]), solution), [timeof(ts[end])])
        compressed::Timeseries = optimalValues(ts, times)
        mse_candidate = mse(compressed, ts)
        if mse_candidate < bestError
            bestError = mse_candidate
            bestSolution = copy(solution)
            bestCompressed = compressed
            bestGeno = falses(length(ts)-2)
            for i in solution
                bestGeno[i-1] = true
            end
        end
        solution = nextsolution!(solution, n)
    end
    if numbersOfTraversedSolutions != binomial(length(ts)-2, remainingPoints-2)
        @warn "expected $(binomial(length(ts)-2, remainingPoints-2)) solutions but only $(numbersOfTraversedSolutions) have been tested!"
    end
    return bestCompressed, bestGeno
end