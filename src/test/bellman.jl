using Test

include("../timeseries.jl")
include("../genetic.jl")
include("../bellman.jl")
include("../datasets.jl")

ts1 = Timeseries(collect(enumerate([1.0, 3.0, 5.0, 3.0, 1.0])))
result = bellman(ts1, 3, true)
@test result[1] == [(1.0, 1.0), (3.0, 5.0), (5.0, 1.0)]
@test result[3] == mse(result[1], ts1)
@test result[3] == 0.0 
@test result[2] == toGenotype(result[1], ts1)
@test bellman_notable(ts1, 3, true)[3] == result[3]

ts2 = Timeseries(collect(enumerate([0, 1.1, 1.9, 3.1, 5.0])))
result = bellman(ts2, 3, true)
@test result[1] == [(1.0, 0.0), (4.0, 3.1), (5.0, 5.0)]
@test result[3] ≈ mse(result[1], ts2)
@test result[2] == toGenotype(result[1], ts2)
@test bellman_notable(ts2, 3, true)[3] == result[3]

ts3 = Timeseries(collect(enumerate([0.0, 2.3, 3.2, 2.7, 2.0, 2.5, 1.0, 0.1])))
result = bellman(ts3, 5, true)
@test result[1] == [(1.0, 0.0), (3.0, 3.2), (5.0, 2.0), (6.0, 2.5), (8.0, 0.1)]
@test result[3] ≈ mse(result[1], ts3)
@test result[2] == toGenotype(result[1], ts3)
@test bellman_notable(ts3, 5, true)[3] == result[3]

ts4 = datasetbyname("apple")
@test bellman_notable(ts4, 100) == bellman(ts4, 100, true)[3]

function L2ErrorBetween(ts::Timeseries, i, j)
    if (j == i+1)
        return 0
    end
    startPoint = ts[i]
    endPoint = ts[j]
    deltaX = timeof(endPoint) - timeof(startPoint)
    deltaY = valueof(endPoint) - valueof(startPoint)
    error::Float64 = sum(i -> (valueof(ts, i) - (valueof(startPoint) + (timeof(ts, i) - timeof(startPoint)) / deltaX * deltaY))^2, i+1:j-1)
    return error
end

function calculate_error_edges_deprecated(ts::Timeseries)           #the less efficient version. remains in code for testing purposes
    numberVertices = length(ts)
    edgeArray = Array{Float64, 2}(undef, numberVertices, numberVertices)
    fill!(edgeArray, Inf)
    for i in 1:numberVertices
        for j in i+1:numberVertices
            edgeArray[i, j] = L2ErrorBetween(ts, i, j)
        end
    end
    return edgeArray
end

ts4 = datasetbyname("apple")[1:1000]
@test round.(calculate_error_edges(ts4), digits=14) ≈ transpose(round.(calculate_error_edges_deprecated(ts4), digits=14))