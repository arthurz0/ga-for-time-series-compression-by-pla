using DataFrames
using CSV
include("timeseries.jl")
include("loop.jl")
include("datasets.jl")

include("bellman.jl")
include("visvalingam.jl")
include("swingfilter.jl")
include("ramer_douglas_peucker.jl")
include("bruteforce.jl")

function createdf_algoeval()
    return DataFrame(
        algorithm = String[],
        dataset = String[],
        consideredLength = Integer[],
        remaining_points = Int[],
        remaining_points_out = Int[],
        mse = Float64[],
        baseTime = Float64[]
    )
end

function compress_log(algorithm::Function, datasetname, consideredLength, remainingPoints, loopIterations::Integer=0)
    df = createdf_algoeval()
    ts = datasetbyname(datasetname)[1:consideredLength]
    result = @timed algorithm(ts, remainingPoints)
    compressed, geno = result[1]
    push!(df, [
        string(algorithm), datasetname, consideredLength, remainingPoints, length(compressed), mse(compressed, ts), result[2]
    ])
    for i in 1:loopIterations
        compressed, geno = loop(compressed, ts) 
        push!(df, [
            string(algorithm) * ", $(i)x looped",
            datasetname, consideredLength, remainingPoints, length(compressed), mse(compressed, ts), result[2]
        ])
    end
    dir_experiment = "../experiments_ppsn/"
    filepath = dir_experiment * "algorithm_eval.csv"
    CSV.write(filepath, df, append=isfile(filepath))
end

function run_evaluation()
    algos = [bellman, bellman, bellman, bellman, bellman]
    dataset_remainingpoints = [("HandOutlines", 10000, 251), ("HandOutlines", 15000, 376), ("HandOutlines", 20000, 501), ("HandOutlines", 25000, 626)]
    for algo in algos
        for pair in dataset_remainingpoints
            compress_log(algo, pair[1], pair[2], pair[3], 0)
            println("$algo on $(pair[1]) with $(pair[2]) points finished")
        end
    end
end