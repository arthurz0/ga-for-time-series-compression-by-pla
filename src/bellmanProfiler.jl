using CSV
using DataFrames
include("datasets.jl")
include("bellman.jl")
include("loop.jl")
include("genetic.jl")

function createdf_bellman()
    return DataFrame(
        dataset = String[],
        length = Int[],
        remaining_points = Int[],
        mse = Float64[],
        time = Float64[],
        #time_tableonly = Float64[],
        mse_opys = Float64[],
        mses_10loops = Array{Float64, 1}[]#,
        #tablePrecision = Int[]
    )
end

dir_experiment = "../experiments2/"
file_bellman = dir_experiment * "bellman_notable.csv"

function profile_bellman(dataset::String, consideredLength::Integer, remainingPoints::Integer)
    println("start profiling $(dataset)[1:$(consideredLength)] down to $(remainingPoints)")
    ts = datasetbyname(dataset)[1:consideredLength]
    #GC.gc()
    #result_table = @timed calculate_error_edges(ts)
    #time_tableonly = result_table[2]
    #println("the table is being computed in $time_tableonly seconds")
    GC.gc()
    result = @timed bellman_notable(ts, remainingPoints)
    println("the complete algorithm is being computed in $(result[2]) seconds")
    compressed = result[1][1]
    mse_compressed = mse(compressed, ts)
    ind = Individual(result[1][2], ts)
    mse_opys = ind.mse
    mses_10loops = [0.0 for i in 1:10]
    for i in 1:10
        compressed, geno = loop(compressed, ts)
        mses_10loops[i] = mse(compressed, ts)
    end
    df = createdf_bellman()
    push!(df, [
        dataset, consideredLength, remainingPoints,
        mse_compressed, result[2], mse_opys, mses_10loops
    ])
    CSV.write(file_bellman, df, append=isfile(file_bellman))
    println("finished profiling")
end

function run_bellman_profiler()
    problems = [["subject103_6", 100000]]
	for problem in problems
		profile_bellman(problem[1], problem[2], floor(Int, problem[2]/100))
	end
end