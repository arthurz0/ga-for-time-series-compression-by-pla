using CSV
using DataFrames
using Base.Threads

include("timeseries.jl")
include("experiment_config.jl")
include("genetic.jl")
include("fitness.jl")
include("createGenes.jl")
include("selection.jl")
include("datasets.jl")
include("ramer_douglas_peucker.jl")
include("swingfilter.jl")
include("visvalingam.jl")

using Dates

dir_experiment = "../paper_eval/"

function createdf_total()
    return hcat(DataFrame(run_id = Int[]),
        createdf_run())
end

function createdf_run()
    return DataFrame(
        generation = Int[],
        fitness_std = Float64[],
        geno_diversity = Float64[],
        best_valid_numberOfPoints = Int[],
        best_valid_mse = Float64[],
        best_valid_fitness = Float64[],
        best_overall_numberOfPoints = Int[],
        best_overall_mse = Float64[],
        best_overall_fitness = Float64[],
        avg_numberOfPoints = Float64[],
        avg_mse = Float64[],
        avg_fitness = Float64[],
        number_valid_individuals = Int[],
        total_time = Float64[]
    )
end

function log_population!(df::DataFrame, population::Array{Individual,1}, generationN::Integer, remainingPoints::Integer, totalTime::Float64)
    best_valid_ind = Individual(trues(0), [(0.0, 0.0), (1.0, 1.0)], Inf, Inf)
    most_fit_ind = population[1]
    for ind in population
        if ind.fitness < most_fit_ind.fitness
            most_fit_ind = ind
        end
        if length(ind.pheno) <= remainingPoints && ind.mse < best_valid_ind.mse
            best_valid_ind = ind
        end
    end
    mean_NPoints = float(sum(ind->length(ind.pheno), population)) / length(population)
    mean_fitness = sum(ind->ind.fitness, population) / length(population)
    mean_mse = sum(ind->ind.mse, population) / length(population)
    numberValid = sum(ind -> length(ind.pheno) <= remainingPoints, population)
    push!(df, [
            generationN, fitness_std(population), genotype_diversity(population),
            length(best_valid_ind.pheno), best_valid_ind.mse, best_valid_ind.fitness, 
            length(most_fit_ind.pheno), most_fit_ind.mse, most_fit_ind.fitness, 
            mean_NPoints, mean_mse, mean_fitness, numberValid, totalTime
        ])
    return df
end

function createOffspringGenes(parents::Array{Individual, 1}, config::ExperimentConfig)
    offspring::Array{BitArray, 1} = Array{BitArray, 1}(undef, length(parents))
    for i in 1:2:length(parents)
        geno1, geno2 = crossover_onepoint(parents[i].geno, parents[i+1].geno)
        offspring[i] = geno1
        offspring[i+1] = geno2
    end
    return offspring
end

function runGeneration(population::Array{Individual,1}, config::ExperimentConfig, se_left::Matrix{Float64}, se_right::Matrix{Float64})
    parents::Array{Individual, 1} = rankingSelection_roulette(population, config.offspring_size)
    offspringGeno::Array{BitArray, 1} = createOffspringGenes(parents, config)
    for geno in offspringGeno
        mutate_bitflip_repair!(geno, config.dataset, config.remainingPoints)
        mutate_shift_harmonic_remaining_points!(geno, config.dataset, config.remainingPoints, 0.8)
    end
    offspring::Array{Individual, 1} = []
    if config.loop
        offspring = map(geno->Individual(geno, config.dataset, se_left, se_right), offspringGeno)
    else
        offspring = map(geno->Individual(geno, config.dataset), offspringGeno)
    end
    return genitor(population, offspring)
end

function runExperiment(config::ExperimentConfig, info=false, info_prefix="")
    remainingPoints::Int64 = config.remainingPoints
    ts::Timeseries = config.dataset
    se_left = Matrix{Float64}(undef, 3, length(ts))
    se_right = Matrix{Float64}(undef, 3, length(ts))
    df = createdf_run()
    totalTime::Float64 = 0.0
    population::Array{Individual,1} = []
    result = @timed begin
        population = createPopulation(ts, config.pop_size - length(config.specialStartingIndivuals), remainingPoints)
        for compression_algo in config.specialStartingIndivuals
            push!(population, Individual(compression_algo(ts, remainingPoints)[2], ts))
        end
    end
    totalTime += result[2]
    log_population!(df, population, 0, remainingPoints, totalTime)
    for i in 1:config.iterations
        result = @timed begin
            population = runGeneration(population, config, se_left, se_right)
        end
        totalTime += result[2]
        log_population!(df, population, i, remainingPoints, totalTime)
        info && println(info_prefix, "gen ", i, "/", config.iterations, " completed")
    end
    return df
end

function run_and_log(config::ExperimentConfig)
    file_meta = dir_experiment * "overview.csv"
    experiment_name = "experiment" * Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS-sss")
    file_experiment = dir_experiment * experiment_name * ".csv"
    df_experiment = createdf_total()
    df_config = createdf_config()
    dfs_run::Array{DataFrame, 1} = [createdf_total() for i in 1:config.runs]
    
    if config.timed
        for run in 1:config.runs
            println("starting run ", run, " at ", Dates.now())
            GC.gc()
            df_run = runExperiment(config, false, "run $(run)/$(config.runs) ")
            df_run.run_id = [run for i in 1:size(df_run)[1]]
            dfs_run[run] = df_run
            println("finished run at ", Dates.now())
        end
    else
        nextJobIdx = Atomic{Int}(1)
        outputFileMutex = Mutex()
        _everythread() do
            _runExperimentThread(config, config.runs, nextJobIdx, dfs_run, outputFileMutex)
        end
    end
    println(experiment_name, " completed at ", Dates.now())
    for df in dfs_run
        CSV.write(file_experiment, df, append=isfile(file_experiment))
    end
    push!(df_config, [
            experiment_name, config.dataset_name, config.consideredLength, config.remainingPoints,
            config.pop_size, config.offspring_size,
            config.loop, config.iterations,
            config.runs, mean(df-> minimum(df[:best_valid_mse]), dfs_run), 
            description(config)
        ])
    CSV.write(file_meta, df_config, append=isfile(file_meta))
end

function _runExperimentThread(config, numberOfRuns, nextJobIdx, outputArray, outputMutex)
    while true
        jobIdx = atomic_add!(nextJobIdx, Int(1))

        if jobIdx > numberOfRuns
            break
        end

        df_run = runExperiment(config, false)
        df_run.run_id = [jobIdx for i in 1:size(df_run)[1]]

        lock(outputMutex)
        outputArray[jobIdx] = df_run
        println("run $(jobIdx) completed")
        unlock(outputMutex)
    end
end


function _everythread(fun)
    ccall(:jl_threading_run, Ref{Cvoid}, (Any,), fun)
end


function run_experiments()
    #run functions once to force the compilation
    problems = [("Phoneme", 20000, 501), ("Phoneme", 25000, 626)]
    produced_offspring = 150000
    runs = 20
    timed = true
    appliesLoop = true
    pop_size = 200
    start_pop = []#[visvalingam_squarederror, swing_filter_number, ramer_douglas_peucker]
    configs::Array{ExperimentConfig, 1} = []
    for problem in problems
        push!(configs, ExperimentConfig(
            runs,
            pop_size,
            pop_size,
            ceil(Int, produced_offspring / pop_size),
            problem[1],
            datasetbyname(problem[1])[1:problem[2]],
            problem[2],
            problem[3],
            start_pop,
            appliesLoop,
            timed
        ))
    end
    for config in configs
        run_and_log(config)
        println("completed config $(config.dataset_name)[1:$(config.consideredLength)]")
    end
end

println("don't forget to 'set JULIA_NUM_THREADS=4'")
println("you currently have $(Threads.nthreads()) threads")
println("or dedicate a single core at real-time priority for profiling")