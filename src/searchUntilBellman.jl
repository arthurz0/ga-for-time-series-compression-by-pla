using CSV
include("experiment.jl")

function findFirstIndividual(config::ExperimentConfig, target_mse::Float64)
    remainingPoints::Int64 = config.remainingPoints
    ts::Timeseries = config.dataset
    se_left = Matrix{Float64}(undef, 3, length(ts))
    se_right = Matrix{Float64}(undef, 3, length(ts))
    population::Array{Individual,1} = []
    population = createPopulation(ts, config.pop_size - length(config.specialStartingIndivuals), remainingPoints)
    for compression_algo in config.specialStartingIndivuals
        push!(population, Individual(compression_algo(ts, remainingPoints)[2], ts))
    end
    for i in 1:config.iterations
        population = runGeneration(population, config, se_left, se_right)
        best_valid_ind = Individual(trues(0), [(0.0, 0.0), (1.0, 1.0)], Inf, Inf)
        for ind in population
            if length(ind.pheno) <= remainingPoints && ind.mse < best_valid_ind.mse
                best_valid_ind = ind
            end
        end
        if best_valid_ind.mse <= target_mse
            return best_valid_ind
        end
    end
    return nothing
end

function firstIndividualBeatingBellman(dataset, consideredLength, remainingPoints)
    bellman_df = CSV.read("experiments2/bellman_notable.csv")
    relevantEntries = filter(row -> row[:dataset] == dataset && row[:length] == consideredLength && row[:remaining_points] == remainingPoints, bellman_df)
    target_mse = relevantEntries[:mse][1]
    config = ExperimentConfig(
        1,
        200,
        200,
        500,
        dataset,
        datasetbyname(dataset)[1:consideredLength],
        consideredLength,
        remainingPoints,
        [visvalingam_squarederror, swing_filter_number, ramer_douglas_peucker],
        true,
        false
    )
    ind = findFirstIndividual(config, target_mse)
    return ind
end