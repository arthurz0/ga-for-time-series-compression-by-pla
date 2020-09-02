using CSV

struct ExperimentConfig
    runs::Integer
    pop_size::Integer
    offspring_size::Integer
    iterations::Integer
    dataset_name::String
    dataset::Timeseries
    consideredLength::Integer
    remainingPoints::Integer
    specialStartingIndivuals::Array{Function, 1}
    loop::Bool
    timed::Bool
end

function description(config::ExperimentConfig)
    desc = "version: 02.07.2019 "
    if isempty(config.specialStartingIndivuals)
        desc *= " the initial population contains no special Individuals"
    else
        desc *= " the initial population contains special Individuals:"
        for f in config.specialStartingIndivuals
            desc *= " $(f) "
        end
    end
    return desc
end

function createdf_config()
    return DataFrame(
        name = String[],
        dataset = String[],
        consideredLength = Int[],
        remaining_points = Int[],
        pop_size = Int[],
        offspring_size = Int[],
        loop = Bool[],
        iterations = Int[],
        runs = Int[],
        mean_best_mse = Float64[],
        description = String[]    
    )
end