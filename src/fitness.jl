include("genetic.jl")

function fitnesscontrol_exponential_01_20(population::Array{Individual, 1}, remainingPoints::Integer, oldPenalty::Float64)
    return fitnesscontrol_exponential(population, remainingPoints, oldPenalty, 0.01, 0.2)
end

function fitnesscontrol_exponential(population::Array{Individual, 1}, remainingPoints::Integer, oldPenalty::Float64, α::Float64, threshold::Float64)
    invalidIndividuals = 0
    for ind in population
        if length(ind.pheno) > remainingPoints
            invalidIndividuals += 1
        end
    end
    return invalidIndividuals / length(population) > threshold ? oldPenalty * (1.0+α) : oldPenalty / (1.0+α)
end

function fitness_penalty_exponential!(offspring::Array{Individual, 1}, remainingPoints::Integer, base::Float64)
    for ind in offspring
        penaltyFactor = length(ind.pheno) > remainingPoints ? base ^ (length(ind.pheno) - remainingPoints) : 1.0
        ind.fitness = ind.mse * penaltyFactor
    end
end

function fitness_penalty_linear!(offspring::Array{Individual, 1}, remainingPoints::Integer, coefficient::Float64)
    for ind in offspring
        penaltyFactor = length(ind.pheno) > remainingPoints ? 1.0 + coefficient * (length(ind.pheno) - remainingPoints) : 1.0
        ind.fitness = ind.mse * penaltyFactor
    end
end