function rouletteWheelSampling(weights::Array{Float64, 1}, n::Integer)
    totalWeights = sum(weights)
    accumulatedWeights = accumulate(+, weights)
    matingPool = Array{Int64, 1}(undef, n)
    for i in 1:n
        r = rand() * totalWeights
        matingPool[i] = searchsortedfirst(accumulatedWeights, r)
    end
    return matingPool
end

function stochasticUniversalSampling(weights::Array{Float64, 1}, n::Integer)
    #would be best to create a random permutation before calling this method
    totalWeights = sum(weights)
    accumulatedWeights = accumulate(+, weights)
    matingPool = zeros(Int, n)
    current = 1
    r = rand() * totalWeights / length(weights)
    i = 1
    while i <= n
        while r <= accumulatedWeights[current] && i <= n
            matingPool[i] = current
            i += 1
            r += totalWeights / length(weights)
        end
        current += 1
    end
    return matingPool
end

function uniformSampling(weights::Array{Float64, 1}, n::Integer)
    return [floor(Int, rand() * length(weights)) + 1 for i in 1:n]
end

#don't use roulette wheel but use Stochastic Universal Sampling
function selectProportionalToFitness(fitnesses)
    total_fitness = sum(fitnesses)
    value = rand() * total_fitness
    i = 1
    while fitnesses[i] < value
        value -= fitnesses[i]
        i += 1
    end
    return i
end

function rankingSelection(population::Array{Individual, 1}, n::Integer, sampling_algorithm::Function)
    sortedPopulation = sort(population, by=x -> x.fitness)
    weights = [Float64(i) for i in 1:length(population)]
    selection = sampling_algorithm(weights, n)
    return map(i -> sortedPopulation[i], selection)
end

function rankingSelection_roulette(population::Array{Individual, 1}, n::Integer)
    return rankingSelection(population, n, rouletteWheelSampling)
end

function rankingSelection_stochasticUniversal(population::Array{Individual, 1}, n::Integer)
    return rankingSelection(population, n, stochasticUniversalSampling)
end

function rankingSelection_uniform(population::Array{Individual, 1}, n::Integer)
    return rankingSelection(population, n, uniformSampling)
end

function tournamentSelection(population::Array{Individual, 1}, n::Integer, tournamentSize::Integer=1)
    matingPool::Array{Individual, 1} = Array{Individual, 1}(undef, n)
    for i in 1:n
        contestant_indexes = rand(1:length(population), tournamentSize)
        winner_index = findmin(map(i -> (population[i].fitness, i), contestant_indexes))[1][2]
        matingPool[i] = population[winner_index]
    end
    return matingPool
end

function tournamentSelection10(population::Array{Individual, 1}, n::Integer)
    return tournamentSelection(population, n, 10)
end

function tournamentSelection5(population::Array{Individual, 1}, n::Integer)
    return tournamentSelection(population, n, 5)
end

function genitor(population::Array{Individual, 1}, offspring::Array{Individual, 1})
    return sort!(vcat(population, offspring), by=ind->ind.mse)[1:length(population)]
end

function bestHalvesOfBoth(population::Array{Individual, 1}, offspring::Array{Individual, 1})
    #assumes offspring to be long enough and population to be even-sized
    k = floor(Int, length(population)/2)
    contestants = vcat(population, offspring)
    return vcat(sort(population, by=x -> x.fitness)[1:k], sort(offspring, by=x -> x.fitness)[1:k])
end

function replaceAll(population::Array{Individual, 1}, offspring::Array{Individual, 1})
    if length(population) != length(offspring)
        @warn "replaceAll has not been made for this configuration"
    end
    return offspring
end

function oneNewVsOneOld(population::Array{Individual, 1}, offspring::Array{Individual, 1})
    survivors = Array{Individual, 1}(undef, length(population))
    shuffle!(offspring)
    for i in 1:length(population)
        if i > length(offspring) || population[i].fitness < offspring[i].fitness
            survivors[i] = population[i]
        else
            survivors[i] = offspring[i]
        end
    end
    return survivors
end