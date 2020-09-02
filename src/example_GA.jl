using Random
using Statistics
using PyPlot

const Individual = Tuple{Float64, Float64}


function fitness(geno::Individual)
    x = geno[1]
    y = geno[2]
    return 0.5*abs(x*y) - x + 0.5y
end

function crossover_onepoint(geno1::Individual, geno2::Individual)
    return (geno1[1], geno2[2]), (geno2[1], geno1[2])
end

function mutate_bitflip!(geno::Individual)
    flipProbability = 1/length(geno)
    draw = rand(length(geno))
    for i in 1:length(geno)
        if draw[i] <= flipProbability
            geno[i] = !geno[i]
        end
    end
    return geno
end

function mutate_normal_distribution(geno::Individual)
    offset = randn() * 0.2
    if rand() < 0.5
        return (max(-5, min(5, geno[1] + offset)), geno[2])
    else
        return (geno[1], max(-5, min(5, geno[2] + offset)))
    end
end

"""receives the weights of the individuals and returns the indexes of the n selected individuals"""
function roulette_wheel_sampling(weights::Array{Float64, 1}, n::Integer)
    totalWeights = sum(weights)
    accumulatedWeights = accumulate(+, weights)
    matingPool = Array{Int64, 1}(undef, n)
    for i in 1:n
        r = rand() * totalWeights
        matingPool[i] = searchsortedfirst(accumulatedWeights, r)
    end
    return matingPool
end

function tournament_selection_2(population::Array{Individual, 1}, λ::Integer)
    parents::Array{Individual, 1} = Array{Individual, 1}(undef, λ)
    for i in 1:λ
        cand1 = population[rand(1:length(population))]
        cand2 = population[rand(1:length(population))]
        if fitness(cand1) > fitness(cand2)
            parents[i] = cand1
        else
            parents[i] = cand2
        end
    end
    return parents
end

function genitor(population::Array{Individual, 1}, offspring::Array{Individual, 1})
    contestants = vcat(population, offspring)
    sort!(contestants, by=ind->fitness(ind), rev=true)
    return contestants[1:length(population)]
end

function create_population_bitarray(μ, len)
    return [bitrand(len) for i in 1:μ]
end

using Printf

function create_population_pairs(μ::Integer)
    return [(rand()*10-5, rand()*10-5) for i in 1:μ]
end

function run_generation(population::Array{Individual, 1}, λ::Integer)
    parents = tournament_selection_2(population, λ)
    offspring = Array{Individual, 1}(undef, λ)
    for i in 1:2:λ
        if i == λ       #uneven number of offspring
            offspring[i] = parents[i]
        else
            off1, off2 = crossover_onepoint(parents[i], parents[i+1])
            offspring[i] = off1
            offspring[i+1] = off2
        end
    end
    offspring = mutate_normal_distribution.(offspring)
    return genitor(population, offspring)
end

function print_stats(population::Array{Individual, 1}, gen::Integer)
    fitnesses = map(fitness, population)
    bestFitness, bestIndex = findmax(fitnesses)
    print(gen, " & ")
    for ind in population
        @printf "\$%0.2f" ind[1] 
        print(", ")
        @printf "%0.2f\$" ind[2] 
        print(" & ")
    end
    @printf "\$%0.2f\$" bestFitness
    println(" \\\\")
    for ind in population
        print(" & ")
        @printf "\$(%0.2f)\$" fitness(ind) 
    end
    println(" & \\\\")
    println("\\hline")
    return bestFitness
end

function run_genetic_algorithm(iterations::Integer, μ::Integer)
    log = Array{Float64, 1}(undef, iterations+1)
    population = create_population_pairs(μ)
    log[1] = print_stats(population, 0)
    for i in 1:iterations
        #println("iteration $i")
        population = run_generation(population, floor(Int, μ/2)*2)
        log[i+1] = print_stats(population, i)
    end
    return log
end

function average_and_plot(runs::Integer, iterations::Integer, μ::Integer)
    logs::Matrix{Float64} = Matrix{Float64}(undef, iterations+1, runs)
    for i in 1:runs
        logs[:,i] = run_genetic_algorithm(iterations, μ)
    end
    plot(0:iterations, mean(logs, dims=2))
    ylabel("average highest fitness")
    xlabel("generation")
    title("average fitness of best individual over $runs runs")
    grid()
end