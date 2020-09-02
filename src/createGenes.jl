using Random
using Statistics

include("loop.jl")

function createGenesWithNPoints(bitarray_length, m, n)
    return [shuffle!(vcat(trues(n), falses(bitarray_length - n))) for i in 1:m]
end

function createPopulation(original::Timeseries, m, remainingPoints::Integer)
    genes = createGenesWithNPoints(length(original) - 2, m, remainingPoints - 2)
    return map(g -> Individual(g, original), genes)
end 

function crossover_onepoint(geno1::BitArray, geno2::BitArray)
    j = floor(Int, rand() * length(geno1))
    return vcat(geno1[1:j], geno2[j+1:end]), vcat(geno2[1:j], geno1[j+1:end])
end

function crossover_onepoint_balanced(geno1::BitArray, geno2::BitArray)
    total1 = sum(geno1)
    total2 = sum(geno2)
    j = floor(Int, rand() * length(geno1))
    total1Left = sum(geno1[1:j])
    total2Left = sum(geno2[1:j])
    totalOff1 = total1Left + (total2-total2Left)
    totalOff2 = total2Left + (total1-total1Left)
    if rand() < 0.5
        while abs(totalOff1 - totalOff2) > 1 && j > 0                   #go left
            if geno1[j]
                totalOff1 -= 1
                totalOff2 += 1
            end
            if geno2[j]
                totalOff1 += 1
                totalOff2 -= 1
            end
            j -= 1
        end
    else
        while abs(totalOff1 - totalOff2) > 1 && j < length(geno1)       #go right
            j += 1
            if geno1[j]
                totalOff1 += 1
                totalOff2 -= 1
            end
            if geno2[j]
                totalOff1 -= 1
                totalOff2 += 1
            end
        end
    end
    return vcat(geno1[1:j], geno2[j+1:end]), vcat(geno2[1:j], geno1[j+1:end])
end

function mutation_control_id(population::Array{Individual, 1}, remainingPoints::Integer, oldParameter::Float64)
    return oldParameter
end

function mutation_linearly_increasing_001(population::Array{Individual, 1}, remainingPoints::Integer, oldParameter::Float64)
    return oldParameter + 0.001
end

function mutation_exponentially_decreasing_999(population::Array{Individual, 1}, remainingPoints::Integer, oldParameter::Float64)
    return oldParameter * 0.999
end

function mutation_exponentially_decreasing_9995(population::Array{Individual, 1}, remainingPoints::Integer, oldParameter::Float64)
    return oldParameter * 0.9995
end

function mutation_control_exponential_05_rechenberg(population::Array{Individual, 1}, newIndsNumber::Integer, oldParameter::Float64)
    return mutation_control_exponential(population, newIndsNumber, oldParameter, 0.05, 0.2)
end

function mutation_control_exponential(population::Array{Individual, 1}, newIndsNumber::Integer, oldParameter::Float64, α::Float64, threshold::Float64)
    return newIndsNumber / length(population) > threshold ? oldParameter * (1.0+α) : oldParameter / (1.0+α)
end

function mutate_bitflip!(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer)
    compressionRatio = Float64(remainingPoints) / length(geno)
    n = length(geno)
    zeroToOne = compressionRatio / n
    oneToZero = (1 - compressionRatio) / n
    draw = rand(n)
    for i in 1:n
        if geno[i]
            if draw[i] <= oneToZero
                geno[i] = false
            end
        else
            if draw[i] <= zeroToOne
                geno[i] = true
            end
        end
    end
    return geno
end

function mutate_bitflip_adaptive!(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer, baseNumber::Float64=1.0)
    truesN = sum(geno)
    falsesN = length(geno) - truesN
    truesOptimal = remainingPoints - 2
    falsesOptimal = length(geno) - truesOptimal
    zeroToOne = truesN < truesOptimal ? Float64(truesOptimal - truesN + baseNumber) / falsesN : baseNumber / falsesN
    oneToZero = falsesN < falsesOptimal ? Float64(falsesOptimal - falsesN + baseNumber) / truesN : baseNumber / truesN
    draw = rand(length(geno))
    for i in 1:length(geno)
        if geno[i]
            if draw[i] <= oneToZero
                geno[i] = false
            end
        else
            if draw[i] <= zeroToOne
                geno[i] = true
            end
        end
    end
    return geno
end

function mutate_bitflip_adaptive_noparam!(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer, baseNumber::Float64=1.0)
	baseNumber = 0.0
    truesN = sum(geno)
    falsesN = length(geno) - truesN
    truesOptimal = remainingPoints - 2
    falsesOptimal = length(geno) - truesOptimal
    zeroToOne = truesN < truesOptimal ? Float64(truesOptimal - truesN + baseNumber) / falsesN : baseNumber / falsesN
    oneToZero = falsesN < falsesOptimal ? Float64(falsesOptimal - falsesN + baseNumber) / truesN : baseNumber / truesN
    draw = rand(length(geno))
    for i in 1:length(geno)
        if geno[i]
            if draw[i] <= oneToZero
                geno[i] = false
            end
        else
            if draw[i] <= zeroToOne
                geno[i] = true
            end
        end
    end
    return geno
end

function mutate_bitflip_repair!(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer)
    truesN = sum(geno)
    falsesN = length(geno) - truesN
    truesOptimal = remainingPoints - 2
    falsesOptimal = length(geno) - truesOptimal
    if truesN == truesOptimal
        return geno
    end
    tooManyOnes = truesN > truesOptimal
    remainToChange = tooManyOnes ? truesN - truesOptimal : falsesN - falsesOptimal
    remainOverall = tooManyOnes ? truesN : falsesN
    draw = rand(length(geno))
    for i in 1:length(geno)
        if geno[i] == tooManyOnes
            if draw[i] <= remainToChange/remainOverall
                geno[i] = !tooManyOnes
                remainToChange -= 1
                if remainToChange == 0
                    break
                end
            end
            remainOverall -= 1
        end
    end
    return geno
end

function mutate_shift_harmonic_remaining_points!(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer, expectedSwaps::Float64=1.0)
    return mutate_shift!(geno, original, remainingPoints, expectedSwaps, harmonic_offset_remainingPoints)
end

function mutate_shift_harmonic!(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer, expectedSwaps::Float64=1.0)
    return mutate_shift!(geno, original, remainingPoints, expectedSwaps, harmonic_offset)
end

function mutate_shift_normal!(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer, expectedSwaps::Float64=1.0)
    return mutate_shift!(geno, original, remainingPoints, expectedSwaps, normal_offset)
end

function mutate_shift!(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer, expectedSwaps::Float64, offsetFunction::Function)
    numberOfPoints = sum(geno)
    changeProbability = expectedSwaps / numberOfPoints
    index = 1
    while true
        index = findnext(identity, geno, index+1)
        if index == nothing
            break
        end
        r = rand()
        if r < changeProbability
            offset = offsetFunction(geno, original, remainingPoints)
            j = min(length(geno), max(1, index + offset))
            geno[j] ⊻= geno[index]
            geno[index] ⊻= geno[j]
        end
    end
    return geno
end

function harmonic_offset(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer)
    sign = rand() < 0.5 ? 1 : -1
    offset = floor(Int, ℯ^(rand()*log(length(geno))))
    return sign*offset
end

function harmonic_offset_remainingPoints(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer)
    sign = rand() < 0.5 ? 1 : -1
    offset = floor(Int, ℯ^(rand()*log(length(geno) / remainingPoints)))
    return sign*offset
end

function normal_offset(geno::BitArray{1}, original::Timeseries, remainingPoints::Integer)
    return trunc(Int, randn() * length(original)/remainingPoints/4)
end

function average_genotype(population::Array{Individual, 1})
    return sum(ind -> ind.geno, population) / length(population)
end

function fitness_std(population::Array{Individual, 1})
    return std(map(ind -> ind.fitness, population))
end

function genotype_diversity(population::Array{Individual, 1})
    return sum(x -> min(x, 1-x), average_genotype(population)) / length(population)
end