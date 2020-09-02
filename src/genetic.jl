include("optimalYs.jl")
include("timeseries.jl")

struct Individual
    geno::BitArray{1}
    pheno::Timeseries
    mse::Float64
    fitness::Float64
    Individual(g::BitArray{1}, original::Timeseries) = length(g) + 2 == length(original) ? _build(g, original) : error("length of genotype and original does not match")
    Individual(g::BitArray{1}, original::Timeseries, se_left::Matrix{Float64}, se_right::Matrix{Float64}) = length(g) + 2 == length(original) ? _buildLoop(g, original, se_left, se_right) : error("length of genotype and original does not match")
    Individual(g::BitArray{1}, p::Timeseries, mse::Float64, f::Float64) = sum(g) + 2 == length(p) ? new(g, p, mse, f) : error("length of phenotype does not match genotype")
end

function _build(geno::BitArray{1}, original::Timeseries)
    times = vcat(timeof(original[1]), map(i -> timeof(original, i), filter(i -> geno[i-1], 2:length(geno)+1)), timeof(original[end]))
    pheno = optimalValues(original, times)
    calculated_mse = mse(pheno, original)
    return Individual(geno, pheno, calculated_mse, Inf)
end

function _buildLoop(geno::BitArray{1}, original::Timeseries, se_left::Matrix{Float64}, se_right::Matrix{Float64})
    times::Array{Float64, 1} = vcat(timeof(original[1]), map(i -> timeof(original, i), filter(i -> geno[i-1], 2:length(geno)+1)), timeof(original[end]))
    pheno::Timeseries = optimalValues(original, times)
    loop!(pheno, original, se_left, se_right)
    pheno = optimalValues(original, timesof(pheno))
    calculated_mse = mse(pheno, original)
    return Individual(geno, pheno, calculated_mse, Inf)
end

function best_maxpoints(population::Array{Individual, 1}, k::Integer)
    best = findfirst(ind -> length(ind.pheno) <= k, population)
    for p in population
        if length(p.pheno) <= k && p.mse < best.mse
            best = p
        end
    end
    return best
end

function toGenotype(compressed::Timeseries, original::Timeseries)
    @assert timeof(original[1]) == timeof(compressed[1])
    @assert timeof(original[end]) == timeof(compressed[end])
    index = 2
    geno = falses(length(original)-2)
    for p in compressed[2:end-1]
        index = findwhile(p_original -> timeof(p_original) <= timeof(p), original, index)
        if timeof(original[index]) != timeof(p)
            @warn "times do not match exactly"
        end
        geno[index - 1] = true
    end
    @assert sum(geno) + 2 == length(compressed)
    return geno
end