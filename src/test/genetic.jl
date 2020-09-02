using Test

include("../genetic.jl")
include("../datasets.jl")
include("../ramer_douglas_peucker.jl")

ts1 = loadApple()
geno1 = ramer_douglas_peucker(ts1, 200)[2]
@test toGenotype(Individual(geno1, ts1).pheno, ts1) == geno1

println("tested genetic")