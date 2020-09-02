using Test
include("../timeseries.jl")
include("../ramer_douglas_peucker.jl")
include("../genetic.jl")

ts1 = Timeseries(collect(enumerate([0, 1, 4, 3.5, 3.25, 2.625, 3.75, 2, 1, -2, -1.5, -0.5, 0])))

@test ramer_douglas_peucker(ts1, 2)[1] == Timeseries([(1, 0), (length(ts1), 0)])
@test ramer_douglas_peucker(ts1, length(ts1))[1] == ts1
@test ramer_douglas_peucker(ts1, 7)[1] == Timeseries([(1, 0), (3, 4), (6, 2.625), (7, 3.75), (9, 1), (10, -2), (length(ts1), 0)])
for i in 2:length(ts1)
    pheno, geno = ramer_douglas_peucker(ts1, i)
    toGenotype(pheno, ts1) == geno
end

println("tested ramer_douglas_peucker")