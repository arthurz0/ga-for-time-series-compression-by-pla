using Test
using Random
using Statistics

include("../timeseries.jl")
include("../genetic.jl")
include("../createGenes.jl")

for i in 1:100
    g1 = bitrand(10)
    g2 = bitrand(10)
    off1, off2 = crossover_onepoint_balanced(g1, g2)
    @test abs(sum(off1) - sum(off2)) < 2 || (g1 == off1 && g2==off2 || g1 ==off2 && g2==off1)
    @test sum(off1) + sum(off2) == sum(g1) + sum(g2)
end

for i in 1:100
    g1 = bitrand(5)
    g2 = bitrand(5)
    off1, off2 = crossover_onepoint_balanced(g1, g2)
    @test abs(sum(off1) - sum(off2)) < 2 || (g1 == off1 && g2==off2 || g1 ==off2 && g2==off1)
    @test sum(off1) + sum(off2) == sum(g1) + sum(g2)
end

for i in 1:100
    g1 = bitrand(100)
    g2 = bitrand(100)
    off1, off2 = crossover_onepoint_balanced(g1, g2)
    @test abs(sum(off1) - sum(off2)) < 2 || (g1 == off1 && g2==off2 || g1 ==off2 && g2==off1)
    @test sum(off1) + sum(off2) == sum(g1) + sum(g2)
end

function testAdaptiveBitflip(genolength, remainingPoints)
    unused_ts = Timeseries(collect(zip(1:genolength, 1:genolength)))
    resultingNumberOfPoints = [0 for i in 1:1000]
    for i in 1:1000
        geno = bitrand(genolength)
        offspring = mutate_bitflip_adaptive!(geno, unused_ts, remainingPoints)
        resultingNumberOfPoints[i] = sum(offspring)
    end
    @test abs(mean(resultingNumberOfPoints) - (remainingPoints - 2)) < 0.5
end

testAdaptiveBitflip(100, 20)
testAdaptiveBitflip(100, 50)
testAdaptiveBitflip(100, 10)
testAdaptiveBitflip(100, 75)

for i in 1:1000
    geno = bitrand(100)
    mutate_bitflip_repair!(geno, [(0.0, 0.0)], 20)
    @test sum(geno) == 18
end

println("createdGenes tested")