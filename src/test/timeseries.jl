using Test
include("../timeseries.jl")

ts1 = Timeseries([(1, 2), (3, 6), (3.5, -99), (7, 38)])
@test timesof(ts1) == [1, 3, 3.5, 7]
@test valuesof(ts1) == [2, 6, -99, 38]
@test all(i -> timeof(ts1, i) == ts1[i][1], 1:length(ts1))
@test all(i -> valueof(ts1, i) == ts1[i][2], 1:length(ts1))

p = Point((3.5, -9.25))
@test timeof(p) == p[1]
@test valueof(p) == p[2]

@test mse(Timeseries([(0.0, 5)]), Timeseries([(-4, 3), (6.3, 5), (7, 4), (8, 5.5), (10, 7)])) == 9.25/5
@test mse(ts1, ts1) == 0
@test mse(Timeseries([(3, 2), (6, 8), (8, 6)]), Timeseries(collect(zip(1.0:10.0, [3, 0, 2.5, 4.25, 6, 7.5, 7, 6, 5.5, 6])))) == 5.8125 / 10

println("tested timeseries")