using Test
include("../timeseries.jl")
include("../optimalYs.jl")

ts1 = Timeseries([(1.0, 1.0), (4.0, 4.0)])
ts2 = Timeseries([(1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (4.0, 4.0)])
times1 = [1.0, 4.0]

@test optimalYsBlock(ts1, times1) == [1.0, 4.0]
@test optimalYsBlock(ts2, times1) ≈ [1.0, 4.0]

println("tested OpYs. but there aren't many tests...")