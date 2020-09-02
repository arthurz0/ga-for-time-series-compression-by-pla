using Test

include("../loop.jl")
include("../loop_deprecated.jl")
include("../datasets.jl")
include("../swingfilter.jl")
include("../genetic.jl")

function Base.:≈(a::Timeseries, b::Timeseries)
    @assert length(a) == length(b)
    return all(i -> a[i][1] ≈ b[i][1] && a[i][2] ≈ b[i][2], 1:length(a))
end

ts1_o = Timeseries(collect(zip(0.0:4.0, [0.0, 1.1, 1.0, 1.0, 0.0])))
ts1_c = Timeseries([(0.0, 0.0), (1.0, 1.0), (4.0, 0.0)])
@test loop(ts1_c, ts1_o)[1] ≈ [(0.0, 0.0), (1.0, 1.3499999999999999), (4.0, 0.2571428571428571)]
@test loop(ts1_c, ts1_o)[2] == toGenotype(loop(ts1_c, ts1_o)[1], ts1_o)

ts2_o = Timeseries(collect(zip(1.0:10.0, [1.0, 2.0, 3.0, 3.5, 5.5, 6.0, 5.0, 4.0, 3.0, 2.0])))
ts2_c = Timeseries([(1.0, 1.0), (pi, pi), (10.0, 2.0)])
@test loop(ts2_c, ts2_o)[1] ≈ [(1.0, 0.9427198817442719), (6.0, 6.032520325203252), (10.0, 1.9891598915989164)]
@test loop(ts2_c, ts2_o)[2] == toGenotype(loop(ts2_c, ts2_o)[1], ts2_o)

ts1 = loadApple()
looped = epsilon_sampling(ts1, 20.0)[1]
looped_2 = epsilon_sampling(ts1, 20.0)[1]
for i in 1:5
    global looped
    global looped_2
    looped = loop_deprecated(ts1, looped)
    looped_2 = loop(looped_2, ts1)[1]
    @test looped_2 ≈ looped
end

println("tested loop")