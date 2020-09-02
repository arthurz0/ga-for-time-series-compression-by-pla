using Test

include("../swingfilter.jl")
include("../timeseries.jl")
include("../genetic.jl")
include("../datasets.jl")

function Base.:≈(a::Timeseries, b::Timeseries)
    @assert length(a) == length(b)
    return all(i -> a[i][1] ≈ b[i][1] && a[i][2] ≈ b[i][2], 1:length(a))
end

ts1 = Timeseries(collect(enumerate([0, 1.1, 1.9, 3.1, 5.0])))
@test epsilon_sampling(ts1, 0.1)[1] == ts1
@test epsilon_sampling(ts1, 0.2)[1] == [(1.0, 0.0), (4.0, 3.1), (5.0, 5.0)]
@test epsilon_sampling(ts1, 0.1)[2] == toGenotype(epsilon_sampling(ts1, 0.1)[1], ts1)
@test epsilon_sampling(ts1, 0.2)[2] == toGenotype(epsilon_sampling(ts1, 0.2)[1], ts1)

@test swing_filter(ts1, 0.1)[1] == [(1, 0), (4, 3), (5, 5)]
@test swing_filter(ts1, 0.2)[1] ≈ [(1.0, 0.0), (4.0, 3.042857143), (5.0, 5.0)]
@test swing_filter(ts1, 0.1)[2] == toGenotype(swing_filter(ts1, 0.1)[1], ts1)
@test swing_filter(ts1, 0.2)[2] == toGenotype(swing_filter(ts1, 0.2)[1], ts1)

@test swing_filter_number(ts1, 5)[1] == ts1
@test swing_filter_number(ts1, 3)[1] ≈ [(1.0, 0.0), (4.0, 3.042857143), (5.0, 5.0)]
@test length(swing_filter_number(loadApple(), 30)[1]) == 30
println("tested epsilon_sampling")