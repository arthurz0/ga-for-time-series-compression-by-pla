using Test

include("../visvalingam.jl")
include("../genetic.jl")

ts1 = Timeseries(collect(enumerate([0, 1, 4, 3.5, 3.25, 2.625, 3.75, 2, 1, -2, -1.5, -0.5, 0])))
ts2 = Timeseries(collect(enumerate([0, 1, 2, 4, 5, 2, 5])))

@test visvalingam(ts2, 7, visvalingam_area_error_helper)[1] == [(1.0, 0), (2.0, 1), (3.0, 2), (4.0, 4), (5.0, 5), (6.0, 2), (7.0, 5)]
@test visvalingam(ts2, 6, visvalingam_area_error_helper)[1] == [(1.0, 0), (3.0, 2), (4.0, 4), (5.0, 5), (6.0, 2), (7.0, 5)]
@test visvalingam(ts2, 5, visvalingam_area_error_helper)[1] == [(1.0, 0), (3.0, 2), (5.0, 5), (6.0, 2), (7.0, 5)]
@test visvalingam(ts2, 4, visvalingam_area_error_helper)[1] == [(1.0, 0), (5.0, 5), (6.0, 2), (7.0, 5)]
@test visvalingam(ts2, 3, visvalingam_area_error_helper)[1] == [(1.0, 0), (5.0, 5), (7.0, 5)]
@test visvalingam(ts2, 2, visvalingam_area_error_helper)[1] == [(1.0, 0), (7.0, 5)]

@test visvalingam(ts2, 7, visvalingam_approximated_se_helper)[1] == [(1.0, 0), (2.0, 1), (3.0, 2), (4.0, 4), (5.0, 5), (6.0, 2), (7.0, 5)]
@test visvalingam(ts2, 6, visvalingam_approximated_se_helper)[1] == [(1.0, 0), (3.0, 2), (4.0, 4), (5.0, 5), (6.0, 2), (7.0, 5)]
@test visvalingam(ts2, 5, visvalingam_approximated_se_helper)[1] == [(1.0, 0), (3.0, 2), (5.0, 5), (6.0, 2), (7.0, 5)]
@test visvalingam(ts2, 4, visvalingam_approximated_se_helper)[1] == [(1.0, 0), (5.0, 5), (6.0, 2), (7.0, 5)]
@test visvalingam(ts2, 3, visvalingam_approximated_se_helper)[1] == [(1.0, 0), (5.0, 5), (7.0, 5)]
@test visvalingam(ts2, 2, visvalingam_approximated_se_helper)[1] == [(1.0, 0), (7.0, 5)]

for i in 2:13
    @test visvalingam(ts1, i, visvalingam_area_error_helper)[2] == toGenotype(visvalingam(ts1, i, visvalingam_area_error_helper)[1], ts1)
end

@test visvalingam_approximated_se_helper((0.0, 0.0), (3.0, 4.0), (5.0, 1.0)) ≈ 19.266666666666666
@test visvalingam_area_error_helper((0.0, 0.0), (3.0, 4.0), (5.0, 1.0)) == 5.0 * 3.4 / 2 

println("tested visvalingam")