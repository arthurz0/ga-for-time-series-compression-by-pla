using CSV
using Dates
include("timeseries.jl")

data_dir = normpath(@__DIR__, "../data/")

function datasetbyname(name)
    if name == "apple"
        return loadApple()
    elseif name == "subject103_5"
        return loadSubject103(5)
    elseif name == "subject103_6"
        return loadSubject103(6)
    elseif name == "Ham"
        return loadUCR("Ham")
    elseif name == "Rock"
        return loadUCR("Rock")
    elseif name == "MixedShapesRegularTrain"
        return loadUCR("MixedShapesRegularTrain")
    elseif name == "HandOutlines"
        return loadUCR("HandOutlines")
    elseif name == "StarLightCurves"
        return loadUCR("StarLightCurves")
    elseif name == "Phoneme"
        return loadUCR("Phoneme")
    elseif name == "UWaveGestureLibraryAll"
        return loadUCR("UWaveGestureLibraryAll")
    elseif name == "UWaveGestureLibraryAll"
        return loadUCR("UWaveGestureLibraryAll")
    elseif name == "MALLAT_"
        return loadTSSA("MALLAT_")
    elseif name == "Donoho-Johnstone"
        return loadTSSA("Donoho-Johnstone")
    end
    @warn "requested data set not available"
    return nothing
end

function paperEvalData()
    return [("Donoho-Johnstone", 2048, 52), ("MALLAT_", 8192, 206), ("subject103_5", 10000, 251), ("subject103_6", 10000, 251), ("Ham", 10000, 251), ("Rock", 10000, 251)]
end

function loadApple()
    data = CSV.read(data_dir * "aapl.csv")
    times_data = map(d -> Dates.value(d), data[:Date])
    startTime = min(times_data...)
    times_data = map(t -> t - startTime, times_data)
    ts::Timeseries = sort(collect(zip(times_data, data[:Open])))
    return ts
end

function loadSubject103(i)
    #good i: 3, 5:Inf
    data = CSV.read(data_dir * "subject103.csv")
    times_data = data[1]
    startTime = min(times_data...)
    times_data = map(t -> t - startTime, times_data)
    ts::Timeseries = sort(collect(zip(times_data, data[i])))
    ts = filter(p -> !isnan(valueof(p)), ts)
    return ts
end

function loadTSSA(name)
    data = CSV.read(data_dir * name * ".txt", delim="\t", header=false)
    ts::Timeseries = sort(collect(zip(data[:Column1], data[:Column2])))
    return ts
end

function loadUCR(name)
    data1 = CSV.read(data_dir * name * "_TRAIN.tsv", delim="\t", header=false)
    data2 = CSV.read(data_dir * name * "_TEST.tsv", delim="\t", header=false)
    data = vcat(data1, data2)
    tssVals::Array{Array{Float64, 1}, 1} = []
    for i in 1:size(data, 1)
        push!(tssVals, data[i,2:end])
    end
    return joinTimeseries(map(vals -> Timeseries(collect(enumerate(vals))), tssVals))
end

function joinTimeseries(tss::Array{Timeseries, 1})
    result = tss[1]
    for ts in tss[2:end]
        dt = timeof(ts[1]) - timeof(result[end])
        dv = valueof(ts[1]) - valueof(result[end])
        append!(result, map(p -> (timeof(p) - dt, valueof(p) - dv), ts[2:end]))
    end
    return result
end