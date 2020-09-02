using Random
include("genetic.jl")
include("ramer_douglas_peucker.jl")
include("visvalingam.jl")
include("swingfilter.jl")
include("loop.jl")
include("datasets.jl")

function improvePoints(original::Timeseries, compressed::Timeseries)
    improved::Timeseries = copy(compressed)
    leftIndex = 1
    rightIndex = 1
    middleIndex = findwhile(p -> timeof(p) <= timeof(compressed, 2), original, 1)
    for i in 2:length(compressed)- 1
        rightIndex = findwhile(p -> timeof(p) <= timeof(compressed, i+1), original, middleIndex+1)
        currentIndex = findwhile(p -> timeof(p) <= timeof(compressed, i), original, middleIndex+1)
        bestIndex, bestY = bestindexandvalue(improved[i-1], improved[i+1], original[leftIndex+1:rightIndex-1], currentIndex-leftIndex:currentIndex-leftIndex)
        bestIndex += leftIndex
        improved[i] = (timeof(original, bestIndex), bestY)
        leftIndex = bestIndex
        middleIndex = rightIndex
    end
    rightIndex = findwhile(p -> timeof(p) <= timeof(improved, 2), original)
    m, b = bestline_fixpoint(original[1:rightIndex-1], improved[2])
    improved[1] = (timeof(improved, 1), m*timeof(improved, 1) + b)
    m, b = bestline_fixpoint(original[leftIndex+1:end], improved[end-1])
    improved[end] = (timeof(improved, length(improved)), m*timeof(improved, length(improved)) + b)
    return improved
end

function challengeOpys(ts::Timeseries, remainingPoints::Integer)
    geno = shuffle!(vcat(trues(remainingPoints - 2), falses(length(ts) - remainingPoints)))
    disconnected = false
    for i in 1:length(geno)-1
        if geno[i] && geno[i+1]
            disconnected = true
        end
    end
    @show remainingPoints
    @show disconnected
    ind = Individual(geno, ts)
    @show ind.mse
    improved = improvePoints(ts, ind.pheno)
    improvedMSE = mse(improved, ts)
    @show improvedMSE
    @show ind.mse - improvedMSE
end

ts = loadUCR("Ham")#loadSubject103(5)#loadUCR("Rock")
for i in 10:50:1000
    for j in 1:5
        challengeOpys(ts, i)
    end
end