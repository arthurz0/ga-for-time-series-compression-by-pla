include("timeseries.jl")
include("loop.jl")

function calculate_error_edges(ts::Timeseries)
    numberVertices = length(ts)
    edgeArray::Matrix{Float64} = Matrix{Float64}(undef, numberVertices, numberVertices)
    reuseStorage::Matrix{Float64} = Matrix{Float64}(undef, 3, numberVertices)
    fill!(edgeArray, Inf)
    for i in 1:numberVertices-1
        polynomials::Matrix{Float64} = mse_polynomials_by_slope(ts[i+1:end], ts[i], reuseStorage)
        for j in i+1:numberVertices
            polynomial1 = polynomials[1, j-i]
            polynomial2 = polynomials[2, j-i]
            polynomial3 = polynomials[3, j-i]
            s = slope(ts[i], ts[j])
            edgeArray[j, i] = polynomial1 * s^2 + polynomial2 * s + polynomial3
        end
    end
    return edgeArray
end

function bellman(ts::Timeseries, remainingPoints::Integer, returnmse=false)
    numberVertices = length(ts)
    edgeArray = calculate_error_edges(ts)
    weights = [Inf for n in 1:numberVertices]
    weights[1] = 0
    iterations = remainingPoints - 1
    predecessor = Matrix{Int32}(undef, numberVertices, iterations)
    fill!(predecessor, 0)
    for it in 1:iterations
        if it > 1
            predecessor[:, it] = predecessor[:, it-1]
        end
        for p in numberVertices:-1:1
            for pointAhead in p:numberVertices
                errorOverEdge = weights[p] + edgeArray[pointAhead, p]
                if errorOverEdge < weights[pointAhead]
                    weights[pointAhead] = errorOverEdge
                    predecessor[pointAhead, it] = p
                end
            end
        end
    end
    compressed::Timeseries = []
    geno = falses(length(ts)-2)
    pre = length(ts)
    for i in iterations:-1:1
        push!(compressed, ts[pre])
        if 1 < pre < length(ts)
            geno[pre - 1] = true
        end
        pre = predecessor[pre, i]
        if pre == 0
            break
        end
    end
    push!(compressed, ts[1])
    reverse!(compressed)
    if returnmse
        return compressed, geno, weights[end]/length(ts)
    else
        return compressed, geno
    end
end

function bellman_notable(ts::Timeseries, remainingPoints::Integer, returnmse=false)
    numberVertices = length(ts)
    weights = [Inf for n in 1:numberVertices]
    weights[1] = 0
    iterations = remainingPoints - 1
    predecessor = Matrix{Int32}(undef, numberVertices, iterations)
    fill!(predecessor, 0)
    for it in 1:iterations
        if it > 1
            predecessor[:, it] = predecessor[:, it-1]
        end
        for p in numberVertices:-1:1
            coeff1::Float64 = 0.0
            coeff2::Float64 = 0.0
            coeff3::Float64 = 0.0
            from::Point = ts[p]
            for pointAhead in p+1:numberVertices
                dx = timeof(ts[pointAhead]) - timeof(from)
                dy = valueof(ts[pointAhead]) - valueof(from)
                coeff1 += dx^2  
                coeff2 += -2.0dy*dx   
                coeff3 += dy^2         
                s = slope(ts[p], ts[pointAhead])
                localError = coeff1 * s^2 + coeff2 * s + coeff3
                errorOverEdge = weights[p] + localError
                if errorOverEdge < weights[pointAhead]
                    weights[pointAhead] = errorOverEdge
                    predecessor[pointAhead, it] = p
                end
            end
        end
    end
    compressed::Timeseries = []
    geno = falses(length(ts)-2)
    pre = length(ts)
    for i in iterations:-1:1
        push!(compressed, ts[pre])
        if 1 < pre < length(ts)
            geno[pre - 1] = true
        end
        pre = predecessor[pre, i]
        if pre == 0
            break
        end
    end
    push!(compressed, ts[1])
    reverse!(compressed)
    if returnmse
        return compressed, geno, weights[end]/length(ts)
    else
        return compressed, geno
    end
end