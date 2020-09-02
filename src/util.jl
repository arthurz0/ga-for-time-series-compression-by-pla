"""
Returns index of last subsequent element fulfilling the predicate. If element i does not fulfill predicate, i-1 will be returned.
Only works for one-dimensional Arrays.
"""
function findwhile(predicate::Function, A, i::Integer = 1)
    #@warn "findwhile is apparently quite inefficient"
    index = findnext(element->!predicate(element), A, i)
    if index == nothing
        return length(A)
    else
        return index - 1
    end
end

"""
A wrapper for sum to deal with empty itr. This method has 0.0 as default value.
"""
function sum_float(f, itr)
    if length(itr) > 0
        return sum(f, itr)
    else
        return 0.0
    end
end