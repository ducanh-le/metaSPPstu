function f1(N::Int)
    x = Array{Int64,1}(undef, N)
    for n = 1:N
        x[n] = n
    end
    return(x)
end

function f2(N::Int)
    x = Array{Int64,1}()
    for n = 1:N
        push!(x, n)
    end
    return(x)
end

N = 100000000
@time f1(N)
@time f2(N)