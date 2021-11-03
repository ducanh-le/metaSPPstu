using JuMP, GLPK

# ---------------------------------------------------------------------------
# modelisation SPP
# T la matrice (nbVar x nbSC) de coefficients de variables dans les contraintes, soit 0 soit 1
# c le vecteur de coefficients de variables dans la fonction objective

function model_spp(C, A)
    m::Model = Model(GLPK.Optimizer)
    set_optimizer_attribute(m, "msg_lev", GLPK.GLP_MSG_ALL)

    # nombre de variable
    nbVar::Int64 = length(C)
    # nombre de contrainte
    nbSC::Int64 = size(A,1)

    # declaration des variables
    @variable(m, x[1:nbVar], Bin)

    # declaration fonction objective
    @objective(m, Max, sum(C[j]x[j] for j = 1:nbVar))

    # declaration contraintes
    @constraint(m, sc[i = 1:nbSC], sum(A[i,j]x[j] for j in 1:nbVar) <= 1)

    return m
end

# ---------------------------------------------------------------------------
# Construction gloutonne d'une solution admissible de SPP

function GreedyConstruction(C, A)
    x = zeros(Int64, length(C))
    S = 1:length(C)
    M = 1:size(A, 1)
    E = Array{Array{Int64,1},1}(undef, length(C))
    Einit = Array{Array{Int64,1},1}(undef, length(C))
    F = Array{Array{Int64,1},1}(undef, length(C))
    Cout = Array{Float64,1}(undef, length(C))
    firstloop = true
    while (length(S) != 0)
        fill!(E, Int64[])
        fill!(F, Int64[])
        fill!(Cout, -Inf)
        if length(S) == 1
            imax = S[1]
        else
            for j in M
                current = intersect(findall(isequal(1), A[j, :]), S)
                for i in current
                    E[i] = union(E[i], setdiff(current, i))
                    F[i] = union(F[i], j)
                end
            end
            for i in S
                if length(E[i]) != 0
                    Cout[i] = C[i] - sum(C[u] for u in E[i])
                end
            end
            if maximum(Cout) == -Inf        #pas de conflit
                imax = S[1]
            else
                imax = findfirst(isequal(maximum(Cout)), Cout)
            end
        end
        if firstloop
            Einit = copy(E)
            firstloop = false
        end
        x[imax] = 1
        S = setdiff(S, union(E[imax], imax))
        M = setdiff(M, F[imax])
    end
    z = sum(C[i]x[i] for i in 1:length(C))
    return x, z, Einit
end

# ---------------------------------------------------------------------------
# Amelioration gloutonne par recherche locale d’une solution de SPP

function findAllCombination(S, E)
    common = Int64[]
    while S != Int64[]
        if (intersect(S, E[S[1]]) != [])
            S1 = setdiff(S,E[S[1]])
            S2 = setdiff(S,S[1])
            res1 = findAllCombination(S1,E)
            for i in 1:length(res1)
                res1[i] = union(res1[i], common)
            end
            res2 = findAllCombination(S2,E)
            for i in 1:length(res2)
                res2[i] = union(res2[i], common)
            end
            return union(res1,res2)
        else
            common = union(common, S[1])
            S = setdiff(S, S[1])
        end
    end
    return [common]
end

function GreedyImprovement(C, A, xInit, zInit, Einit)
    meilleur = true
    xMax = xInit
    zMax = zInit
    while meilleur
        meilleur = false
        Y = findall(isequal(1), xMax)
        zCurrent = zMax
        xCurrent = copy(xMax)
        for i in Y
            Ek = Int64[]
            for k in setdiff(Y,i)
                Ek = union(Ek,Einit[k])
            end
            S = setdiff(Einit[i],Ek)
            Combis = findAllCombination(S,Einit)
            if Combis != [Int64[]]
                for u in 1:length(Combis)
                    Cv = sum(C[v] for v in Combis[u])
                    Ci = C[i]
                    zPrime = zCurrent - Ci + Cv
                    if zPrime > zMax
                        meilleur = true
                        zMax = zPrime
                        xMax = copy(xCurrent)
                        for v in Combis[u]
                            xMax[v] = 1
                        end
                        xMax[i] = 0
                    end
                end
            end
        end
    end
    return xMax, zMax
end
