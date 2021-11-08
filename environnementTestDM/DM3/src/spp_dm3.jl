# ---------------------------------------------------------------------------
# Construction gloutonne d'une solution admissible de SPP - porter de DM1

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
# le metaheuristique Tabou

function Tabou(C,A,iterMax,iterPenaliser,length_tm,xInit,zInit,Einit)
    iter  = 1
    nIter = 1
    x = copy(xInit)
    z = zInit
    xBest = copy(x)
    zBest = z
    tb_list = Array{Tuple{Int64,Int64},1}(undef, length_tm)
    memory_longtime = zeros(Int64,2,length(C))          #1ere ligne est 0, 2eme ligne est 1
    while iter < iterMax
        eqZero = findall(isequal(0),x)
        eqOne  = findall(isequal(1),x)
        for j in eqOne
            setdiff!(eqZero,Einit[j])
        end
        changeable = union(eqZero,eqOne)    #table d'indice de x[j] qui peut être fait add ou drop (neighbor search)
        neighborBest = (0,-Inf)
        foundAdd = false                    #var boolean indiqué si on a trouvé un voisin valid avec mouvement add    
        for j in changeable
            if nIter < iterPenaliser        #pas de penalisation
                if x[j] == 0                #mouvement add
                    zPrime = z + C[j]
                    if (zPrime > zBest) || (!((j,1) in tb_list) && zPrime > neighborBest[2])           #override tabou or not tabou
                        neighborBest = (j,zPrime)
                        foundAdd = true
                    end
                elseif !foundAdd            #pas besoin de faire le mouvement drop si on a deja trouvé un voisin valid avec mouvement add en cas de sans penaliser
                    zPrime = z - C[j]
                    if zPrime >= neighborBest[2] && !((j,0) in tb_list)                                #not tabou, ajouter l'equal dans la comparaison pour le cas Σ(x_j) = 1
                        neighborBest = (j,zPrime)
                    end
                end
            else                            #commencer la penalisation pour chercher sur un autre voisinage
                if x[j] == 0
                    zPrime = z + C[j] - memory_longtime[1,j]
                    if (zPrime + memory_longtime[1,j] > zBest)                                         #override tabou and memory longtime
                        neighborBest = (j,zPrime + memory_longtime[1,j])
                    elseif (zPrime > neighborBest[2] && !((j,1) in tb_list))                           #not tabou
                        neighborBest = (j,zPrime)
                    end
                else
                    zPrime = z - C[j] - memory_longtime[2,j]
                    if zPrime >= neighborBest[2] && !((j,0) in tb_list)                                #not tabou, ajouter l'equal dans la comparaison pour le cas Σ(x_j) = 1
                        neighborBest = (j,zPrime)
                    end
                end
            end
        end
        #modifier x, z, taboulist, memoire longtemps, xBest, zBest, nIter, iter
        if neighborBest[1] == 0             #le cas de blocage car |TM| trop grande
            println("Blocked ! Every mouvement is TABOU ! Reduce |TM| by 1 ...")
            xBest2, zBest2, length_tm2 = Tabou(C, A, iterMax, iterMax, length_tm - 1, xInit, zInit, Einit)
            if zBest2 > zBest       #test pour choisir la meilleure solution entre avant et apres modifier |TM|
                return xBest2, zBest2, length_tm2
            else
                return xBest, zBest, length_tm
            end
        else
            before = x[neighborBest[1]]
            after  = (before == 0 ? 1 : 0)
            if (nIter >= iterPenaliser)
                if (neighborBest[2] < zBest)
                    neighborBest = (neighborBest[1], neighborBest[2] + memory_longtime[before+1,neighborBest[1]])        #retouner à la valeur sans penaliser, sauf si le statut tabou est dominé
                end
                nIter = 0       #reinitialiser nIter apres changer de voisinage
            end
            x[neighborBest[1]] = after
            z = neighborBest[2]
            tb_list = push!(tb_list[2:length_tm],(neighborBest[1],before))
            memory_longtime[before+1,neighborBest[1]] += 1
            if z > zBest
                zBest = z
                xBest = copy(x)
                nIter = 1       #reinitialiser nIter si trouver la nouvelle meilleure solution
                println("New zBest = ",z," found at iter = ",iter," !")
            else
                nIter += 1
            end
            iter += 1
        end
    end
    println("|TM| final = ", length_tm)
    return xBest, zBest, length_tm
end