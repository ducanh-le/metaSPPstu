include("./spp_dm1.jl")

# ---------------------------------------------------------------------------
# GRASP

function greedyRandomizedConstruction(C, A, a)
    x = zeros(Int64, length(C))
    S = 1:length(C)
    M = 1:size(A, 1)
    E = Array{Array{Int64,1},1}(undef, length(C))
    Einit = Array{Array{Int64,1},1}(undef, length(C))
    F = Array{Array{Int64,1},1}(undef, length(C))
    U = Array{Float64,1}(undef, length(C))     # profit
    firstloop = true
    while (length(S) != 0)
        fill!(E, Int64[])
        fill!(F, Int64[])
        fill!(U, -Inf)
        if length(S) == 1
            iSel = S[1]
        else
            #mettre à jour E et F
            for j in M
                current = intersect(findall(isequal(1), A[j, :]), S)
                for i in current
                    E[i] = union(E[i], setdiff(current, i))
                    F[i] = union(F[i], j)
                end
            end
            #calculer uMax, uMin et Ui
            uMax = -Inf
            uMin = Inf
            for i in S
                if length(E[i]) != 0
                    U[i] = C[i] - sum(C[v] for v in E[i])
                else
                    U[i] = C[i]
                end
                if U[i] < uMin
                    uMin = U[i]
                end
                if U[i] > uMax
                    uMax = U[i]
                end
            end
            #calculer limit, RCL et choisir iSel aléatoire
            limit = uMin + a*(uMax - uMin)
            RCL = findall(>=(limit),U)
            iSel = RCL[rand(1:length(RCL))]
        end
        if firstloop
            Einit = copy(E)
            firstloop = false
        end
        #mettre a jour le problème
        x[iSel] = 1
        S = setdiff(S, union(E[iSel], iSel))
        M = setdiff(M, F[iSel])
    end
    z = sum(C[i]x[i] for i in 1:length(C))
    return x, z, Einit
end

function GRASP(C, A, a, finish)
    zMax = 0
    xMax = zeros(Int64, length(C))
    t = 0
    while (t <= finish)
        t = t + @elapsed begin
            xInit, zInit, Einit = greedyRandomizedConstruction(C, A, a)
            x, z = GreedyImprovement(C, A, xInit, zInit, Einit)
        end
        if z > zMax
            zMax = z
            xMax = x
        end
    end
    return xMax, zMax
end

# ---------------------------------------------------------------------------
# ReactiveGRASP

# Fonction tirer aléatoire un élément suivant son poids
function draw(a,p)
    k = findfirst(cumsum(p) .> rand())
    return k
end

function reactiveGRASP(C, A, a, finish, N)
    zMin = 0
    m = length(a)
    zTotal = Array{Pair{Int64, Int64}}(undef, m)        #pair la sommes de z[k] avec nb fois utilise a[k]
    fill!(zTotal, 0 => 0)
    zMax = 0
    xMax = zeros(Int64, length(C))
    p = Array{Float64,1}(undef, m)
    fill!(p, 1/m)
    q = Array{Float64,1}(undef, m)
    nbIter = 1
    t = 0
    while (t <= finish)
        #modifier p
        if (nbIter - 1)%N == 0 && nbIter != 1
            if zMax == zMin             #optimal found
                break
            else
                sumQ = 0
                for k in 1:m
                    q[k] = ((zTotal[k].first / zTotal[k].second) - zMin)/(zMax - zMin)
                    sumQ = sumQ + q[k]
                end
                for k in 1:m
                    p[k] = q[k]/sumQ
                end
            end
        end
        #tirer a[k]
        if nbIter <= m  #assurer que chaque a[k] a utilisé au moins une fois
            k = nbIter
        else
            k = draw(a,p)
        end
        #GRASP
        t = t + @elapsed begin
            xInit, zInit, Einit = greedyRandomizedConstruction(C, A, a[k])
            x, z = GreedyImprovement(C, A, xInit, zInit, Einit)
        end
        #mettre à jours les valeurs de x, z et nbIter
        if nbIter == 1
            zMax = z
            xMax = x
            zMin = z
        elseif z > zMax
            zMax = z
            xMax = x
        elseif z < zMin
            zMin = z
        end
        zTotal[k] = zTotal[k].first + z => zTotal[k].second + 1
        nbIter += 1
    end
    return xMax, zMax
end

# -------------------------------------------------------------------------------------------
# Version GRASP et ReactiveGRASP avec la condition d'arrêt est nombre des itérations
# L'ultiliser pour le graphique

function graspSPP(fname, a, nbIterationGrasp)
    target = "../../../Data"
    C, A = loadSPP(string(target,"/",fname))

    zconstruction = zeros(Int64,nbIterationGrasp)
    zamelioration = zeros(Int64,nbIterationGrasp)
    zbest = zeros(Int64,nbIterationGrasp)
    zbetter = 0

    for i=1:nbIterationGrasp
        xInit, zInit, Einit = greedyRandomizedConstruction(C, A, a)
        x, z = GreedyImprovement(C, A, xInit, zInit, Einit)
        zconstruction[i] = zInit
        zamelioration[i] = z
        zbetter = max(zbetter, zamelioration[i])
        zbest[i] = zbetter
    end
    return zconstruction, zamelioration, zbest
end

function regraspSPP(fname, a, nbIterationGrasp, N)
    target = "../../../Data"
    C, A = loadSPP(string(target,"/",fname))

    zconstruction = zeros(Int64,nbIterationGrasp)
    zamelioration = zeros(Int64,nbIterationGrasp)
    zbest = zeros(Int64,nbIterationGrasp)
    zbetter = 0

    zMin = 0
    m = length(a)
    zTotal = Array{Pair{Int64, Int64}}(undef, m)        #pair la sommes de z[k] avec nb fois utilise a[k]
    fill!(zTotal, 0 => 0)
    zMax = 0
    xMax = zeros(Int64, length(C))
    p = Array{Float64,1}(undef, m)
    fill!(p, 1/m)
    q = Array{Float64,1}(undef, m)
    nbIter = 1
    while (nbIter <= nbIterationGrasp)
        #modifier p
        if (nbIter - 1)%N == 0 && nbIter != 1
            if zMax == zMin             #optimal found
                break
            else
                sumQ = 0
                for k in 1:m
                    q[k] = ((zTotal[k].first / zTotal[k].second) - zMin)/(zMax - zMin)
                    sumQ = sumQ + q[k]
                end
                for k in 1:m
                    p[k] = q[k]/sumQ
                end
            end
        end
        #tirer a[k]
        if nbIter <= m  #assurer que chaque a[k] a utilisé au moins une fois
            k = nbIter
        else
            k = draw(a,p)
        end
        #GRASP
        xInit, zInit, Einit = greedyRandomizedConstruction(C, A, a[k])
        x, z = GreedyImprovement(C, A, xInit, zInit, Einit)
        zconstruction[nbIter] = zInit
        zamelioration[nbIter] = z
        zbetter = max(zbetter, zamelioration[nbIter])
        zbest[nbIter] = zbetter
        #mettre à jours les valeurs de x, z et nbIter
        if nbIter == 1
            zMax = z
            xMax = x
            zMin = z
        elseif z > zMax
            zMax = z
            xMax = x
        elseif z < zMin
            zMin = z
        end
        zTotal[k] = zTotal[k].first + z => zTotal[k].second + 1
        nbIter += 1
    end
    return zconstruction, zamelioration, zbest, p
end
