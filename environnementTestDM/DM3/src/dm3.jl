include("../../LibSPP/librarySPP.jl")
include("./spp_dm3.jl")
include("../../../setSPP.jl")


function main()
    println("Etudiant(e)s : LE et OUSMAN")

    # Collecting the names of instances to solve located in the folder Data ----
    target = "../dat"
    fnames = getfname(target)
    println()

    fres = splitdir(splitdir(pwd())[end-1])[end]
    io = open("../res/"*fres*".res", "w")
    for instance = 1:length(fnames)

        C, A = loadSPP(string(target,"/",fnames[instance]))

        iterMax,iterPenaliser,length_tm = 100000,rand(35:45),rand(5:9)

        println("Instance : ",fnames[instance])
        println("iterMax = ", iterMax, " | iterPenaliser = ", iterPenaliser, " | |TM| init = ", length_tm)
        sleep(5)        #pause 5s pour lire les informations

        t1 = @elapsed xInit, zInit, Einit = GreedyConstruction(C,A)
        println("z(x0) = ",zInit)
        t2 = @elapsed xBest, zBest, length_tm2 = Tabou(C,A,iterMax,iterPenaliser,length_tm,xInit,zInit,Einit)

        println(io, "Instance : ",fnames[instance])
        println(io, "iterMax = ", iterMax, " | iterPenaliser = ", iterPenaliser, " | |TM| init = ", length_tm, " | |TM| final = ", length_tm2)
        if isAdmissible(C,A,xBest)
            println(io, "Feasible : yes | Σ(x_j) = ", length(findall(isequal(1),xBest)), " | z(x0) = ", zInit, " | z(x) = ", zBest)
        else
            println(io, "Feasible : no")
        end
        println(io, "CPU time = ",t1 + t2)
        println(io)
        println()
    end
    close(io)

end

main()
