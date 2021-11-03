#= ATTENTION:
   your own folder is considered as the current working directory
   for running your solver
=#

include("../../LibSPP/librarySPP.jl")
include("./spp_dm1.jl")
include("../../../setSPP.jl")


function main()
    println("Etudiant(e)s : LE et OUSMAN")

    # Collecting the names of instances to solve located in the folder Data ----
    target = "../../Data"
    fnames = getfname(target)

    fres = splitdir(splitdir(pwd())[end-1])[end]
    io = open("../res/"*fres*".res", "w")
    for instance = 1:length(fnames)

        # Load one numerical instance ------------------------------------------
        C, A = loadSPP(string(target,"/",fnames[instance]))

        t1 = @elapsed x, zInit, Einit = GreedyConstruction(C, A)
        t2 = @elapsed xBest, zBest = GreedyImprovement(C, A, x, zInit, Einit)

        # Saving results -------------------------------------------------------
        println(io, fnames[instance], " ", zInit, " ", zBest, " ", t1, " ", t2, " ", t1+t2)
    end
    close(io)

end

main()
