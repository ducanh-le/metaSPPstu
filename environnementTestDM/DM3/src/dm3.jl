include("../../LibSPP/librarySPP.jl")
include("./spp_dm3.jl")
include("../../../setSPP.jl")


function main()
    println("Etudiant(e)s : LE et OUSMAN")

    # Collecting the names of instances to solve located in the folder Data ----
    target = "../dat"
    fnames = getfname(target)

    fres = splitdir(splitdir(pwd())[end-1])[end]
    io = open("../res/"*fres*".res", "w")
    for instance = 1:length(fnames)

        C, A = loadSPP(string(target,"/",fnames[instance]))

        iterMax,iterPenaliser,length_tb_list = 70000,45,8

        xBest, zBest = Tabou(C,A,iterMax,iterPenaliser,length_tb_list)
        isAdmissible(C,A,xBest)
        
        println(io, "xBest = ",xBest)
        println(io, "zBest = ",zBest)
    end
    close(io)

end

main()
