include("../../LibSPP/librarySPP.jl")
include("./spp_dm2.jl")
include("../../../setSPP.jl")


function main()
    println("Etudiant(e)s : LE et OUSMAN")

    # Collecting the names of instances to solve located in the folder Data ----
    target = "../../Data"
    fnames = getfname(target)

    fres = splitdir(splitdir(pwd())[end-1])[end]
    io = open("../res/"*fres*".res", "w")
    for instance = 1:length(fnames)

        C, A = loadSPP(string(target,"/",fnames[instance]))

        resGRASP   = Array{Int64}(undef, 10)
        resReGRASP = Array{Int64}(undef, 10)

        for iter in 1:10
            t1 = @elapsed xGRASP, zGRASP = GRASP(C, A, 0.7, 10)
            t2 = @elapsed xReGRASP, zReGRASP = reactiveGRASP(C, A, [0.2,0.5,0.7,0.8,0.9,1], 120, 17)
            resGRASP[iter]   = zGRASP
            resReGRASP[iter] = zReGRASP
            #setSPP(C,A)
            println(io, fnames[instance], " | iter = ", iter, " | zGRASP = ", zGRASP, " | zReGRASP = ", zReGRASP, " | t1 = ", t1, " | t2 = ", t2)
        end
        println(io, "GRASP         : zMin = ", minimum(resGRASP), " | zMax = ", maximum(resGRASP), " | zMoy = ", sum(resGRASP[k] for k in 1:10)/10)
        println(io, "reactiveGRASP : zMin = ", minimum(resReGRASP), " | zMax = ", maximum(resReGRASP), " | zMoy = ", sum(resReGRASP[k] for k in 1:10)/10)
        println(io)
    end
    close(io)

end

main()
