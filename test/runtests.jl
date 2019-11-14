using bayNormJL
using Test

@testset "Getopt" begin

    Daa=[1:10 1:10 3:12 4:13] 
    Bea=[0.1,0.5,0.3,0.4]
    outt=bayNormJL.bayNorm(Data=Daa, BETA_vec=Bea, Conditions = nothing,UMI_sffl = nothing,Prior_type = nothing,mode_version =false,mean_version=true,S = 20,FIX_MU = true,BB_SIZE_par = true, verbose = true)

	@test size(outt["Bay_out"]) ==(10,4)

end



# add Distributed
# add ProgressMeter
# add DataFrames
# add GLM
# add SharedArrays
# add Distributions
# add Random
# add Optim
# add SpecialFunctions
# add BlackBoxOptim