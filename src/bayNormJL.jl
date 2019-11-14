module bayNormJL

    using Distributed
    using ProgressMeter
    using DataFrames
    using GLM
    using SharedArrays
    using Distributions
    using Random
    using Optim
    using SpecialFunctions
    using BlackBoxOptim
    #using OptimPack
    #using JuliaDB
    #using RData
    #using CSV

    include("bayNorm_main.jl")
    include("PRIOR_FUNCTIONS.jl")
    include("bayNorm_wrapper.jl")

    export 
    bayNorm, 
    AdjustSIZE_fun,
    EstPrior,
    BB_Fun,
    Prior_fun,
    DownSampling,
    GradientFun_NB_1D,
    GradientFun_NB_2D,
    MarginalF_NB_1D,
    MarginalF_NB_2D,
    Main_NB_Bay,
    Main_mean_NB_Bay,
    Main_mode_NB_Bay

end # module