##' @title  A wrapper function of prior estimation and bayNorm function

function bayNorm(;Data, BETA_vec=nothing, Conditions = nothing,UMI_sffl = nothing,Prior_type = nothing,mode_version =false,mean_version=false,S = 20,FIX_MU = true,BB_SIZE_par = true, verbose = true)

    if (mode_version & mean_version)
        show("Only one of mode_version and mean_version should be specified to be TRUE, otherwise both should be set to FALSE so that 3D array normalized data will be returned.")
    end

    if (!mode_version & !mean_version)
        myFunc = Main_NB_Bay
    elseif (mode_version & !mean_version)
        myFunc = Main_mode_NB_Bay
    elseif (!mode_version & mean_version)
        myFunc =  Main_mean_NB_Bay
    end

    input_params=Dict("BETA_vec"=>BETA_vec,"Conditions"=>Conditions,"UMI_sffl"=>UMI_sffl,"Prior_type"=>Prior_type,"FIX_MU"=>FIX_MU,"BB_SIZE"=>BB_SIZE_par)


    if (isnothing(BETA_vec)) 
        BETA_vec = sum(Data,dims=1)./mean(sum(Data,dims=1)).*0.06
        BETA_vec[BETA_vec.<=0].=minimum(BETA_vec[BETA_vec.>0])
        BETA_vec[BETA_vec.>=1].=maximum(BETA_vec[BETA_vec.<1])
        BETA_vec=BETA_vec[1,:]
    end


    if (isnothing(Conditions))
        if (isnothing(UMI_sffl))
            Data_sr=Data
        else
            Data_sr = round.(Data./UMI_sffl)
            if (isnothing(BETA_vec)) 
                BETA_vec = sum(Data_sr,dims=1)./mean(sum(Data_sr,dims=1)).*0.06
                BETA_vec[BETA_vec.<=0].=minimum(BETA_vec[BETA_vec.>0])
                BETA_vec[BETA_vec.>=1].=maximum(BETA_vec[BETA_vec.<1])
                BETA_vec=BETA_vec[1,:]
            end

        end

        PRIORS = Prior_fun(Data=Data_sr, BETA_vec=BETA_vec, FIX_MU = FIX_MU,BB_SIZE_par = BB_SIZE_par, verbose = verbose) 

        MU_input = PRIORS.MME_MU

        if (BB_SIZE_par) 
            SIZE_input = PRIORS.MME_SIZE_adjust
        else 
            SIZE_input = PRIORS.MME_SIZE
        end

        Bay_out = myFunc(
            Data = Data_sr,
            BETA_vec = BETA_vec,
            input_size = SIZE_input,
            input_mu = MU_input, S = S)

        bayNorm_output=Dict("Bay_out"=>Bay_out,"PRIORS"=>PRIORS,"input_params"=>input_params)

    else
        # multiple groups
        if (size(Data)[2] != length(Conditions)) 
            show("Number of columns in expression matrix must match length of conditions vector!")
        end
        if (isnothing(Prior_type)) 
            show("Prior_type needs to be specified when Conditions are specified, now Prior_type is set to be LL")
            Prior_type = "LL"
        end
        # if (is(names(Conditions))) {
        #     names(Conditions) <- colnames(Data)
        # }

        #Conditions=repeat([1,2],inner=50)

        Levels = unique(Conditions)

        DataList=Dict()
        BETAList=Dict()
        for i in range(1, stop = length(Levels))
            DataList[Levels[i]]=Data[:, Conditions .== Levels[i]]
            BETAList[Levels[i]]=BETA_vec[Conditions .== Levels[i]]
        end

        
        if (isnothing(UMI_sffl)) 
            # UMI
            DataList_sr = copy(DataList)
        else 
            # non-UMI
            DataList_sr=Dict()
            for i in range(1, stop = length(Levels))
                DataList_sr[Levels[i]]=round.(DataList[Levels[i]]./UMI_sffl[i])
            end

            
            if (isnothing(BETA_vec)) 
                BETAList=Dict()
                for i in range(1, stop = length(Levels))
                    Betatemp = sum(DataList_sr[Levels[i]],dims=1)./mean(sum(DataList_sr[Levels[i]],dims=1)).*0.06
                    Betatemp[Betatemp .<=0] .=minimum(Betatemp[Betatemp .>0])
                    Betatemp[Betatemp .>=1] .=maximum(Betatemp[Betatemp .<1])
                    BETAList[Levels[i]]=Betatemp[1,:]

                end
            end
                
        end

        PRIORS_LIST=Dict()
        if (Prior_type == "LL") 
            for i in range(1, stop = length(Levels)) 
                PRIORS_LIST[Levels[i]] = Prior_fun(Data=DataList_sr[Levels[i]], BETA_vec=BETAList[Levels[i]], FIX_MU = FIX_MU,BB_SIZE_par = BB_SIZE_par, verbose = verbose)
            end
        elseif (Prior_type == "GG") 
            Prior_fun(Data=DataList_sr[Levels[i]], BETA_vec=BETAList[Levels[i]], FIX_MU = FIX_MU,BB_SIZE_par = BB_SIZE_par, verbose = verbose)

            
            dat_qq=copy(DataList_sr)
            for i in range(2, stop = length(Levels)) 
                dat_qq[Levels[1]]=hcat(dat_qq[Levels[1]], dat_qq[Levels[i]])
            end

            beta_qq=copy(BETAList)
            for i in range(2, stop = length(Levels)) 
                beta_qq[Levels[1]]=vcat(beta_qq[Levels[1]], beta_qq[Levels[i]])
            end

            PROPRS_TEMP = Prior_fun(Data = dat_qq[Levels[1]],BETA_vec =beta_qq[Levels[1]], FIX_MU = FIX_MU,BB_SIZE_par = BB_SIZE_par, verbose = verbose)
            for i in range(1, stop = length(Levels))
                PRIORS_LIST[Levels[i]] = PROPRS_TEMP

            end

        end

        Bay_out_list = Dict()
        for i in range(1, stop = length(Levels))

            MU_input = PRIORS_LIST[Levels[i]].MME_MU

            if (BB_SIZE_par) 
                SIZE_input = PRIORS_LIST[Levels[i]].MME_SIZE_adjust
            else 
                SIZE_input = PRIORS_LIST[Levels[i]].MME_SIZE
            end

            Bay_out_list[Levels[i]] = myFunc(
                Data = DataList_sr[Levels[i]],
                BETA_vec = BETAList[Levels[i]],
                input_size = SIZE_input,
                input_mu = MU_input, S = S)

        end


        bayNorm_output=Dict("Bay_out_list"=>Bay_out_list,"PRIORS_LIST"=>PRIORS_LIST,"input_params"=>input_params)

    end

    if (verbose) 
        show("bayNorm has completed!")
    end
    
    return(bayNorm_output)
end


# Conditions=repeat([1,2],inner=50)
# outt=bayNorm(Data=Data, BETA_vec=BETA_vec, Conditions = Conditions,UMI_sffl = nothing,Prior_type = nothing,mode_version =false,mean_version=true,S = 20,FIX_MU = true,BB_SIZE_par = true, verbose = true)


