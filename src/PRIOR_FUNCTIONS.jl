

function EstPrior(;Data, verbose = true) 

    #add progress bar
    
    # p = Progress(size(Data)[1], dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
    # channel = RemoteChannel(()->Channel{Bool}(size(Data)[1]), 1)


    temp_result=SharedArray{Float64,2}((size(Data)[1],2))


    # @sync begin
    #     # this task prints the progress bar
    #     @async while take!(channel)
    #         next!(p)
    #     end
    #     # this task does the computation
    #     @async begin
            @distributed for gene_ind=1:size(Data)[1]
            #println(gene_ind)
            vals=Data[gene_ind,:]
            vals=convert(Vector,vals)
            n = length(vals)
            m = mean(vals)
            v = (n - 1)/n*var(vals)
            #sleep(3)
            if (v > m) 
                mme_size=m^2/(v - m)
            else 
                mme_size=NaN
            end
            estimate=[mme_size, m]
            temp_result[gene_ind,:]=estimate
            #put!(channel, true)
        #end
        #     put!(channel, false) # this tells the printing task to finish
            end

        if verbose
            show("Priors estimation based on MME method has completed.")
        end
        
    #end

    return(temp_result)
end

function AdjustSIZE_fun(;BB_SIZE, MME_MU, MME_SIZE)

    
    fitind = findall((BB_SIZE .< maximum(BB_SIZE)*0.8) .* (BB_SIZE .> minimum(BB_SIZE)))


    #lmfit = lm(log.(BB_SIZE)[fitind] ~ log.(MME_SIZE)[fitind])
    Y=log.(BB_SIZE)[fitind]
    X=log.(MME_SIZE)[fitind]
    Z=reshape(X, length(X), 1)
    data = DataFrame(X=X, Y=Y)
    
    lmfit=lm(@formula(Y ~ X), data)
    #lmfit=lm(Z, Y)


    MME_SIZE_adjust = coef(lmfit)[1] .+ coef(lmfit)[2] .* log.(MME_SIZE)
    #MME_SIZE_adjust = coef(lmfit)[1] .* log.(MME_SIZE)
    MME_SIZE_adjust = exp.(MME_SIZE_adjust)
    return(MME_SIZE_adjust)
end

function t2(d1,d2)
    d1=[d1 d2]
    d1
end

function BB_Fun(;Data, BETA_vec, INITIAL_MU_vec, INITIAL_SIZE_vec, MU_lower = 0.01, MU_upper = 1000, SIZE_lower = 0.01,SIZE_upper = 100,FIX_MU = true)

    Priors=[INITIAL_SIZE_vec INITIAL_MU_vec]
    Priors=convert(SharedArray,Priors)
    Data=convert(SharedArray,Data)

    # p = Progress(size(Data)[1], dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
    # channel = RemoteChannel(()->Channel{Bool}(size(Data)[1]), 1)
    nCol=size(Data)[1]


    temp_result=SharedArray{Float64,2}((size(Data)[1],2))
    #sum(temp_result[:,1].==0)
    if FIX_MU
        # @sync begin    
        #     # this task prints the progress bar
        #     @async while take!(channel)
        #         next!(p)
        #     end
        #     # this task does the computation
        #     @async begin
               @distributed (t2) for gene_ind in range(1, stop = nCol)
                    println(gene_ind)
                    #gene_ind=4086
                    input_obs=Data[gene_ind,:]

                    input_mme=Priors[gene_ind,2]
                    input_mme_size=Priors[gene_ind,1]

                    function temp_fun(sss)
                        temp_val=MarginalF_NB_1D(SIZE=sss[1],MU=input_mme, m_observed=input_obs,BETA=BETA_vec)
                        return(temp_val)
                    end
                    
                    function temp_fun_g!(G,sss)        
                        G[1]=GradientFun_NB_1D(SIZE=sss[1],MU=input_mme, m_observed=input_obs,BETA=BETA_vec)
                    end
                    #dd=Uniform(SIZE_lower, SIZE_upper)

                    #random initial parameters
                    # r_r=rand(dd,10)
                    # storrr=SharedArray{Float64,1}(10)
                    # for i = 1:length(r_r)
                    #     inner_optimizer =  LBFGS()
                    #     res_1 =optimize(temp_fun,temp_fun_g!,[SIZE_lower],[SIZE_upper],[r_r[i]], Fminbox(inner_optimizer))
                    #     storrr[i]=res_1.minimizer[1]
                    # end


                    # function fg!(x, g)
                    #     g = GradientFun_NB_1D(SIZE=x[1],MU=input_mme, m_observed=input_obs,BETA=BETA_vec)
                    #     return MarginalF_NB_1D(SIZE=x[1],MU=input_mme, m_observed=input_obs,BETA=BETA_vec)
                    #   end

                    # function projection_of(x)
                    #     if x<SIZE_lower
                    #         x=SIZE_lower
                    #     end
                    #     if x>SIZE_upper
                    #         x=SIZE_upper
                    #     end
                    # end

                    # function prj!(xp, x)
                    # xp = projection_of(x)
                    # end
                    # x0=DenseArray{Float64,1}
                    # ws = spg2(fg!, prj!,x0 , 1)

                    # inner_optimizer =  LBFGS()
                    # res_m =optimize(temp_fun,temp_fun_g!,[SIZE_lower],[SIZE_upper],[input_mme_size], Fminbox(inner_optimizer))

                    res = bboptimize(temp_fun; SearchRange = [(SIZE_lower, SIZE_upper)], NumDimensions = 1,Method= :xnes,MaxFuncEvals=1000,TraceMode=:silent)


                    #mean([SIZE_upper,SIZE_lower])
                    #pmap(x -> QQ(x):x,r_r)
              
                    #rr_val=res_m.minimizer[1]
                    rr_val=best_candidate(res)[1]
                    
                    #res =optimize(SIZE -> MarginalF_NB_1D(SIZE,MU=input_mme, m_observed=input_obs,BETA=BETA_vec),SIZE -> GradientFun_NB_1D(SIZE,MU=input_mme, m_observed=input_obs,BETA=BETA_vec),[0.01],[100.], [Priors[gene_ind,1]], Fminbox(inner_optimizer))
                    
                    temp_result[gene_ind,:]=[rr_val,Priors[gene_ind,2]]
                    
                    #put!(channel, true)
                end
        #         put!(channel, false) # this tells the printing task to finish
        #     end

        # end 

    else  #######################################2D ##############
        # @sync begin
        #     # this task prints the progress bar
        #     @async while take!(channel)
        #         next!(p)
        #     end
        #     # this task does the computation
        #     @async begin
                @distributed (t2) for gene_ind in range(1, stop = nCol)
                    println(gene_ind)
                    #gene_ind=4086
                    input_obs=Data[gene_ind,:]

                    input_mme=Priors[gene_ind,:]
                    if input_mme[2]!=0
                

                        function temp_fun(sss)
                            temp_val=MarginalF_NB_2D(SIZE_MU=sss, m_observed=input_obs,BETA=BETA_vec)
                            return(temp_val)
                        end
                        
                        function temp_fun_g!(G,sss)   
                            teee=GradientFun_NB_2D(SIZE_MU=sss, m_observed=input_obs,BETA=BETA_vec)
                            G[1]=teee[1]
                            G[2]=teee[2]
                        end

                        #dd_1=Uniform(SIZE_lower, SIZE_upper)
                        #dd_2=Uniform(MU_lower, MU_upper)

                        # inner_optimizer =  LBFGS()
                        # res =optimize(temp_fun,temp_fun_g!,[SIZE_lower,MU_lower],[SIZE_upper,MU_upper], input_mme, Fminbox(inner_optimizer))

                        res = bboptimize(temp_fun; SearchRange = [(SIZE_lower, SIZE_upper),(MU_lower, MU_upper)], NumDimensions = 2,Method= :xnes,MaxFuncEvals=1000,TraceMode=:silent)
                        #res =optimize(temp_fun,temp_fun_g!,[SIZE_lower,MU_lower],[SIZE_upper,MU_upper], [rand(dd_1,1)[1],rand(dd_2,1)[1]], Fminbox(inner_optimizer))

                        
                        #temp_result[gene_ind,:]=res.minimizer
                        temp_result[gene_ind,:]=best_candidate(res)
                    else
                        temp_result[gene_ind,:]=input_mme

                    end
            #         put!(channel, true)
            #     end
            #     put!(channel, false) # this tells the printing task to finish
            # end

        end 
    end #end of parallel for loop

    return(temp_result)

end


#FIX_MU = false is not suggested, use FIX_MU = true instead.
function Prior_fun(;Data, BETA_vec, FIX_MU = true,BB_SIZE_par = true, verbose = true) 


    normcount_N = transpose( (1 ./BETA_vec) .*  transpose(Data))
    Priors_MME = EstPrior(Data=normcount_N, verbose = true)
    Priors_MME[isnan.(Priors_MME[:,1]),1].=minimum(Priors_MME[.!isnan.(Priors_MME[:,1]),1])


    M_ave_ori = Priors_MME[:,2]
    size_est = Priors_MME[:,1]

    # BB_SIZE_par = true
    # FIX_MU = false
    # verbose = true

    if BB_SIZE_par
        if verbose
            show("Start optimization using spg from BB package. This part may be time-consuming.")
        end
        ooot=BB_Fun(Data=Data, BETA_vec=BETA_vec, INITIAL_MU_vec=Priors_MME[:,2], INITIAL_SIZE_vec=Priors_MME[:,1], MU_lower = minimum(M_ave_ori), MU_upper = maximum(M_ave_ori), SIZE_lower =minimum(size_est),SIZE_upper = maximum(size_est),FIX_MU = FIX_MU)

    end 
    #sum(ooot[:,1].<=0)
    ooot[ooot[:,1].<=0,1].=minimum(ooot[ooot[:,1].>0,1])
    MME_SIZE_adjust=AdjustSIZE_fun(BB_SIZE=ooot[:,1], MME_MU=M_ave_ori, MME_SIZE=size_est)
    if BB_SIZE_par
        if FIX_MU
            output=DataFrame(MME_MU=M_ave_ori,MME_SIZE=size_est,BB_SIZE=ooot[:,1],MME_SIZE_adjust=MME_SIZE_adjust)
        else
            output=DataFrame(MME_MU=M_ave_ori,MME_SIZE=size_est,BB_SIZE=ooot[:,1],BB_MU=ooot[:,2],MME_SIZE_adjust=MME_SIZE_adjust)
        end
    else
        output=DataFrame(MME_MU=M_ave_ori,MME_SIZE=size_est)
    end
    #CSV.write("D:/RNAseqProject/bayNorm_dev/bayNorm_Julia/temp_result.csv",output)
    return(output)
end



#oot=Prior_fun(Data=Data, BETA_vec=BETA_vec, FIX_MU = true,BB_SIZE_par = true, verbose = true)









#script############

# AdjustSIZE_fun(BB_SIZE, MME_MU, MME_SIZE) 
# AdjustSIZE_fun(BB_SIZE=BB_SIZE, MME_MU=MME_MU, MME_SIZE=MME_SIZE) 


#G=Array{Float64,1}(undef, 1)

#res = bboptimize(temp_fun; SearchRange = [(0.001, 100.)], NumDimensions = 1,Method= :xnes,MaxFuncEvals=1000,TraceMode=:silent)
#inner_optimizer = GradientDescent()


