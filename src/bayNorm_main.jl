##' @title Binomial downsampling

function DownSampling(;Data,BETA_vec)
    Counts_downsampling=Array{Int64,2}(undef, size(Data))

    for i in 1:size(Data)[1]
        for j in 1:size(Data)[2]
            d=Binomial(Int(Data[i,j]),BETA_vec[j])
            Counts_downsampling[i,j]= rand(d, 1)[1]
        end
    end

    return(Counts_downsampling)
end

#dd=DownSampling(;Data=RAW_DAT,BETA_vec=BETA_vec)
# m_observed=RAW_DAT[1,:]
# BETA=BETA_vec
# SIZE=1
# MU=3
# S=5

function GradientFun_NB_1D(;SIZE,MU, m_observed,BETA)
    m=m_observed
    nCells=length(m_observed)
    temp_vec_2_SIZE=Array{Float64,1}(undef, length(m_observed))

    # dd=NegativeBinomial(SIZE,SIZE/(SIZE+BETA[i]*MU))
    # temp_vec_2[i]=pdf(dd,m_observed[i])
    for i=1:nCells
        temp_vec_2_SIZE[i]=digamma(m_observed[i]+SIZE)-digamma(SIZE)+log(SIZE/(SIZE+MU*BETA[i]))+(BETA[i]*MU-m_observed[i])/(BETA[i]*MU+SIZE);
    end
    Gradd=sum(temp_vec_2_SIZE);
    return Gradd;
end




function GradientFun_NB_2D(;SIZE_MU, m_observed,BETA)
    m=m_observed
    nCells=length(m_observed)

    temp_mu=Array{Float64,1}(undef, length(m_observed))
    temp_size=Array{Float64,1}(undef, length(m_observed))
    Gradd_vec=Array{Float64,1}(undef, 2)

    # dd=NegativeBinomial(SIZE,SIZE/(SIZE+BETA[i]*MU))
    # temp_vec_2[i]=pdf(dd,m_observed[i])
    for i=1:nCells
        temp_mu[i]=(m_observed[i]*SIZE_MU[1]-SIZE_MU[2]*BETA[i]*SIZE_MU[1])/(SIZE_MU[2]*(SIZE_MU[2]*BETA[i]+SIZE_MU[1]));

        temp_size[i]=digamma(m_observed[i]+SIZE_MU[1])-digamma(SIZE_MU[1])+log(SIZE_MU[1]/(SIZE_MU[1]+SIZE_MU[2]*BETA[i]))+(BETA[i]*SIZE_MU[2]-m_observed[i])/(BETA[i]*SIZE_MU[2]+SIZE_MU[1]);

    end
    Gradd_vec[1]=sum(temp_size);
    Gradd_vec[2]=sum(temp_mu);

    return Gradd_vec
end
#GradientFun_NB_2D(SIZE_MU=[1,3], m_observed=m_observed[1:100],BETA=BETA_vec)




function MarginalF_NB_1D(;SIZE,MU, m_observed,BETA)
    temp_vec_2=Array{Float64,1}(undef, length(m_observed))
    #CC = Nemo.ComplexField(64)

    for i = 1:length(m_observed)
        dd=NegativeBinomial(SIZE,SIZE/(SIZE+BETA[i]*MU))
        temp_vec_2[i]=pdf(dd,m_observed[i])
    end
    
    MarginalVal=-sum(log.(temp_vec_2))
    #MarginalVal=-Nemo.sum(Nemo.log.(CC.(temp_vec_2)))
    return MarginalVal
    #return convert.(Float64,real(MarginalVal))
end    


function MarginalF_NB_2D(;SIZE_MU,m_observed,BETA) 

    temp_vec_2=Array{Float64,1}(undef, length(m_observed))

    for i = 1:length(m_observed)
        dd=NegativeBinomial(SIZE_MU[1],SIZE_MU[1]/(SIZE_MU[1]+BETA[i]*SIZE_MU[2]))
        temp_vec_2[i]=pdf(dd,m_observed[i])
    end

    MarginalVal=-sum(log.(temp_vec_2));
    return MarginalVal;
end





function Main_NB_Bay(;Data,BETA_vec,input_size,input_mu,S)

    M=Data
    M_ave=input_mu

    nrow=size(Data)[1]
    ncol=size(Data)[2]
    Final_mat=Array{Int64,3}(undef, nrow, ncol, S)

    for i=1:ncol
        for j=1:nrow
            if M[j,i]==NaN
                for q=1:S
                    Final_mat[j,i,q]=NaN;
                end
            else

                tempmu=(M[j,i]+input_size[j])*(M_ave[j]-M_ave[j]*BETA_vec[i])/(input_size[j]+M_ave[j]*BETA_vec[i])

                dd=NegativeBinomial(input_size[j]+M[j,i],(input_size[j]+M[j,i])/(input_size[j]+M[j,i]+tempmu))
                S_temp=rand(dd, S) .+ M[j,i]
                Final_mat[j,i,:]=S_temp
            end

        end
    end


    return(Final_mat)
end

function Main_mean_NB_Bay(;Data,BETA_vec,input_size,input_mu,S)
    M=Data
    M_ave=input_mu

    nrow=size(Data)[1]
    ncol=size(Data)[2]
    Final_mat=Array{Float64,2}(undef, nrow, ncol)


    for i=1:ncol
        for j=1:nrow
            if M[j,i]==NaN
                Final_mat[j,i]=NaN;
            else
                tempmu=(M[j,i]+input_size[j])*(M_ave[j]-M_ave[j]*BETA_vec[i])/(input_size[j]+M_ave[j]*BETA_vec[i])
                Final_mat[j,i]=M[j,i]+tempmu
            end
        end
    end

    return(Final_mat)
end


function Main_mode_NB_Bay(;Data,BETA_vec,input_size,input_mu,S)
    M=Data
    M_ave=input_mu

    nrow=size(Data)[1]
    ncol=size(Data)[2]
    Final_mat=Array{Float64,2}(undef, nrow, ncol)


    for i=1:ncol
        for j=1:nrow
            if M[j,i]==NaN
                Final_mat[j,i]=NaN;
            elseif (input_size[j]+M[j,i])<=1
                Final_mat[j,i]=0
            elseif (input_size[j]+M[j,i])>1
                meann=(M[j,i]+input_size[j])*(M_ave[j]-M_ave[j]*BETA_vec[i])/(input_size[j]+M_ave[j]*BETA_vec[i])

                Final_mat[j,i]= floor(meann/(input_size[j]+M[j,i])*((input_size[j]+M[j,i])-1))+M[j,i];
            end
        end
    end

    return(Final_mat)
end


# oo=Main_NB_Bay(;Data=RAW_DAT[1:3,1:5],BETA_vec=BETA_vec[1:5],input_size=[1,2,3],input_mu=[4,5,6],S=5)
# oo2=Main_mean_NB_Bay(;Data=RAW_DAT[1:3,1:5],BETA_vec=BETA_vec[1:5],input_size=[1,2,3],input_mu=[4,5,6])
# oo3=Main_mode_NB_Bay(;Data=RAW_DAT[1:3,1:5],BETA_vec=BETA_vec[1:5],input_size=[1,2,3],input_mu=[4,5,6])