# bayNormJL.jl: bayNorm in Julia
bayNorm in Julia, super fast! The R package can be found here (more detailed about parameter setting): `https://github.com/WT215/bayNorm`.

bayNormJL is a Julia package which is used to normalize single-cell RNA-seq data. 

## Code for bayNorm paper
The code for producing figures in bayNorm paper [1] can be found [here](https://github.com/WT215/bayNorm_papercode)

To install the package under Julia 1.2+ using `Pkg`, simply enter its interactive mode by pressing `]` in the Julia REPL and run
`add https://github.com/WT215/bayNormJL.jl.git`

The output of bayNorm is of type `Dict`.

## Quick start
```julia
using bayNormJL

#A toy scRNA-seq data
Data=[1:10 1:10 3:12 4:13] 
#Cell-specific capture efficiency. Default is nothing.
Beta=[0.1,0.5,0.3,0.4]

outt=bayNormJL.bayNorm(Data=Data, BETA_vec=Beta, Conditions = nothing,UMI_sffl = nothing,Prior_type = nothing,mode_version =false,mean_version=true,S = 20,FIX_MU = true,BB_SIZE_par = true, verbose = true)

#Access to the normalized data
outt["Bay_out"]
```
## Parallel computing
```julia
using Distributed
#Using 8 cores
addprocs(8)
@everwhere using bayNormJL

#A toy scRNA-seq data
@everywhere Data=[1:10 1:10 3:12 4:13] 
#Cell-specific capture efficiency. Default is nothing.
@everywhere Beta=[0.1,0.5,0.3,0.4]

outt=bayNormJL.bayNorm(Data=Data, BETA_vec=Beta, Conditions = nothing,UMI_sffl = nothing,Prior_type = nothing,mode_version =false,mean_version=true,S = 20,FIX_MU = true,BB_SIZE_par = true, verbose = true)

#Access to the normalized data
outt["Bay_out"]
```


## References

- [1] <a href="https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz726/5581401">Tang <em>et al.</em> (2019). Bioinformatics. </a>
