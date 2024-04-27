# Replication details for [TO FILL](TO FILL)

This is the replication code for providing lower bounds for the spectral gaps of the cohomological Laplacians $\Delta_1$ for the symplectic groups $\text{Sp}_4(\mathbb{Z})$ and $\text{Sp}_6(\mathbb{Z})$. 

The details concerning the definitions of $\Delta_1$ and the spectral gap can be found in the [arXiv preprint](TO FILL) concerning this repository, as well as in the README file of the [LowCohomologySOS](https://github.com/piotrmizerka/LowCohomologySOS.git) repository which provides some of the necessary procedures for this repository as well.

For the computations we used julia in version `1.9.4` but in principle any later version should work.

## Obtaining code
To obtain the code for the replication, you can either download it directly from [Zenodo](TO FILL), or use git for this. In the latter case, first clone this repository via
```bash
git clone https://github.com/jtszy/SP_4_Cohomology.git
```
and checkout to the correct branch
```bash
cd SP_4_Cohomology
git checkout TO FILL
```

## Setting up the environment
First, run julia in the terminal in `SP_4_Cohomology` folder
```bash
julia --project=.
```
Next, to set up the proper environment for the replication run in julia REPL
```julia
julia> using Pkg; Pkg.instantiate()
```
This command installs and precompiles, if needed, all the necessary dependences,
so it may take a while.
Note that this step needs to be executed only once per installation.

## Running actual replication
We wish to prove that for the Steinberg presentations of $\text{Sp}_4(\mathbb{Z})$ and $\text{Sp}_6(\mathbb{Z})$
on $12$ and $18$ generators respectively (as defined in TO FILL SECTION [TO FILL](TO FILL))
$\Delta_1-\lambda I_{12}$ and $\Delta_1-\lambda I_{18}$ is a sum of squares for $\lambda=0.0833$ and $\lambda=0.0302$ respectively.

In addition, we provide a srcipt for the estimation of the spectral gap of $\Delta_1$ for the Behr presentation on six generators of $\text{Sp}_4(\mathbb{Z})$ (details to be found in TO FILL SECTION of [TO FILL](TO FILL)). The script certifies that $\Delta_1-\lambda I_6$ is a sum of squares for $\lambda=0.0789$.

Our scripts perform the necessary optimizations to find such sums of squares decomposition.

### Steinberg presentations of $\text{Sp}_4(\mathbb{Z})$ and $\text{Sp}_6(\mathbb{Z})$
The following command needs to be executed in the terminal in `SP_4_Cohomology` folder:
```bash
julia --project=. ./scripts/SP_2N_Steinberg.jl n
```
for running the whole replication script for $\text{Sp}_{2n}(\mathbb{Z})$ for $n=2$ and $n=3$.

The running time of the script will be approximately `25` minutes and `115` hours on a standard laptop computer for the cases $n=2$ and $n=3$ respectively. Therefore, in the latter case, we encourage to use the precomputed solution, focusing on the part providing the rigorous proof only.

To run the script which uses the precomputed solution for $\text{Sp}_6(\mathbb{Z})$ (stored in the file "Steinberg_Solution_Sp_6.sjl", where $n=2,3$) and provides rigorous proof (certification) of the result, the following command needs to be executed in the terminal in `SP_4_Cohomology` folder:
```bash
julia --project=. ./scripts/SP_6_replication_precomputed/SP_6_Steinberg_precomputed.jl
```
The running time of the above script will be approximately `70` minutes on a standard laptop computer.

### Behr presentation of $\text{Sp}_4(\mathbb{Z})$
The following command needs to be executed in the terminal in `SP_4_Cohomology` folder:
```bash
julia --project=. ./scripts/SP_4_Behr.jl
```
for running the whole replication script for $\text{Sp}_{4}(\mathbb{Z})$.

The running time of the script will be approximately `8` minutes on a standard laptop computer.
