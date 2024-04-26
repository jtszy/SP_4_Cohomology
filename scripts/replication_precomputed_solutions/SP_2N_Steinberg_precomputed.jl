using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS ÷ 2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS ÷ 2)

using Groups
using IntervalArithmetic
using JuMP
using LowCohomologySOS
using SCS
using Serialization
using SP_4_Cohomology
using SparseArrays

N = Int8(ARGS[1])

Sp_2N = MatrixGroups.SymplecticGroup{2*N}(Int8)

F_Sp_2N_Steinberg = FreeGroup(alphabet(Sp_2N))

S = gens(Sp_2N)

quotient_hom_Steinberg = let source = F_Sp_2N_Steinberg, target = Sp_2N
    Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
end

for i in eachindex(S)
    @assert quotient_hom_Steinberg(gens(F_Sp_2N_Steinberg,i)) == S[i]
    @assert quotient_hom_Steinberg(gens(F_Sp_2N_Steinberg,i)^(-1)) == S[i]^(-1)
end

support_jacobian, min_support = SP_4_Cohomology.symplectic_min_supports(quotient_hom_Steinberg, S)

Steinberg_relations = SP_4_Cohomology.relations_St(F_Sp_2N_Steinberg, S, N)

for r in Steinberg_relations
    @assert quotient_hom_Steinberg(r) == one(Sp_2N)
end

Δ₁, I_N = LowCohomologySOS.spectral_gap_elements(quotient_hom_Steinberg, Steinberg_relations, support_jacobian);

RG = LowCohomologySOS.group_ring(Sp_2N, min_support, star_multiplication = true)

Δ₁ = LowCohomologySOS.embed.(identity, Δ₁, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N, Ref(RG))

# Load the precomputed solution
Solution = deserialize("./replication_precomputed_solutions/Steinberg_Solution_Sp_"*string(2*N)*".sjl")
λ, Q = Solution["lambda"], Solution["Q"]

# Certify the numerical estimate
LowCohomologySOS.certify_sos_decomposition(Δ₁, I_N, λ, Q, min_support)
