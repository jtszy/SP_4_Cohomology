using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS ÷ 2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS ÷ 2)

using Groups
using IntervalArithmetic
using JuMP
using LowCohomologySOS
using PermutationGroups
using SCS
using Serialization
using SP_4_Cohomology
using SparseArrays
using SymbolicWedderburn

N = 2

Sp_N = MatrixGroups.SymplecticGroup{2*N}(Int8)

F_Sp_N_Steinberg = FreeGroup(alphabet(Sp_N))

S = gens(Sp_N)

quotient_hom_Steinberg = let source = F_Sp_N_Steinberg, target = Sp_N
    Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
end

for i in eachindex(S)
    @assert quotient_hom_Steinberg(gens(F_Sp_N_Steinberg,i)) == S[i]
    @assert quotient_hom_Steinberg(gens(F_Sp_N_Steinberg,i)^(-1)) == S[i]^(-1)
end

support_jacobian, min_support = SP_4_Cohomology.symplectic_min_supports(quotient_hom_Steinberg, S)

Steinberg_relations = SP_4_Cohomology.relations_St(F_Sp_N_Steinberg, S, N)

Δ₁, I_N = LowCohomologySOS.spectral_gap_elements(quotient_hom_Steinberg, Steinberg_relations, support_jacobian);

RG = LowCohomologySOS.group_ring(Sp_N, min_support, star_multiplication = true)

Δ₁ = LowCohomologySOS.embed.(identity, Δ₁, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N, Ref(RG))

constraints_basis, psd_basis, Σ, action = SP_4_Cohomology.wedderburn_data(RG.basis, min_support, S);

# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    @assert LowCohomologySOS.act_on_matrix(Δ₁, σ, action.alphabet_perm, S) == Δ₁
    @assert LowCohomologySOS.act_on_matrix(I_N, σ, action.alphabet_perm, S) == I_N
end

@time begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

@time begin
    sos_problem, P = LowCohomologySOS.sos_problem(
        Δ₁,
        I_N,
        w_dec_matrix
        # 0.7 / 0.05
    )
end

# Find a numerical spectral gap
JuMP.set_optimizer(sos_problem, SP_4_Cohomology.scs_opt(eps = 1e-5, max_iters = 500))
JuMP.optimize!(sos_problem)

# Certify the numerical estimate
λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)
LowCohomologySOS.certify_sos_decomposition(M, I, λ, Q, support_jacobian)