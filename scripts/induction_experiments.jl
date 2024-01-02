using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using LowCohomologySOS
using SP_4_Cohomology
using StarAlgebras
using Groups
using PermutationGroups

# TODO: adjust the whole code below to sp2ns

N = 3
M = 4

i_emb = SP_4_Cohomology.sp2n_sp2m_embedding(2*N, 2*M)

Sp2N = i_emb.source
Sp2M = i_emb.target

const half_radius = 2

Sp2N_S = gens(Sp2N);
Sp2N_S_inv = let s = Sp2N_S
    [s; inv.(s)]
end;
Sp2M_S = gens(Sp2M);
Sp2M_S_inv = let s = Sp2M_S
    [s; inv.(s)]
end;

Sp2N_S_2 = []
Sp2M_S_2 = []

for s in Sp2N_S 
    _, i, j, _ = SP_4_Cohomology.elementary_element_info(MatrixGroups.matrix_repr(s))
    if i != j
        push!(Sp2N_S_2, s)
    end
end

for s in Sp2M_S 
    _, i, j, _ = SP_4_Cohomology.elementary_element_info(MatrixGroups.matrix_repr(s))
    if i != j
        push!(Sp2M_S_2, s)
    end
end


Sp2N_half_basis, Sp2N_sizes = Groups.wlmetric_ball(Sp2N_S_inv, radius = half_radius);
Sp2M_half_basis, Sp2M_sizes = Groups.wlmetric_ball(Sp2M_S_inv, radius = half_radius);

# Sp2N_Δ₁, Sp2N_Iₙ, Sp2N_Δ₁⁺, Sp2N_Δ₁⁻ = LowCohomologySOS.laplacians(Sp2N, Sp2N_half_basis, Sp2N_S, sq_adj_ = "adj");
# Sp2M_Δ₁, Sp2M_Iₙ, Sp2M_Δ₁⁺, Sp2M_Δ₁⁻ = LowCohomologySOS.laplacians(Sp2M, Sp2M_half_basis, Sp2M_S, sq_adj_ = "adj");

Sp2N_Δ₁⁻ = SP_4_Cohomology.lower_laplacian(Sp2N, Sp2N_half_basis, Sp2N_S_2);
Sp2M_Δ₁⁻ = SP_4_Cohomology.lower_laplacian(Sp2M, Sp2M_half_basis, Sp2M_S_2);

# Sp2N_mono, Sp2N_sq, Sp2N_adj, Sp2N_op = SP_4_Cohomology.mono_sq_adj_op(Sp2N_Δ₁⁻, Sp2N_S)
# Sp2M_mono, Sp2M_sq, Sp2M_adj, Sp2M_op = SP_4_Cohomology.mono_sq_adj_op(Sp2M_Δ₁⁻, Sp2M_S)
Sp2N_mono, Sp2N_sq, Sp2N_adj, Sp2N_op = SP_4_Cohomology.mono_sq_adj_op(Sp2N_Δ₁⁻, Sp2N_S_2)
Sp2M_mono, Sp2M_sq, Sp2M_adj, Sp2M_op = SP_4_Cohomology.mono_sq_adj_op(Sp2M_Δ₁⁻, Sp2M_S_2)

# RG_prime = parent(first(Sp2M_Δ₁⁺))

# Δ₁⁺_emb = LowCohomologySOS.embed_matrix(Sp2N_Δ₁⁺, i_emb, RG_prime)
# adj_emb = LowCohomologySOS.embed_matrix(Sp2N_adj, i_emb, RG_prime)
RG = parent(first(Sp2M_adj))
adj_emb = LowCohomologySOS.embed_matrix(Sp2N_adj, i_emb, RG, Sp2N_S_2, Sp2M_S_2)

# @assert parent(first(Δ₁⁺_emb)) == parent(first(adj_emb)) == parent(first(Sp2M_Δ₁⁺)) == parent(first(Sp2M_adj))

# Δ₁⁺_emb_symmetrized = let
#     Σ = PermutationGroups.SymmetricGroup(3)
#     LowCohomologySOS.weyl_symmetrize_matrix(Δ₁⁺_emb, Σ, SP_4_Cohomology.conjugation, Sp2M_S)
# end

adj_emb_symmetrized = let
    Σ = PermutationGroups.SymmetricGroup(4)
    LowCohomologySOS.weyl_symmetrize_matrix(adj_emb, Σ, SP_4_Cohomology.conjugation, Sp2M_S_2)
end

# TODO change scalars below:
3*Sp2M_Δ₁⁺-Δ₁⁺_emb_symmetrized # it looks like symmetrization works for upper Laplacians!
zero_matrix = [zero(RG) for i in Sp2M_S_2, j in Sp2M_S_2]
6*Sp2M_adj-adj_emb_symmetrized == zero_matrix
@info (6*Sp2M_adj-adj_emb_symmetrized) # looks like adj symmetrizes if we remove Z from gens

@info Sp2M_Δ₁⁺[1,1]
@info Δ₁⁺_emb_symmetrized[1,1]

@info Sp2M_adj[1,2]
@info adj_emb_symmetrized[1,2]

A = gens(Sp2N)[10]
Σ = PermutationGroups.SymmetricGroup(3)
sigma = gens(Σ)[1]

SP_4_Cohomology.conjugation(A^(-1), sigma)
SP_4_Cohomology.elementary_conjugation(A, sigma)

for s in gens(Sp2N)
    @info s, SP_4_Cohomology.conjugation(s, sigma)
    @info s, SP_4_Cohomology.conjugation(s^(-1), sigma)
end

@info gens(Sp2N)[10], SP_4_Cohomology.conjugation(gens(Sp2N)[10], sigma)
@info gens(Sp2N)[13], SP_4_Cohomology.conjugation(gens(Sp2N)[13], sigma)

B = gens(Sp2M)[19]
B_M = MatrixGroups.matrix_repr(B)

SP_4_Cohomology.elementary_element_info(B_M)
print(B_M)

MatrixGroups.matrix_repr(Groups.MatrixGroups.ElementarySymplectic{6}(:B, 2, 4))
MatrixGroups.matrix_repr(Groups.MatrixGroups.ElementarySymplectic{6}(:B, 1, 5))

B = Groups.MatrixGroups.ElementarySymplectic{6}(:B, 2, 4)
typeof(B)

gens(Sp2N)[7]