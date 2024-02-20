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
M = 5

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

Sp2N_Δ₁, Sp2N_Iₙ, Sp2N_Δ₁⁺, Sp2N_Δ₁⁻ = LowCohomologySOS.laplacians(Sp2N, Sp2N_half_basis, Sp2N_S, sq_adj_ = "adj");
Sp2M_Δ₁, Sp2M_Iₙ, Sp2M_Δ₁⁺, Sp2M_Δ₁⁻ = LowCohomologySOS.laplacians(Sp2M, Sp2M_half_basis, Sp2M_S, sq_adj_ = "adj");

# Sp2N_Δ₁, Sp2N_Iₙ, Sp2N_Δ₁⁺, Sp2N_Δ₁⁻ = SP_4_Cohomology.non_mono_laplacians(Sp2N, Sp2N_half_basis, Sp2N_S_2, twist_coeffs = false);
# Sp2M_Δ₁, Sp2M_Iₙ, Sp2M_Δ₁⁺, Sp2M_Δ₁⁻ = SP_4_Cohomology.non_mono_laplacians(Sp2M, Sp2M_half_basis, Sp2M_S_2, twist_coeffs = false);

Sp2N_mono, Sp2N_sq_mono, Sp2N_sq_mix, Sp2N_sq_double, Sp2N_adj_mix, Sp2N_adj_double, Sp2N_op = SP_4_Cohomology.mono_sq_adj_op(Sp2N_Δ₁⁻, Sp2N_S)
Sp2M_mono, Sp2M_sq_mono, Sp2M_sq_mix, Sp2M_sq_double, Sp2M_adj_mix, Sp2M_adj_double, Sp2M_op = SP_4_Cohomology.mono_sq_adj_op(Sp2M_Δ₁⁻, Sp2M_S)

# Sp2N_mono, Sp2N_sq, Sp2N_adj, Sp2N_op = SP_4_Cohomology.mono_sq_adj_op(Sp2N_Δ₁⁻, Sp2N_S_2)
# Sp2M_mono, Sp2M_sq, Sp2M_adj, Sp2M_op = SP_4_Cohomology.mono_sq_adj_op(Sp2M_Δ₁⁻, Sp2M_S_2)

RG = parent(first(Sp2M_Δ₁⁺))

# For non-mono part 

# Δ₁⁺_emb = SP_4_Cohomology.non_mono_embed_matrix(Sp2N_Δ₁⁺, i_emb, RG, Sp2N_S_2, Sp2M_S_2)
# adj_emb = SP_4_Cohomology.non_mono_embed_matrix(Sp2N_adj, i_emb, RG, Sp2N_S_2, Sp2M_S_2)

# @assert parent(first(Δ₁⁺_emb)) == parent(first(adj_emb)) == parent(first(Sp2M_Δ₁⁺)) == parent(first(Sp2M_adj))

# Δ₁⁺_emb_symmetrized = let
#     Σ = PermutationGroups.SymmetricGroup(4)
#     LowCohomologySOS.weyl_symmetrize_matrix(Δ₁⁺_emb, Σ, SP_4_Cohomology.conjugation, Sp2M_S_2)
# end

# adj_emb_symmetrized = let
#     Σ = PermutationGroups.SymmetricGroup(4)
#     LowCohomologySOS.weyl_symmetrize_matrix(adj_emb, Σ, SP_4_Cohomology.conjugation, Sp2M_S_2)
# end

Δ₁⁺_emb = LowCohomologySOS.embed_matrix(Sp2N_Δ₁⁺, i_emb, RG)
adj_double_emb = LowCohomologySOS.embed_matrix(Sp2N_adj_double, i_emb, RG)
adj_mono_emb = LowCohomologySOS.embed_matrix(Sp2N_adj_mono, i_emb, RG)
sq_mix_emb = LowCohomologySOS.embed_matrix(Sp2N_sq_mix, i_emb, RG)

@assert parent(first(Δ₁⁺_emb)) == parent(first(adj_double_emb)) == parent(first(Sp2M_Δ₁⁺))

sq_mix_emb_symmetrized = let
    Σ = PermutationGroups.SymmetricGroup(M)
    LowCohomologySOS.weyl_symmetrize_matrix(sq_mix_emb, Σ, SP_4_Cohomology.conjugation, Sp2M_S)
end

Δ₁⁺_emb_symmetrized = let
    Σ = PermutationGroups.SymmetricGroup(M)
    LowCohomologySOS.weyl_symmetrize_matrix(Δ₁⁺_emb, Σ, SP_4_Cohomology.conjugation, Sp2M_S)
end

adj_mono_emb_symmetrized = let
    Σ = PermutationGroups.SymmetricGroup(M)
    LowCohomologySOS.weyl_symmetrize_matrix(adj_mono_emb, Σ, SP_4_Cohomology.conjugation, Sp2M_S)
end

adj_double_emb_symmetrized = let
    Σ = PermutationGroups.SymmetricGroup(M)
    LowCohomologySOS.weyl_symmetrize_matrix(adj_double_emb, Σ, SP_4_Cohomology.conjugation, Sp2M_S)
end


# TODO change scalars below:
zero_matrix = [zero(RG) for i in Sp2M_S, j in Sp2M_S]

12*Sp2M_Δ₁⁺-Δ₁⁺_emb_symmetrized == zero_matrix# it looks like symmetrization works for upper Laplacians!

12*Sp2M_adj_double-adj_double_emb_symmetrized == zero_matrix

12*Sp2M_adj_mono - adj_mono_emb_symmetrized == zero_matrix

36*Sp2M_sq_mix - sq_mix_emb_symmetrized == zero_matrix