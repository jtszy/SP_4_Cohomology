@testset "sp2n_sp2m_embedding" begin
    N = 3
    M = 4

    emb = SP_4_Cohomology.sp2n_sp2m_embedding(2*N, 2*M)
    Sp2N = emb.source
    Sp2M = emb.target

    for s in gens(Sp2N)
        emb_s = MatrixGroups.matrix_repr(emb(s))

        s_matrix = MatrixGroups.matrix_repr(s)
        type_s, s_i, s_j, s_transpose = SP_4_Cohomology.elementary_element_info(s_matrix)
        @info s_i, s_j, s_transpose
        @info s

        if type_s == :A 
            @test emb_s == MatrixGroups.matrix_repr(Groups.MatrixGroups.ElementarySymplectic{2*M}(type_s, s_i, s_j))

        elseif s_transpose == 0
            @test emb_s == MatrixGroups.matrix_repr(Groups.MatrixGroups.ElementarySymplectic{2*M}(type_s, max(s_i, s_j), min(s_i, s_j) + M))
        
        else
            @test emb_s == MatrixGroups.matrix_repr(Groups.MatrixGroups.ElementarySymplectic{2*M}(type_s, min(s_i, s_j) + M, max(s_i, s_j)))
        end
    end
end

# In Sp_4:

# Y_1_2 -> B_2_3

# Y_1_2^T -> B_3_2

@testset "lower_laplacians" begin
    N = 6
    Sp_N = MatrixGroups.SymplecticGroup{N}(Int8)
    generators_Sp_N = gens(Sp_N)
    sq_adj_ = "adj"
    half_basis, sizes = Groups.wlmetric_ball([generators_Sp_N; inv.(generators_Sp_N)], radius = 2)

    Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.laplacians(
        Sp_N,
        half_basis,
        generators_Sp_N,
        sq_adj_ = sq_adj_,
        twist_coeffs = false
    )

    @test size(Δ₁⁻) == size(Δ₁⁺) == size(Δ₁) == 
        size(Iₙ) == (length(generators_Sp_N), length(generators_Sp_N))

    @test parent(first(Δ₁⁻)) == parent(first(Δ₁⁺)) == parent(first(Δ₁)) == parent(first(Iₙ))
    
    RG = parent(first(Δ₁))

    for i in eachindex(generators_Sp_N)
        for j in eachindex(generators_Sp_N)
            s = generators_Sp_N[i]
            t = generators_Sp_N[j]
            @test (one(RG) - RG(s))*(one(RG) - RG(t))' == Δ₁⁻[i,j]
        end

    end

end

@testset "relations" begin

    N = 6
    Sp_N = MatrixGroups.SymplecticGroup{N}(Int8)
    S = gens(Sp_N)

    F = FreeGroup(alphabet(Sp_N))

    gen_dict = Dict(LowCohomologySOS.determine_letter(S[i]) => gens(F, i) for i in eachindex(S))

    x(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{N}(:A,i,j)]
    y(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{N}(:B,max(i,j),min(i,j) + div(N,2))]
    yt(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{N}(:B,min(i,j) + div(N,2),max(i,j))]
    z(i) = gen_dict[MatrixGroups.ElementarySymplectic{N}(:B,i,i + div(N,2))]
    zt(i) = gen_dict[MatrixGroups.ElementarySymplectic{N}(:B,i + div(N,2),i)]

    com(s,t) = SP_4_Cohomology.com(s,t)

    relations = [
        com(x(1,2), x(2,3)) * x(1,3)^(-1),
        com(x(1,3), x(3,2)) * x(1,2)^(-1),
        com(x(2,1), x(1,3)) * x(2,3)^(-1),
        com(x(2,3), x(3,1)) * x(2,1)^(-1),
        com(x(3,1), x(1,2)) * x(3,2)^(-1),
        com(x(3,2), x(2,1)) * x(3,1)^(-1),
        com(x(1,2), y(2,3)) * y(1,3)^(-1),
        com(x(1,3), y(2,3)) * y(1,2)^(-1),
        com(x(2,1), y(1,3)) * y(2,3)^(-1),
        com(x(2,3), y(1,3)) * y(1,2)^(-1),
        com(x(3,1), y(1,2)) * y(2,3)^(-1),
        com(x(3,2), y(1,2)) * y(1,3)^(-1),
        com(x(1,2), yt(1,3)) * yt(2,3),
        com(x(1,3), yt(1,2)) * yt(2,3),
        com(x(2,1), yt(2,3)) * yt(1,3),
        com(x(2,3), yt(1,2)) * yt(1,3),
        com(x(3,1), yt(2,3)) * yt(1,2),
        com(x(3,2), yt(1,3)) * yt(1,2)
    ]

    relations_LowCohomology = LowCohomologySOS.relations(F,S,div(N,2), "adj")

    @test size(relations) == size(relations_LowCohomology)

    for r in relations
        @test r ∈ relations_LowCohomology
    end

end

