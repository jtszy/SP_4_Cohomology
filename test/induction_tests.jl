@testset "sp2n_sp2m_embedding" begin
    N = 2
    M = 3

    emb = SP_4_Cohomology.sp2n_sp2m_embedding(2*N, 2*M)
    Sp2N = emb.source
    Sp2M = emb.target

    for s in gens(Sp2N)
        emb_s = emb(s)

        s_matrix = MatrixGroups.matrix_repr(s)
        type_s, s_i, s_j, s_transpose = SP_4_Cohomology.elementary_element_info(s_matrix)
        @info s_i, s_j, s_transpose

        type_s ,s_i_M, s_j_M = SP_4_Cohomology.ind_embed((type_s, s_i + s_transpose*N, s_j + (1-s_transpose)*N), 2*N, 2*M)

        @info s
        @test emb_s == Groups.MatrixGroups.ElementarySymplectic{2*M}(type_s, s_i_M, s_j_M)

    end

end
