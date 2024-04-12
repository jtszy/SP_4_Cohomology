@testset "relations_St" begin

    Sp_4 = MatrixGroups.SymplecticGroup{4}(Int8)

    S_4 = gens(Sp_4)

    F_8 = FreeGroup(8)

    gen_dict = Dict(LowCohomologySOS.determine_letter(S_4[i]) => gens(F_8, i) for i in eachindex(S_4))

    x(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:A,i,j)]
    y(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:B,max(i,j),min(i,j) + 2)]
    yt(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:B,min(i,j) + 2,max(i,j))]
    z(i) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:B,i,i + 2)]
    zt(i) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:B,i + 2,i)]

    com(x,y) = SP_4_Cohomology.com(x,y)

    test_relators_Sp_4 = [
        com(x(1, 2), y(1, 2)) * z(1)^(-2),
        com(x(1, 2), yt(1, 2)) * zt(2)^2,
        com(x(1, 2), z(2)) * z(1)^(-1) * y(1, 2)^(-1),
        com(z(1), y(1, 2)),
        com(x(1, 2), zt(1)) * yt(1, 2) * zt(2)^(-1),
        com(zt(2), yt(1, 2)^(-1)),
        com(y(1, 2), zt(1)) * z(2) * x(2, 1)^(-1),
        com(x(2, 1), z(2)^(-1)),
        com(yt(1, 2), z(1)) * zt(2) * x(1, 2),
        com(x(1, 2), zt(2)),
        com(x(2, 1), y(2, 1)) * z(2)^(-2),
        com(x(2, 1), yt(2, 1)) * zt(1)^2,
        com(x(2, 1), z(1))* z(2)^(-1) * y(2, 1)^(-1),
        com(z(2), y(2, 1)),
        com(x(2, 1), zt(2)) * yt(2, 1) * zt(1)^(-1),
        com(zt(1), yt(2, 1)^(-1)),
        com(y(2, 1), zt(2)) * z(1) * x(1, 2)^(-1),
        com(x(1, 2), z(1)^(-1)),
        com(yt(2, 1), z(2)) * zt(1) * x(2, 1),
        com(x(2, 1), zt(1)),
    ]

    computed_relators_Sp_4 = SP_4_Cohomology.relations_St(
        F_8,
        S_4,
        2
    )

    @test Set(test_relators_Sp_4) == Set(computed_relators_Sp_4)

    Sp_6 = MatrixGroups.SymplecticGroup{6}(Int8)

    S_6 = gens(Sp_6)

    F_18 = FreeGroup(18)

    gen_dict = Dict(LowCohomologySOS.determine_letter(S_6[i]) => gens(F_18, i) for i in eachindex(S_6))

    x(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{6}(:A,i,j)]
    y(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{6}(:B,max(i,j),min(i,j) + 3)]
    yt(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{6}(:B,min(i,j) + 3,max(i,j))]

    test_relators_Sp_6 = [
        com(x(1,2), x(2,3)) * x(1,3) ^ (-1),
        com(x(1,2), y(2,3)) * y(1,3) ^ (-1),
        com(x(1,2), yt(1,3)) * yt(2,3),
        com(x(1,3), x(3,2)) * x(1,2) ^ (-1),
        com(x(1,3), y(3,2)) * y(1,2) ^ (-1),
        com(x(1,3), yt(1,2)) * yt(3,2),
        com(x(2,1), x(1,3)) * x(2,3) ^ (-1),
        com(x(2,1), y(1,3)) * y(2,3) ^ (-1),
        com(x(2,1), yt(2,3)) * yt(1,3),
        com(x(2,3), x(3,1)) * x(2,1) ^ (-1),
        com(x(2,3), y(3,1)) * y(2,1) ^ (-1),
        com(x(2,3), yt(2,1)) * yt(3,1),
        com(x(3,1), x(1,2)) * x(3,2) ^ (-1),
        com(x(3,1), y(1,2)) * y(3,2) ^ (-1),
        com(x(3,1), yt(3,2)) * yt(1,2),
        com(x(3,2), x(2,1)) * x(3,1) ^ (-1),
        com(x(3,2), y(2,1)) * y(3,1) ^ (-1),
        com(x(3,2), yt(3,1)) * yt(2,1),
    ]

    computed_relators_Sp_6 = SP_4_Cohomology.relations_St(
        F_18,
        S_6,
        3,
        sq_adj_ = "adj"
    )

    @test Set(test_relators_Sp_6) == Set(computed_relators_Sp_6)

end