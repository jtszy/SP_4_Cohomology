@testset "quotient_homomorphism" begin
    F_2 = Groups.FreeGroup(2)
    a, b = Groups.gens(F_2)
    Z2 = Groups.FPGroup(F_2, [a * b => b * a])
    A, B = Groups.gens(Z2)

    test_group = Dict()
    test_group[a] = A
    test_group[b] = B

    hom_Z2 = SP_4_Cohomology.quotient_homomorphism(F_2, Z2, test_group)

    @test hom_Z2(a) == A
    @test hom_Z2(b) == B
    @test hom_Z2(a * b) == A * B
    @test hom_Z2(b * a) == A * B
end

@testset "com" begin
    F_2 = Groups.FreeGroup(2)
    a, b = Groups.gens(F_2)
    Z2 = Groups.FPGroup(F_2, [a * b => b * a])
    A, B = Groups.gens(Z2)

    @test SP_4_Cohomology.com(A, B) == one(Z2)

    sl_3 = MatrixGroups.SpecialLinearGroup{3}(Int8)
    e12, e13, e21, e23, e31, e32 = Groups.gens(sl_3)

    @test SP_4_Cohomology.com(e12, e13) == one(sl_3)
    @test SP_4_Cohomology.com(e12, e32) == oen(sl_3)
    @test SP_4_Cohomology.com(e13, e23) == one(sl_3)
    @test SP_4_Cohomology.com(e21, e23) == one(sl_3)
    @test SP_4_Cohomology.com(e21, e31) == one(sl_3)
    @test SP_4_Cohomology.com(e31, e32) == one(sl_3)
    @test SP_4_Cohomology.com(e12, e23) == e13
    @test SP_4_Cohomology.com(e13, e32) == e12
    @test SP_4_Cohomology.com(e21, e13) == e23
    @test SP_4_Cohomology.com(e23, e31) == e21
    @test SP_4_Cohomology.com(e31, e12) == e32
    @test SP_4_Cohomology.com(e32, e21) == e31
end

@testset "symplectic_min_supports" begin
    Sp_4 = MatrixGroups.SymplecticGroup{4}(Int8)

    F_Sp_4_Steinberg = FreeGroup(alphabet(Sp_4))

    S = gens(Sp_4)

    quotient_hom_Steinberg = let source = F_Sp_4_Steinberg, target = Sp_4
        Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
    end

    support_jacobian, min_support =
        SP_4_Cohomology.symplectic_min_supports(quotient_hom_Steinberg, S)

    I = MatrixGroups.matrix_repr(one(Sp_4))

    x(i, j) = MatrixGroups.matrix_repr(MatrixGroups.ElementarySymplectic{4}(:A, i, j))
    y(i, j) = MatrixGroups.matrix_repr(
        MatrixGroups.ElementarySymplectic{4}(:B, max(i, j), min(i, j) + 2),
    )
    yt(i, j) = MatrixGroups.matrix_repr(
        MatrixGroups.ElementarySymplectic{4}(:B, min(i, j) + 2, max(i, j)),
    )
    z(i) = MatrixGroups.matrix_repr(MatrixGroups.ElementarySymplectic{4}(:B, i, i + 2))
    zt(i) = MatrixGroups.matrix_repr(MatrixGroups.ElementarySymplectic{4}(:B, i + 2, i))

    test_support_jacobian = vcat(
        [
            x(1, 2),
            x(1, 2) * y(1, 2),
            x(1, 2) * y(1, 2) * x(1, 2)^(-1),
            x(1, 2) * y(1, 2) * x(1, 2)^(-1) * y(1, 2)^(-1),
            x(1, 2) * y(1, 2) * x(1, 2)^(-1) * y(1, 2)^(-1) * z(1)^(-1),
            x(1, 2) * y(1, 2) * x(1, 2)^(-1) * y(1, 2)^(-1) * z(1)^(-2),
        ],
        [
            x(1, 2),
            x(1, 2) * yt(1, 2),
            x(1, 2) * yt(1, 2) * x(1, 2)^(-1),
            x(1, 2) * yt(1, 2) * x(1, 2)^(-1) * yt(1, 2)^(-1),
            x(1, 2) * yt(1, 2) * x(1, 2)^(-1) * yt(1, 2)^(-1) * zt(2),
            x(1, 2) * yt(1, 2) * x(1, 2)^(-1) * yt(1, 2)^(-1) * zt(2)^2,
        ],
        [
            x(1, 2),
            x(1, 2) * z(2),
            x(1, 2) * z(2) * x(1, 2)^(-1),
            x(1, 2) * z(2) * x(1, 2)^(-1) * z(2)^(-1),
            x(1, 2) * z(2) * x(1, 2)^(-1) * z(2)^(-1) * z(1)^(-1),
            x(1, 2) * z(2) * x(1, 2)^(-1) * z(2)^(-1) * z(1)^(-1) * y(1, 2)^(-1),
        ],
        [
            z(1),
            z(1) * y(1, 2),
            z(1) * y(1, 2) * z(1)^(-1),
            z(1) * y(1, 2) * z(1)^(-1) * y(1, 2)^(-1),
        ],
        [
            x(1, 2),
            x(1, 2) * zt(1),
            x(1, 2) * zt(1) * x(1, 2)^(-1),
            x(1, 2) * zt(1) * x(1, 2)^(-1) * zt(1)^(-1),
            x(1, 2) * zt(1) * x(1, 2)^(-1) * zt(1)^(-1) * yt(1, 2),
            x(1, 2) * zt(1) * x(1, 2)^(-1) * zt(1)^(-1) * yt(1, 2) * zt(2)^(-1),
        ],
        [
            zt(2),
            zt(2) * yt(1, 2)^(-1),
            zt(2) * yt(1, 2)^(-1) * zt(2)^(-1),
            zt(2) * yt(1, 2)^(-1) * zt(2)^(-1) * yt(1, 2),
        ],
        [
            y(1, 2),
            y(1, 2) * zt(1),
            y(1, 2) * zt(1) * y(1, 2)^(-1),
            y(1, 2) * zt(1) * y(1, 2)^(-1) * zt(1)^(-1),
            y(1, 2) * zt(1) * y(1, 2)^(-1) * zt(1)^(-1) * z(2),
            y(1, 2) * zt(1) * y(1, 2)^(-1) * zt(1)^(-1) * z(2) * x(2, 1)^(-1),
        ],
        [
            x(2, 1),
            x(2, 1) * z(2)^(-1),
            x(2, 1) * z(2)^(-1) * x(2, 1)^(-1),
            x(2, 1) * z(2)^(-1) * x(2, 1)^(-1) * z(2),
        ],
        [
            yt(1, 2),
            yt(1, 2) * z(1),
            yt(1, 2) * z(1) * yt(1, 2)^(-1),
            yt(1, 2) * z(1) * yt(1, 2)^(-1) * z(1)^(-1),
            yt(1, 2) * z(1) * yt(1, 2)^(-1) * z(1)^(-1) * zt(2),
            yt(1, 2) * z(1) * yt(1, 2)^(-1) * z(1)^(-1) * zt(2) * x(1, 2),
        ],
        [
            x(1, 2),
            x(1, 2) * zt(2),
            x(1, 2) * zt(2) * x(1, 2)^(-1),
            x(1, 2) * zt(2) * x(1, 2)^(-1) * zt(2)^(-1),
        ],
        [
            x(2, 1),
            x(2, 1) * y(2, 1),
            x(2, 1) * y(2, 1) * x(2, 1)^(-1),
            x(2, 1) * y(2, 1) * x(2, 1)^(-1) * y(2, 1)^(-1),
            x(2, 1) * y(2, 1) * x(2, 1)^(-1) * y(2, 1)^(-1) * z(2)^(-1),
            x(2, 1) * y(2, 1) * x(2, 1)^(-1) * y(2, 1)^(-1) * z(2)^(-2),
        ],
        [
            x(2, 1),
            x(2, 1) * yt(2, 1),
            x(2, 1) * yt(2, 1) * x(2, 1)^(-1),
            x(2, 1) * yt(2, 1) * x(2, 1)^(-1) * yt(2, 1)^(-1),
            x(2, 1) * yt(2, 1) * x(2, 1)^(-1) * yt(2, 1)^(-1) * zt(1),
            x(2, 1) * yt(2, 1) * x(2, 1)^(-1) * yt(2, 1)^(-1) * zt(1)^2,
        ],
        [
            x(2, 1),
            x(2, 1) * z(1),
            x(2, 1) * z(1) * x(2, 1)^(-1),
            x(2, 1) * z(1) * x(2, 1)^(-1) * z(1)^(-1),
            x(2, 1) * z(1) * x(2, 1)^(-1) * z(1)^(-1) * z(2)^(-1),
            x(2, 1) * z(1) * x(2, 1)^(-1) * z(1)^(-1) * z(2)^(-1) * y(2, 1)^(-1),
        ],
        [
            z(2),
            z(2) * y(2, 1),
            z(2) * y(2, 1) * z(2)^(-1),
            z(2) * y(2, 1) * z(2)^(-1) * y(2, 1)^(-1),
        ],
        [
            x(2, 1),
            x(2, 1) * zt(2),
            x(2, 1) * zt(2) * x(2, 1)^(-1),
            x(2, 1) * zt(2) * x(2, 1)^(-1) * zt(2)^(-1),
            x(2, 1) * zt(2) * x(2, 1)^(-1) * zt(2)^(-1) * yt(2, 1),
            x(2, 1) * zt(2) * x(2, 1)^(-1) * zt(2)^(-1) * yt(2, 1) * zt(1)^(-1),
        ],
        [
            zt(1),
            zt(1) * yt(2, 1)^(-1),
            zt(1) * yt(2, 1)^(-1) * zt(1)^(-1),
            zt(1) * yt(2, 1)^(-1) * zt(1)^(-1) * yt(2, 1),
        ],
        [
            y(2, 1),
            y(2, 1) * zt(2),
            y(2, 1) * zt(2) * y(2, 1)^(-1),
            y(2, 1) * zt(2) * y(2, 1)^(-1) * zt(2)^(-1),
            y(2, 1) * zt(2) * y(2, 1)^(-1) * zt(2)^(-1) * z(1),
            y(2, 1) * zt(2) * y(2, 1)^(-1) * zt(2)^(-1) * z(1) * x(1, 2)^(-1),
        ],
        [
            x(1, 2),
            x(1, 2) * z(1)^(-1),
            x(1, 2) * z(1)^(-1) * x(1, 2)^(-1),
            x(1, 2) * z(1)^(-1) * x(1, 2)^(-1) * z(1),
        ],
        [
            yt(2, 1),
            yt(2, 1) * z(2),
            yt(2, 1) * z(2) * yt(2, 1)^(-1),
            yt(2, 1) * z(2) * yt(2, 1)^(-1) * z(2)^(-1),
            yt(2, 1) * z(2) * yt(2, 1)^(-1) * z(2)^(-1) * zt(1),
            yt(2, 1) * z(2) * yt(2, 1)^(-1) * z(2)^(-1) * zt(1) * x(2, 1),
        ],
        [
            x(2, 1),
            x(2, 1) * zt(1),
            x(2, 1) * zt(1) * x(2, 1)^(-1),
            x(2, 1) * zt(1) * x(2, 1)^(-1) * zt(1)^(-1),
        ],
    )
    unique!(test_support_jacobian)

    @test length(support_jacobian) == length(test_support_jacobian)
    for el in support_jacobian
        @test MatrixGroups.matrix_repr(el) in test_support_jacobian
    end

    sup_sq_relations_derived(i, j) = vcat(
        [z(i)^2 * y(i, j), z(i)^2],
        [(zt(j)^(-2)) * yt(i, j), (zt(j)^(-2))],
        [y(i, j) * z(i) * z(j), y(i, j) * z(i)],
        [zt(j) * (yt(i, j))^(-1) * zt(i), zt(j) * (yt(i, j))^(-1)],
        [zt(j) * (yt(i, j))^(-1)],
        [x(j, i) * z(j)^(-1) * zt(i), x(j, i) * z(j)^(-1)],
        [x(j, i) * z(j)^(-1)],
        [x(i, j)^(-1) * zt(j)^(-1) * z(i), x(i, j)^(-1) * zt(j)^(-1)],
    )

    test_min_support = vcat(
        sup_sq_relations_derived(1, 2),
        sup_sq_relations_derived(2, 1),
        inv.(sup_sq_relations_derived(1, 2)),
        inv.(sup_sq_relations_derived(2, 1)),
        [I, x(1, 2), x(2, 1), y(1, 2), yt(1, 2), z(1), z(2), zt(1), zt(2)],
        inv.([x(1, 2), x(2, 1), y(1, 2), yt(1, 2), z(1), z(2), zt(1), zt(2)]),
    )

    unique!(test_min_support)

    @test length(min_support) == length(test_min_support)
    for el in min_support
        @test MatrixGroups.matrix_repr(el) in test_min_support
    end

end
