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

    @test SP_4_Cohomology.com(A,B) == one(Z2)

    sl_3 = MatrixGroups.SpecialLinearGroup{3}(Int8)
    e12, e13, e21, e23, e31, e32 = Groups.gens(sl_3)

    @test SP_4_Cohomology.com(e12,e13) == one(sl_3)
    @test SP_4_Cohomology.com(e12,e32) == oen(sl_3)
    @test SP_4_Cohomology.com(e13,e23) == one(sl_3)
    @test SP_4_Cohomology.com(e21,e23) == one(sl_3)
    @test SP_4_Cohomology.com(e21,e31) == one(sl_3)
    @test SP_4_Cohomology.com(e31,e32) == one(sl_3)
    @test SP_4_Cohomology.com(e12,e23) == e13
    @test SP_4_Cohomology.com(e13,e32) == e12
    @test SP_4_Cohomology.com(e21,e13) == e23
    @test SP_4_Cohomology.com(e23,e31) == e21
    @test SP_4_Cohomology.com(e31,e12) == e32
    @test SP_4_Cohomology.com(e32,e21) == e31
end

@testset "symplectic_min_supports" begin
    # todo
end