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
end