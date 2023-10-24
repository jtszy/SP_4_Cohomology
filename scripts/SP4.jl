using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
using SparseArrays
using JuMP
using SCS
using IntervalArithmetic
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Groups
using LowCohomologySOS

n = 4
sp4 = MatrixGroups.SymplecticGroup{n}(Int8)
a12, a21, b31, b42, b32, b13, b24, b23 = gens(sp4)

# alphabet(Sp4)
# F_8 = FreeGroup(8)
F_sp4 = FreeGroup(alphabet(sp4))
A12, A21, B31, B42, B32, B13, B24, B23 = gens(F_sp4)

quotient_hom = let source = F_sp4, target = sp4
    Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
end

quotient_hom(A12) == a12
quotient_hom(A21) == a21
quotient_hom(B31) == b31
quotient_hom(B42) == b42
quotient_hom(B32) == b32
quotient_hom(B13) == b13
quotient_hom(B24) == b24
quotient_hom(B23) == b23

F_sp4_Behr = FreeGroup(6)
x_a, x_b, x_ab, x_2ab, w_a, w_b = gens(F_sp4_Behr)

x_a_matrix = [
    1 1 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 -1 1;
]

x_b_matrix = [
    1 0 0 0;
    0 1 0 1;
    0 0 1 0;
    0 0 0 1; 
]

x_ab_matrix = [
    1 0 0 1;
    0 1 1 0;
    0 0 1 0;
    0 0 0 1;
]

x_2ab_matrix = [
    1 0 1 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1;
]

w_a_matrix = [
    0 -1 0 0;
    1 0 0 0;
    0 0 0 -1;
    0 0 1 0;
]

w_b_matrix = [
    1 0 0 0;
    0 0 0 -1;
    0 0 1 0;
    0 1 0 0;
]

gens_sp4 = [gens(sp4); inv.(gens(sp4))]

Ball_4, sizes = Groups.wlmetric_ball(gens_sp4, radius = 4)
# Ball_4[18000]

x_a_group = []
x_b_group = []
x_ab_group = []
x_2ab_group = []
w_a_group = []
w_b_group = []

for g in Ball_4
    if MatrixGroups.matrix_repr(g) == x_a_matrix
        x_a_group = g
    elseif MatrixGroups.matrix_repr(g) == x_b_matrix
        x_b_group = g
    elseif MatrixGroups.matrix_repr(g) == x_ab_matrix
        x_ab_group = g
    elseif MatrixGroups.matrix_repr(g) == x_2ab_matrix
        x_2ab_group = g
    elseif MatrixGroups.matrix_repr(g) == w_a_matrix
        w_a_group = g
    elseif MatrixGroups.matrix_repr(g) == w_b_matrix
        w_b_group = g
    end         
end 

gens(F_sp4_Behr)
F_sp4_Behr([1])

quotient_hom_Behr = let source = F_sp4_Behr, target = sp4
    function F(i,source, target)
        if source([i]) == x_a
            return word(x_a_group)
        elseif source([i]) == x_a^(-1)
            return word(x_a_group^(-1))
        elseif  source([i]) == x_b
            return word(x_b_group)
        elseif source([i]) == x_b^(-1)
            return word(x_b_group^(-1))
        elseif  source([i]) == x_ab
            return word(x_ab_group)
        elseif source([i]) == x_ab^(-1)
            return word(x_ab_group^(-1))
        elseif  source([i]) == x_2ab
            return word(x_2ab_group)
        elseif source([i]) == x_2ab^(-1)
            return word(x_2ab_group^(-1))
        elseif  source([i]) == w_a
            return word(w_a_group)
        elseif source([i]) == w_a^(-1)
            return word(w_a_group^(-1))
        elseif  source([i]) == w_b
            return word(w_b_group)
        elseif source([i]) == w_b^(-1)
            return word(w_b_group^(-1))
        end     
    end
    Groups.Homomorphism(F, source, target)
end

quotient_hom_Behr(x_a) == x_a_group
quotient_hom_Behr(x_b) == x_b_group
quotient_hom_Behr(x_ab) == x_ab_group
quotient_hom_Behr(x_2ab) == x_2ab_group
quotient_hom_Behr(w_a) == w_a_group
quotient_hom_Behr(w_b) == w_b_group

Behr_relations = [
    x_a * x_b * x_a^(-1) * x_b^(-1) * x_2ab^(-1) * x_ab^(-1),
    x_a * x_ab * x_a^(-1) * x_ab^(-1) * x_2ab^(-2),
    x_a * x_2ab * x_a^(-1) * x_2ab^(-1),
    x_b * x_ab * x_b^(-1) * x_ab^(-1),
    x_b * x_2ab * x_b^(-1) * x_2ab^(-1),
    x_ab * x_2ab * x_ab^(-1) * x_2ab^(-1),
    w_a * w_b^2 * w_a * w_b^(-2),
    w_b * w_a^2 * w_b^(-1) * w_a^(-2),
    (w_a * w_b)^2 * (w_b * w_a)^(-2),
    w_b^4,
    w_a * x_b * w_a^(-1) * x_2ab^(-1),
    w_a * x_2ab * w_a^(-1) * x_b^(-1),
    w_a * x_ab * w_a^(-1) * x_ab,
    w_b * x_a * w_b^(-1) * x_ab^(-1),
    w_b * x_ab * w_b^(-1) * x_a,
    w_b * x_2ab * w_b^(-1) * x_2ab^(-1),
    w_a * x_a * w_a^(-1) *x_a * w_a * x_a,
    w_b * x_b * w_b^(-1) * x_b * w_b * x_b
]

# for r in Behr_relations
#     @assert quotient_hom_Behr(r) == one(sp4)
# end

jacobian_matrix_Behr = LowCohomologySOS.jacobian_matrix(Behr_relations)

RF_n = parent(first(jacobian_matrix_Behr))
i = RF_n(RF_n.basis[100]).coeffs
RF_n.basis[100]
SparseArrays.nonzeroinds(i)

support_jacobian = union!(
    [one(sp4)],
    [x_a_group, x_b_group, x_ab_group, x_2ab_group, w_a_group, w_b_group],
    inv.([x_a_group, x_b_group, x_ab_group, x_2ab_group, w_a_group, w_b_group])
)

for entry in jacobian_matrix_Behr
    for i in SparseArrays.nonzeroinds(entry.coeffs)
        push!(support_jacobian, quotient_hom_Behr(RF_n.basis[i]))
    end
end

support_jacobian = unique!([support_jacobian; inv.(support_jacobian)])

Δ₁, I = LowCohomologySOS.spectral_gap_elements(
    quotient_hom_Behr,
    Behr_relations,
    support_jacobian
)

sos_problem_Behr = LowCohomologySOS.sos_problem(
    Δ₁,
    I
)

function scs_opt(;
    accel = 10,
    alpha = 1.5,
    eps = 1e-9,
    max_iters = 6_000,
    verbose = true,
)
    return JuMP.optimizer_with_attributes(
        SCS.Optimizer,
        "acceleration_lookback" => accel,
        "acceleration_interval" => max(abs(accel), 1),
        "alpha" => alpha,
        "eps_abs" => eps,
        "eps_rel" => eps,
        "linear_solver" => SCS.DirectSolver,
        "max_iters" => max_iters,
        "warm_start" => true,
        "verbose" => verbose,
    )
end

JuMP.set_optimizer(sos_problem_Behr, scs_opt(eps = 1e-7, max_iters = 15000))

JuMP.optimize!(sos_problem_Behr)

λ, Q =  LowCohomologySOS.get_solution(sos_problem_Behr)

result_bool, result = LowCohomologySOS.certify_sos_decomposition(
    Δ₁,
    I,
    λ,
    Q,
    support_jacobian
)

result