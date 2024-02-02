using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
using SparseArrays
using JuMP
using SCS
using IntervalArithmetic
using Serialization
using PermutationGroups
using SymbolicWedderburn
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Groups
using LowCohomologySOS

using SP_4_Cohomology

N = 3
Sp_N = MatrixGroups.SymplecticGroup{2*N}(Int8)
S = gens(Sp_N)

F_Sp_N_Steinberg = FreeGroup(2*N^2)
Gens = gens(F_Sp_N_Steinberg)

id = Gens[1] * Gens[1] ^ -1

X = [id for i in 1:N, j in 1:N]
Y = [id for i in 1:N, j in 1:N]
Y_t = [id for i in 1:N, j in 1:N]
Z = [id for i in 1:N]
Z_t = [id for i in 1:N]

# Assignment of generators to matrices X, Y, Y_t, Z, Z_t for code simplicity
iter = 1

for i in 1:N
    for j in 1:N
        if i != j
            X[i,j] = Gens[iter]
            iter += 1
            # for limiting the number of generators we consider symmetric Y and Y_t
            # it can be done thanks to relations y_i_j = y_j_i and y_t_i_j = y_t_j_i
            if i < j
                Y[i,j] = Gens[iter]
                iter += 1
                Y_t[i,j] = Gens[iter]
                iter += 1
            else
                Y[i,j] = Y[j,i]
                Y_t[i,j] = Y_t[j,i]
            end
        end
    end
    Z[i] = Gens[iter]
    iter += 1
    Z_t[i] = Gens[iter]
    iter += 1
end

Steinberg_matrix = Dict()

function E(i,j)
    res = [0 for i in 1:2*N, j in 1:2*N]
    res[i,j] = 1
    return res
end

for i in 1:N
    for j in 1:N
        if i != j
            Steinberg_matrix[X[i,j]] = 1*Matrix(I, 2*N, 2*N) + E(i,j) - E(j+N,i+N)
            if i < j
                # Y and Y_t are symmetric so we conduct assignment for upper halfs only
                Steinberg_matrix[Y[i,j]] = 1*Matrix(I, 2*N, 2*N) + E(i,j+N) + E(j,i+N)
                Steinberg_matrix[Y_t[i,j]] = transpose(Steinberg_matrix[Y[i,j]])
            end
        end
    end
    Steinberg_matrix[Z[i]] = 1*Matrix(I, 2*N, 2*N) + E(i,i+N)
    Steinberg_matrix[Z_t[i]] = transpose(Steinberg_matrix[Z[i]])
end


gens_Sp_N = [S; inv.(S)]

Ball_4, sizes = Groups.wlmetric_ball(gens_Sp_N, radius = 3)

Steinberg_group = Dict()

for g in Ball_4
    for gen in keys(Steinberg_matrix)
        if MatrixGroups.matrix_repr(g) == Steinberg_matrix[gen]
            Steinberg_group[gen] = g
        end
    end
end

hom_Steinberg = SP_4_Cohomology.quotient_homomorphism(F_Sp_N_Steinberg, Sp_N, Steinberg_group)

for gen in keys(Steinberg_matrix)
    @assert MatrixGroups.matrix_repr(hom_Steinberg(gen)) == Steinberg_matrix[gen]
end

Steinberg_relations = []

function comm(x,y)
    x*y*x^(-1)*y^(-1)
end

for i in 1:N
    for j in 1:N
        for k in 1:N
            if i != j && j != k && i != k
                #Adj relatons
                push!(Steinberg_relations, comm(X[i,j],X[j,k]) * X[i,k]^(-1))
                push!(Steinberg_relations, comm(X[i,j],Y[j,k]) * Y[i,k]^(-1))
                push!(Steinberg_relations, comm(X[i,j],Y_t[i,k]) * Y_t[j,k])
            end
        end
        if i != j
            #Sq relations

            # we dont need these two becouse Y and Y_t are already symmetric
            # push!(Steinberg_relations, Y[i,j]*Y[j,i]^(-1))
            # push!(Steinberg_relations, Y_t[i,j]*Y_t[j,i]^(-1))

            push!(Steinberg_relations, comm(X[i,j],Y[i,j]) * Z[i] ^ (-2))
            push!(Steinberg_relations, comm(X[i,j],Y_t[i,j]) * Z_t[j] ^ 2)
            push!(Steinberg_relations, comm(X[i,j],Z[j]) * (Y[i,j] * Z[i]) ^ (-1))
            push!(Steinberg_relations, comm(Z[i],Y[i,j]))
            push!(Steinberg_relations, comm(X[i,j],Z_t[i]) * Z_t[j] ^ (-1) * Y_t[i,j])
            push!(Steinberg_relations, comm(Z_t[j],Y_t[i,j]^(-1)))
            push!(Steinberg_relations, comm(Y[i,j],Z_t[i]) * X[j,i] ^ (-1) * Z[j])
            push!(Steinberg_relations, comm(X[j,i],Z[j]^(-1)))
            push!(Steinberg_relations, comm(Y_t[i,j],Z[i]) * X[i,j] * Z_t[j])
            push!(Steinberg_relations, comm(X[i,j],Z_t[j]))
        end
    end
end

for r in Steinberg_relations
    @info r, hom_Steinberg(r)
    @assert hom_Steinberg(r) == one(Sp_N)
end

jacobian_matrix_Steinberg = LowCohomologySOS.jacobian_matrix(Steinberg_relations)

min_support = SP_4_Cohomology.minimalistic_support(jacobian_matrix_Steinberg, hom_Steinberg, Steinberg_group)

function wedderburn_data(basis, half_basis, S)
    @time begin
        N = div(size(first(S))[1], 2)
        
        # Z_2_wr_S(n) = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(n))
        # Σ = Z_2_wr_S(N)

        Σ = PermutationGroups.SymmetricGroup(N)

        actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, SP_4_Cohomology.conjugation, S, basis)
        constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)
    end

    return constraints_basis, psd_basis, Σ, actions
end

# Δ₁, I_N = LowCohomologySOS.spectral_gap_elements(
#     hom_Steinberg,
#     Steinberg_relations,
#     support_jacobian,
#     # twist_coeffs = false
# )

half_basis, sizes = Groups.wlmetric_ball(gens_Sp_N, radius = 3);

Δ₁, I_N, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.laplacians(Sp_N, half_basis, S, sq_adj_ = "all");

RG = LowCohomologySOS.group_ring(Sp_N, min_support, star_multiplication = true)

Δ₁ = LowCohomologySOS.embed.(identity, Δ₁, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N, Ref(RG))
Δ₁⁺ = LowCohomologySOS.embed.(identity, Δ₁⁺, Ref(RG))
Δ₁⁻ = LowCohomologySOS.embed.(identity, Δ₁⁻, Ref(RG))

basis = RG.basis

constraints_basis, psd_basis, Σ, action = wedderburn_data(basis, support_jacobian, S);

M = Δ₁

# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    @info σ
    @assert LowCohomologySOS.act_on_matrix(M, σ, action.alphabet_perm, S) == M
    @assert LowCohomologySOS.act_on_matrix(I_N, σ, action.alphabet_perm, S) == I_N
end

A = LowCohomologySOS.act_on_matrix(M, collect(Σ)[2], action.alphabet_perm, S)

@time begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

@time begin
    sos_problem, P = LowCohomologySOS.sos_problem(
        M,
        I_N,
        w_dec_matrix
        # 0.7 / 0.05
    )
end 

# Find a numerical spectral gap
JuMP.set_optimizer(sos_problem, SP_4_Cohomology.scs_opt(eps = 1e-5, max_iters = 20000))
JuMP.optimize!(sos_problem)

# Certify the numerical estimate
λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)
LowCohomologySOS.certify_sos_decomposition(M, I, λ, Q, support_jacobian)

# JuMP.set_optimizer(sos_problem_Steinberg, SP_4_Cohomology.scs_opt(eps = 1e-7, max_iters = 20000))

# JuMP.optimize!(sos_problem_Steinberg)

# λ, Q =  LowCohomologySOS.get_solution(sos_problem_Steinberg)

# result_bool, result = LowCohomologySOS.certify_sos_decomposition(
#     Δ₁,
#     I_N,
#     λ,
#     Q,
#     support_jacobian
# )

# Solution = Dict("lambda" => λ, "Q" => Q)

# serialize("./Steinberg_Solution.sjl", Solution)

# Solution_read = deserialize("./Steinberg_Solution.sjl")

# LowCohomologySOS.certify_sos_decomposition(
#     Δ₁,
#     I_4,
#     Solution_read["lambda"],
#     Solution_read["Q"],
#     support_jacobian
# )