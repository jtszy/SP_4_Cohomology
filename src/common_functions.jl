function quotient_homomorphism(free_group, sp_n, gen_dict)
    function F(i, source, target)
        if source([i]) == one(free_group)
            return one(sp_n)
        end
        for gen in keys(gen_dict)
            if source([i]) == gen
                return word(gen_dict[gen])
            elseif source([i]) == gen^(-1)
                return word(gen_dict[gen]^(-1))
            end
        end
    end
    Groups.Homomorphism(F, free_group, sp_n)
end

function minimalistic_support(jacobian_matrix, quotient_homomorphism, gen_dict)
    RF_n = parent(first(jacobian_matrix))
    F_n = parent(first(RF_n.basis))
    sp_n = parent(quotient_homomorphism(one(F_n)))

    support_jacobian = union!(
        [one(sp_n)],
        values(gen_dict),
        inv.(values(gen_dict))
    )

    for entry in jacobian_matrix
        for i in SparseArrays.nonzeroinds(entry.coeffs)
            push!(support_jacobian, quotient_homomorphism(RF_n.basis[i]))
        end
    end

    support_jacobian = unique!([support_jacobian; inv.(support_jacobian)])
    return support_jacobian
end

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
