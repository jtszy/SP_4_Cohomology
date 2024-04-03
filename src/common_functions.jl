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

function com(x,y)
    x*y*x^(-1)*y^(-1)
end

function minimalistic_support(relations, quotient_hom)
    jacobian_matrix = LowCohomologySOS.jacobian_matrix(relations)

    RF_G = parent(first(jacobian_matrix))
    Sp_N = quotient_hom.target

    min_support = union!(
        [one(Sp_N)],
        gens(Sp_N),
        inv.(gens(Sp_N))
    )

    for entry in jacobian_matrix
        for i in SparseArrays.nonzeroinds(entry.coeffs)
            push!(min_support, quotient_hom(RF_G.basis[i]))
        end
    end

    min_support = unique!([min_support; inv.(min_support)])
    return min_support
end

function support_jacobian(relations, quotient_hom)
    Sp_N = quotient_hom.target
    F_G = quotient_hom.source

    support_jacobian = [one(Sp_N)]
    for r in relations
        current_factor = one(Sp_N)
        for i in 1:length(word(r))
            g = quotient_hom(F_G(word(r)[i:i]))
            push!(support_jacobian,current_factor*g)
            push!(support_jacobian,current_factor*g^(-1))
            current_factor *= g
        end
    end
    support_jacobian = unique(support_jacobian)

    return support_jacobian
end

function symplectic_min_supports(
    quotient_hom,
    S
)   
    Sp_N = quotient_hom.target
    F_Sp_N_Steinberg = quotient_hom.source

    N = div(size(MatrixGroups.matrix_repr(first(S)))[1],2)

    Steinberg_relations = relations_St(F_Sp_N_Steinberg, S, N)

    for r in Steinberg_relations
        @assert quotient_hom(r) == one(Sp_N)
    end 

    sup_jacobian = support_jacobian(Steinberg_relations, quotient_hom)

    min_support = minimalistic_support(Steinberg_relations, quotient_hom)

    return sup_jacobian, min_support

end