Base.adjoint(X::AlgebraElement) = StarAlgebras.star(X)
StarAlgebras.star(g::Groups.AbstractFPGroupElement) = inv(g)

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

function elementary_element_info(l_matrix)
    N = size(l_matrix)[1]
    n = div(N,2)
    type = :B
    l_transpose = 1
    l_i, l_j = -1, -1

    for i in 1:N
        for j in 1:N

            if l_matrix[i,j] < 0
                type = :A
            end

            if l_matrix[i,j] == 1 && i != j
                l_i = i % n
                l_j = j % n

                l_i = (l_i == 0) ? n : l_i
                l_j = (l_j == 0) ? n : l_j

                if i < j
                    l_transpose = 0
                end
            end
        end
    end

    return (type, l_i, l_j, l_transpose)
end