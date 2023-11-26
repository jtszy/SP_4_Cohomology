function sp2n_sp2m_embedding(N::Integer,M::Integer)
    @assert N <= M
    @assert mod(N, 2) == 0
    @assert mod(M, 2) == 0

    Sp_N = MatrixGroups.SymplecticGroup{N}(Int8)
    Sp_M = MatrixGroups.SymplecticGroup{M}(Int8)

    ind_A(k) = [(i,j) for i in 1:k for j in 1:k if i != j]

    function ind_B(k)
        res = [(i,i+k) for i in 1:k]
        res = [res; [(i+k, i) for i in 1:k]]
        res = [res; [(i+k, j) for i in 1:k for j in 1:k if i<j]]
        res = [res; [(j, i+k) for i in 1:k for j in 1:k if i<j]]
        return res

    end

    inds_S_Sp_N = let 
        res = Dict(
           let Aij = MatrixGroups.ElementarySymplectic{N}(:A, i, j)
               Sp_N([alphabet(Sp_N)[Aij]])
           end => ("A", i, j)
           for (i,j) in ind_A(Int8(N/2))
        )

        res_ = Dict(
            let Bij = MatrixGroups.ElementarySymplectic{N}(:B, i, j)
                Sp_N([alphabet(Sp_N)[Bij]])
            end => ("B", i, j)
            for (i,j) in ind_B(Int8(N/2))
        )

        merge!(res, res_)
        res

    end

    S_Sp_M = let 
        res = Dict(("A",i,j) => 
        let Aij = MatrixGroups.ElementarySymplectic{M}(:A, i, j)
            Sp_M([alphabet(Sp_M)[Aij]])

        end
        for (i,j) in ind_A(Int8(M/2))
        )

        res_ =Dict(("B",i,j) =>
        let Bij = MatrixGroups.ElementarySymplectic{M}(:B, i , j)
            Sp_M([alphabet(Sp_M)[Bij]])

        end
        for (i,j) in ind_B(Int8(M/2))
        )

        merge!(res, res_)
        res

    end

    function ind_embed((type, i, j))
            if type == "A"
                return ("A", i, j)

            elseif i < j
                return ("B", i, j + Int8((M - N)/2))

            else
                return ("B", i + Int8((M - N)/2), j)
            
            end

    end    

    function f(letter_id, Sp_N, G)
        if letter_id <= length(gens(Sp_N))
            N_elem = Sp_N([letter_id])
            N_ind = inds_S_Sp_N[N_elem]
            M_ind = ind_embed(N_ind)
            M_elem = S_Sp_M[M_ind]
            return word(M_elem)

        else 
            N_elem = inv(Sp_N([letter_id]))
            N_ind = inds_S_Sp_N[N_elem]
            M_ind = ind_embed(N_ind)
            M_elem = inv(S_Sp_M[M_ind])
            return word(M_elem)

        end

    end

    result = let source = Sp_N, target = Sp_M
        Groups.Homomorphism(f, source, target, check = false)
    end

    return result
 
end

# function LowCohomologySOS.laplacians(
#     G, # either SL(n,ℤ) or SAut(Fₙ)
#     half_basis,
#     S; # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
#     twist_coeffs = true,
#     sq_adj_op_ = "all"
# )
#     # TODO
# end

# function LowCohomologySOS.sq_adj_op(matrix, gen_set)
#     # TODO
# end

# function LowCohomologySOS._conj(
#     t::fill in the type,
#     σ::PermutationGroups.AbstractPerm,
# )
#     # TODO
# end