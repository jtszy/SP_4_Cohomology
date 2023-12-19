#aux functions for sp2n_sp2m_embedding

function ind_embed((type, i, j), N, M)
    if type == :A
        return (:A, i, j)

    elseif i < j
        return (:B, i, j + Int8((M - N)/2))

    else
        return (:B, i + Int8((M - N)/2), j)
    
    end

end

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
           end => (:A, i, j)
           for (i,j) in ind_A(Int8(N/2))
        )

        res_ = Dict(
            let Bij = MatrixGroups.ElementarySymplectic{N}(:B, i, j)
                Sp_N([alphabet(Sp_N)[Bij]])
            end => (:B, i, j)
            for (i,j) in ind_B(Int8(N/2))
        )

        merge!(res, res_)
        res

    end

    S_Sp_M = let 
        res = Dict((:A,i,j) => 
        let Aij = MatrixGroups.ElementarySymplectic{M}(:A, i, j)
            Sp_M([alphabet(Sp_M)[Aij]])

        end
        for (i,j) in ind_A(Int8(M/2))
        )

        res_ =Dict((:B,i,j) =>
        let Bij = MatrixGroups.ElementarySymplectic{M}(:B, i , j)
            Sp_M([alphabet(Sp_M)[Bij]])

        end
        for (i,j) in ind_B(Int8(M/2))
        )

        merge!(res, res_)
        res

    end    

    function f(letter_id, Sp_N, G)
        if letter_id <= length(gens(Sp_N))
            N_elem = Sp_N([letter_id])
            N_ind = inds_S_Sp_N[N_elem]
            M_ind = ind_embed(N_ind, N, M)
            M_elem = S_Sp_M[M_ind]
            return word(M_elem)

        else 
            N_elem = inv(Sp_N([letter_id]))
            N_ind = inds_S_Sp_N[N_elem]
            M_ind = ind_embed(N_ind, N, M)
            M_elem = inv(S_Sp_M[M_ind])
            return word(M_elem)

        end

    end

    result = let source = Sp_N, target = Sp_M
        Groups.Homomorphism(f, source, target, check = false)
    end

    return result
    
end

SP_2N = MatrixGroups.SymplecticGroup{6}(Int8)


A = alphabet(SP_2N)

S = gens(SP_2N)

length(S)

A[word(S[1])[1]].j



function LowCohomologySOS.laplacians(
    G, # either SL(n,ℤ) or SAut(Fₙ)
    half_basis,
    S; # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
    twist_coeffs = true,
    sq_adj_ = "adj"
)
    N = div(size(MatrixGroups.matrix_repr(gens(G)[1]))[1],2)

    F_G = FreeGroup(alphabet(G))
    quotient_hom = let source = F_G, target = G
        Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
    end

    # check if the quotient homomorphism is defined properly
    @assert length(gens(F_G)) == length(S)
    for i in eachindex(S)
        @assert quotient_hom(gens(F_G,i)) == S[i]
        @assert quotient_hom(gens(F_G,i)^(-1)) == S[i]^(-1)
    end

    relationsx = LowCohomologySOS.relations(G, F_G, S, N, sq_adj_)
    return LowCohomologySOS.spectral_gap_elements(quotient_hom, relationsx, half_basis, twist_coeffs = twist_coeffs)
end

function mono_sq_adj_op(
    Δ₁⁻,
    S # generating set indexing Δ₁⁻
)
    RG = parent(first(Δ₁⁻))
    Sp2N = parent(first(RG.basis))
    N = Int8(sqrt(length(gens(Sp2N))/2))
    mono_pairs = []
    sq_pairs = []
    adj_pairs = []
    op_pairs = []
    A = alphabet(Sp2N)
    for s in eachindex(S)
        for t in eachindex(S)
            s_i, s_j = mod(A[word(S[s])[1]].i,N), mod(A[word(S[s])[1]].j,N)
            t_i, t_j = mod(A[word(S[t])[1]].i,N), mod(A[word(S[t])[1]].j,N)
            if sort([s_i,s_j]) == sort([t_i, t_j])
                if s_i == s_j
                    push!(mono_pairs,(s,t))
                else
                    push!(sq_pairs,(s,t))
                end
            elseif length(intersect!([s_i,s_j],[t_i,t_j])) == 1
                push!(adj_pairs,(s,t))
            else
                push!(op_pairs,(s,t))
            end
        end
    end
    mono = [(i,j) in mono_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    sq = [(i,j) in sq_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    adj = [(i,j) in adj_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    op = [(i,j) in op_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]

    @assert mono+sq+adj+op == Δ₁⁻

    return mono, sq, adj, op
end

function elementary_conjugation(
    # l,
    l::MatrixGroups.ElementarySymplectic,
    σ::PermutationGroups.AbstractPerm 
)
    l_matrix = MatrixGroups.matrix_repr(l)
    N = size(l_matrix)[1]
    n = div(N,2)

    type, l_i, l_j, l_transpose = elementary_element_info(l_matrix)

    if type == :A
        # @info l_matrix
        # @info l_i, l_i^σ, l_j, l_j^σ
        return Groups.MatrixGroups.ElementarySymplectic{N}(type, l_i^σ, l_j^σ)
    end

    return Groups.MatrixGroups.ElementarySymplectic{N}(type, l_i^σ + (1 - l_transpose)*n, l_j^σ + l_transpose*n)
end

function conjugation(
    # l,
    l::MatrixGroups.ElementarySymplectic,
    σ::PermutationGroups.AbstractPerm
)
    l_matrix = MatrixGroups.matrix_repr(l)
    N = size(l_matrix)[1]
    n = div(N,2)

    for i in 1:N
        for j in 1:N
            if l_matrix[i,j] < 0 && (i <= n || j <= n)
                return elementary_conjugation(l^(-1), σ)^(-1)
            end
        end
    end

    return elementary_conjugation(l,σ)
    
end

function LowCohomologySOS.relations(
    G,
    F_G::Groups.FreeGroup,
    S, # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
    N::Integer,
    sq_adj_ = "all"
)
    gen_dict = Dict(LowCohomologySOS.determine_letter(S[i]) => gens(F_G, i) for i in eachindex(S))

    range_as_list = [i for i in 1:N]
    pairs = [(i,j) for i ∈ 1:N for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    triples = [(i,j,k) for i ∈ range_as_list
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]
    
    x(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:A,i,j)]
    y(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i+N,j)]
    yt(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,j,i+N)]
    z(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i+N,i)]
    zt(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i,i+N)]

    relations_adj = vcat(
        [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),y(min(j,k),max(j,k)))*y(min(i,k),max(i,k))^(-1) for (i,j,k) in triples],
        [com(x(i,j),yt(min(i,k),max(i,k)))*yt(min(j,k),max(j,k)) for (i,j,k) in triples]                
    )

    relations_sq = vcat(
        [com(x(i,j),y(min(i,j),max(i,j)))*z(i)^(-2) for (i,j) in pairs],
        [com(x(i,j),yt(min(i,j),max(i,j)))*zt(j)^2 for (i,j) in pairs],
        [com(x(i,j),z(j))*(z(i)*y(min(i,j),max(i,j)))^(-1) for (i,j) in pairs],
        [com(z(i),y(min(i,j),max(i,j))) for (i,j) in pairs],
        [com(x(i,j),zt(i))*yt(min(i,j),max(i,j))*zt(j)^(-1) for (i,j) in pairs],
        [com(zt(j),yt(min(i,j),max(i,j))^(-1)) for (i,j) in pairs],
        [com(y(min(i,j),max(i,j)),zt(i))*z(j)*x(j,i)^(-1) for (i,j) in pairs],
        [com(x(j,i),zt(j)^(-1)) for (i,j) in pairs],
        [com(yt(min(i,j),max(i,j)), z(i))*zt(j)*x(i,j) for (i,j) in pairs],
        [com(x(i,j),z(j)) for (i,j) in pairs]
    )

    if sq_adj_ == "sq" 
        return relations_sq
    elseif sq_adj_ == "adj"
        return relations_adj
    end

    return vcat(relations_sq, relations_adj)
end

N = 4
range_as_list = [i for i in 1:N]
pairs = [(i,j) for i ∈ 1:N for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]


