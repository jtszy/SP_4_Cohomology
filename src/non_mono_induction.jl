function non_mono_laplacians(
    G, # either SL(n,ℤ) or SAut(Fₙ)
    half_basis,
    S; # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
    sq_adj_ = "adj",
    twist_coeffs = true
)
    N = div(size(MatrixGroups.matrix_repr(gens(G)[1]))[1],2)

    F_G = FreeGroup(alphabet(G))
    quotient_hom = let source = F_G, target = G
        Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
    end

    # # check if the quotient homomorphism is defined properly
    # @assert length(gens(F_G)) == length(S)
    # for i in eachindex(S)
    #     @assert quotient_hom(gens(F_G,i)) == S[i]
    #     @assert quotient_hom(gens(F_G,i)^(-1)) == S[i]^(-1)
    # end

    relationsx = non_mono_relations(F_G, S, N, sq_adj_)

    embedding = non_mono_gens_embedding(S)

    S_FG = []

    for i in embedding
        push!(S_FG, gens(F_G, i))
    end

    return LowCohomologySOS.spectral_gap_elements(quotient_hom, relationsx, half_basis, S = S_FG, twist_coeffs = twist_coeffs)
end

function non_mono_gens_embedding(
    S # subset of group generators in non-changed order
)
    G = parent(S[1])
    G_gen = gens(G)

    res = []
    ind = 1
    diff = 0
    

    while ind <= length(S)
        if S[ind] == G_gen[ind + diff]
            push!(res, ind + diff)
            ind += 1
        else
            diff += 1
        end
    end

    return res

end

com = SP_4_Cohomology.com

function non_mono_relations(
    F_G::Groups.FreeGroup,
    S, # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
    N::Integer,
    sq_adj_ = "all"
)
    embedding = non_mono_gens_embedding(S)
    gen_dict = Dict(LowCohomologySOS.determine_letter(S[i]) => gens(F_G, embedding[i]) for i in eachindex(S))

    range_as_list = [i for i in 1:N]
    pairs = [(i,j) for i ∈ 1:N for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    triples = [(i,j,k) for i ∈ range_as_list
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]
    
    x(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:A,i,j)]
    y(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,max(i,j),min(i,j) + N)]
    yt(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,min(i,j) + N,max(i,j))]
    z(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i,i + N)]
    zt(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i + N,i)]

    relations_adj = vcat(
        [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),y(j,k))*y(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),yt(i,k))*yt(j,k) for (i,j,k) in triples]                
    )

    relations_sq = vcat(
        # [com(x(i,j),y(i,j))*z(i)^(-2) for (i,j) in pairs],
        # [com(x(i,j),yt(i,j))*zt(j)^2 for (i,j) in pairs],
        # [com(x(i,j),z(j))*(z(i)*y(i,j))^(-1) for (i,j) in pairs],
        # [com(z(i),y(i,j)) for (i,j) in pairs],
        # [com(x(i,j),zt(i))*yt(i,j)*zt(j)^(-1) for (i,j) in pairs],
        # [com(zt(j),yt(i,j)^(-1)) for (i,j) in pairs],
        # [com(y(i,j),zt(i))*z(j)*x(j,i)^(-1) for (i,j) in pairs],
        # [com(x(j,i),zt(j)^(-1)) for (i,j) in pairs],
        # [com(yt(i,j), z(i))*zt(j)*x(i,j) for (i,j) in pairs],
        # [com(x(i,j),z(j)) for (i,j) in pairs]
    )

    if sq_adj_ == "sq" 
        return relations_sq
    elseif sq_adj_ == "adj"
        return relations_adj
    end

    return vcat(relations_sq, relations_adj)
end

function non_mono_embed_matrix(
    M::AbstractMatrix{<:AlgebraElement},
    i::Groups.Homomorphism, # i is intended to be an embedding
    RG_prime::StarAlgebra, # we must provide a suitable underlying group ring (it has to be the group ring of i.target)
    S,
    S_prime
)
    G = i.source

    @assert size(M) == (length(S), length(S))

    S_idies = Dict(S[i] => i for i in eachindex(S))
    S_prime_idies = Dict(S_prime[i] => i for i in eachindex(S_prime))

    RG = parent(first(M))
    @assert all(x -> parent(x) === RG, M)
    @assert G == parent(first(basis(RG)))

    result = [zero(RG_prime) for i in eachindex(S_prime), j in eachindex(S_prime)]

    for s in S
        for t in S
            result[S_prime_idies[i(s)],S_prime_idies[i(t)]] = LowCohomologySOS.embed(i, M[S_idies[s],S_idies[t]], RG_prime)
        end
    end

    return result
end