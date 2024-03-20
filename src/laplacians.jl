function LowCohomologySOS.laplacians(
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

    # check if the quotient homomorphism is defined properly
    @assert length(gens(F_G)) == length(S)
    for i in eachindex(S)
        @assert quotient_hom(gens(F_G,i)) == S[i]
        @assert quotient_hom(gens(F_G,i)^(-1)) == S[i]^(-1)
    end

    relationsx = LowCohomologySOS.relations(F_G, S, N, sq_adj_)
    print(relationsx)
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
    sq_pairs_mono = []
    sq_pairs_mix = []
    sq_pairs_double = []
    adj_pairs_mix = []
    adj_pairs_double = []
    op_pairs = []
    A = alphabet(Sp2N)
    for s in eachindex(S)
        for t in eachindex(S)
            s_i, s_j = mod(A[word(S[s])[1]].i,N), mod(A[word(S[s])[1]].j,N)
            t_i, t_j = mod(A[word(S[t])[1]].i,N), mod(A[word(S[t])[1]].j,N)
            if sort([s_i,s_j]) == sort([t_i, t_j])
                if s_i == s_j
                    push!(mono_pairs,(s,t)) #Mono
                else
                    push!(sq_pairs_double,(s,t)) #Sq-double
                end
            elseif length(intersect!([s_i,s_j],[t_i,t_j])) == 1
                if s_i == s_j || t_i == t_j
                    push!(sq_pairs_mix,(s,t)) #Sq-mixed
                else
                    push!(adj_pairs_double,(s,t)) #Adj-double
                end
            elseif s_i == s_j && t_i == t_j
                push!(sq_pairs_mono,(s,t)) #Sq-mono
            elseif s_i == s_j || t_i == t_j
                push!(adj_pairs_mix,(s,t)) #Adj-mixed
            elseif s_i != s_j && t_i != t_j
                push!(op_pairs,(s,t))
            end
        end
    end
    mono = [(i,j) in mono_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    sq_mono = [(i,j) in sq_pairs_mono ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    sq_mix = [(i,j) in sq_pairs_mix ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    sq_double = [(i,j) in sq_pairs_double ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    adj_double = [(i,j) in adj_pairs_double ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    adj_mix = [(i,j) in adj_pairs_mix ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    op = [(i,j) in op_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]

    @assert mono + sq_mono + sq_mix + sq_double + adj_double + adj_mix + op == Δ₁⁻

    return mono, sq_mono, sq_mix, sq_double, adj_mix, adj_double, op
end

function LowCohomologySOS.relations(
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
        [com(x(i,j),y(i,j))*z(i)^(-2) for (i,j) in pairs],
        [com(x(i,j),yt(i,j))*zt(j)^2 for (i,j) in pairs],
        [com(x(i,j),z(j))*(z(i)*y(i,j))^(-1) for (i,j) in pairs],
        [com(z(i),y(i,j)) for (i,j) in pairs],
        [com(x(i,j),zt(i))*yt(i,j)*zt(j)^(-1) for (i,j) in pairs],
        [com(zt(j),yt(i,j)^(-1)) for (i,j) in pairs],
        [com(y(i,j),zt(i))*z(j)*x(j,i)^(-1) for (i,j) in pairs],
        [com(x(j,i),zt(j)^(-1)) for (i,j) in pairs],
        [com(yt(i,j), z(i))*zt(j)*x(i,j) for (i,j) in pairs],
        [com(x(i,j),z(j)) for (i,j) in pairs]
    )

    if sq_adj_ == "sq" 
        return relations_sq
    elseif sq_adj_ == "adj"
        return relations_adj
    elseif sq_adj_ == "all"
        return vcat(relations_sq, relations_adj)
    end
end
