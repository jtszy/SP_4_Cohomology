function sp2n_sp2m_embedding(N::Integer,M::Integer)
    # TODO
end

function LowCohomologySOS.laplacians(
    G, # either SL(n,ℤ) or SAut(Fₙ)
    half_basis,
    S; # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
    twist_coeffs = true,
    sq_adj_op_ = "all"
)
    # TODO
end

function LowCohomologySOS.sq_adj_op(matrix, gen_set)
    # TODO
end

function LowCohomologySOS._conj(
    t::fill in the type,
    σ::PermutationGroups.AbstractPerm,
)
    # TODO
end