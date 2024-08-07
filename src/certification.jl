# We provide more subtle version of certification - it turns out that one can 
# add the maximum ℓ₁-norm among all rows of the matrix instead of the whole
# ℓ₁-norm of the matrix.
function LowCohomologySOS.certify_sos_decomposition(
    X,
    order_unit,
    λ::Number,
    Q::AbstractMatrix,
    support,
)
    λ_interval = @interval(λ)
    eoi = LowCohomologySOS._eoi(X, λ_interval, order_unit)
    
    residual = eoi - LowCohomologySOS.sos_from_matrix(parent(first(X)), Q, support)
    max_norm = 0
    for i in 1:size(residual)[1]
        current_norm = @interval(0)
        for j in 1:size(residual)[2]
            current_norm += norm(residual[i,j],1)
        end
        if max_norm < current_norm.hi
            max_norm = current_norm.hi
        end
    end
    l1_norm = @interval(max_norm)
    
    @info "l₁ norm of the error in interval arithmetic:" l1_norm radius(l1_norm)
    
    result = λ_interval - l1_norm
    
    @info "Is a positive spectral gap certidied? ", result.lo > 0
    @info "Cetified value: ", result.lo
end
