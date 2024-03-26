module SP_4_Cohomology

using SparseArrays
using Groups
using JuMP
using SCS
using LowCohomologySOS
using LinearAlgebra
using PermutationGroups
using StarAlgebras
using IntervalArithmetic

Base.adjoint(X::AlgebraElement) = StarAlgebras.star(X)
StarAlgebras.star(g::Groups.AbstractFPGroupElement) = inv(g)

include("certification.jl")
include("common_functions.jl")
include("optimizers.jl")
include("relations.jl")
include("wedderburn.jl")

end