# GraphEikonal.jl

Eikonal equations on graphs following Desquenes et al. 2013. The fundamental function is
```
graph_eikonal(g::AbstractGraph,
    srcs::Vector{U},
    distmx::AbstractMatrix{T} = weights(g);
    slowv::Vector{T} = ones(T, nv(g)),
    norm::LPNorm = L2Norm()) where T <: Real where U <: Integer 
```
with three norms supported - L1Norm(), L2Norm(), LInfNorm()