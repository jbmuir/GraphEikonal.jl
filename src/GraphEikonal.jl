module GraphEikonal

export graph_eikonal
export L1Norm, L2Norm, LInfNorm

using LightGraphs: AbstractGraph, outneighbors, nv, weights
using DataStructures: PriorityQueue, dequeue!

abstract type LPNorm end
struct L1Norm <: LPNorm end
struct L2Norm <: LPNorm end
struct LInfNorm <: LPNorm end

function eik_solve(v::Integer, 
          g::AbstractGraph,
          dists::Vector{T},
          finalized::Vector{Bool},
          distmx::AbstractMatrix{T},
          slowv::Vector{T},
          norm::L1Norm) where T <: Real

    distv = 1/zero(T)
    w = outneighbors(g, v)
    w = w[finalized[w].==true] #note that this is not in desquesnes but seems like it should be true...only use finalized values when computing solutions, so you can never accidentally dequeue a value that could rely on unfinalized values
    pdw = sortperm(dists[w])
    a = push!(dists[w[pdw]], 1/zero(T))
    w = w[pdw]
    h = 1 ./ sqrt.(distmx[v, w])
    c = slowv[v]

    for m in 1:length(w)
        distv = (sum(a[i]/h[i] for i in 1:m) + c)/sum(1/h[i] for i in 1:m)
        if distv <= a[m+1]
            break
        end
    end

    return distv
end

function eik_solve(v::Integer, 
          g::AbstractGraph,
          dists::Vector{T},
          finalized::Vector{Bool},
          distmx::AbstractMatrix{T},
          slowv::Vector{T},
          norm::L2Norm) where T <: Real
        
    distv = 1/zero(T)
    w = outneighbors(g, v)
    w = w[finalized[w].==true] #not that this is not in desquesnes but seems like it should be true...only use finalized values when computing solutions, so you can never accidentally dequeue a value that could rely on unfinalized values
    pdw = sortperm(dists[w])
    a = push!(dists[w[pdw]], 1/zero(T))
    w = w[pdw]
    h = 1 ./ sqrt.(distmx[v, w])
    c = slowv[v]
      

    for m in 1:length(w) # desquesnes formula for this never seems to work (ie. i guess I can't code properly but whatver)
        H1 = sum(1/h[i]^2 for i in 1:m)
        H2 = sum(a[i]/h[i]^2 for i in 1:m)
        H3 = sum(a[i]^2/h[i]^2 for i in 1:m)
        distv = (H2 + sqrt(H1*(c^2-H3)+H2^2))/H1
        if distv <= a[m+1]
            break
        end
    end

    return distv
end


function eik_solve(v::Integer, 
          g::AbstractGraph,
          dists::Vector{T},
          finalized::Vector{Bool},
          distmx::AbstractMatrix{T},
          slowv::Vector{T},
          norm::LInfNorm) where T <: Real

          distv= 1 / zero(T)

          for w in outneighbors(g, v)
            distvp = dists[w] + slowv[v] / sqrt(distmx[v,w])

            if distvp < distv
                distv = distvp
            end

          end

    return distv
end

function graph_eikonal(g::AbstractGraph,
    srcs::Vector{U},
    distmx::AbstractMatrix{T} = weights(g);
    slowv::Vector{T} = ones(T, nv(g)), # as a seismologist, this is the slowness (P in Desquesnes, Elmoataz & Lezoray 2013)
    norm::LPNorm = L2Norm()) where T <: Real where U <: Integer 

    nvg = nv(g)
    dists = fill(typemax(T), nvg)
    visited = zeros(Bool, nvg)
    finalized = zeros(Bool, nvg)
    H = PriorityQueue{U,T}()
    # fill creates only one array.

    for src in srcs
        dists[src] = zero(T)
        visited[src] = true
        H[src] = zero(T)
    end

    while !isempty(H)
        u = dequeue!(H)
        finalized[u] = true
        
        for v in outneighbors(g, u)
            if !finalized[v]
                alt = eik_solve(v, g, dists, finalized, distmx, slowv, norm)

                if !visited[v] 
                    visited[v] = true
                    dists[v] = alt
                    H[v] = alt
                elseif alt < dists[v]
                    dists[v] = alt
                    H[v] = alt
                end

            end
            
        end

    end

    return dists
end

end 
