struct DepondtOrdinary <: OrdinaryLLGAlgorithm end
struct DepondtStochastic <: StochasticLLGAlgorithm end
isfsal(alg::DepondtOrdinary) = false
alg_order(alg::SIBOrdinary) = 1


@cache struct DepondtCache{uType} <: OrdinaryDiffEqMutableCache
    up::uType
    Acache::uType
    N::Int64
    αkbT::Float64
    detM::Vector{Float64}
    detMx::Vector{Float64}
    detMy::Vector{Float64}
    detMz::Vector{Float64}
end

@fastmath @muladd function OrdinaryDiffEq.perform_step!(integrator, cache::DepondtCache, repeat_step = false)
    @unpack t, dt, u, f, p = integrator
    @unpack Acache, N, up = cache

    

end