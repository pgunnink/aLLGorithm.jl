struct LLGEulerHeun <: OrdinaryLLGAlgorithm end

isfsal(alg::LLGEulerHeun) = false
alg_order(alg::LLGEulerHeun) = 1


@cache struct EulerHeunCache{uType} <: OrdinaryDiffEqMutableCache
    up::uType
    Acache::uType
    N::Int64
end


@fastmath @muladd function OrdinaryDiffEq.perform_step!(integrator, cache::EulerHeunCache, repeat_step = false)
    @unpack Acache, up, N = cache
    f(Acache, up, t)
    @. up = u + integrator.dt
end
