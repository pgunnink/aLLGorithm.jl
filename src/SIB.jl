struct SIBOrdinary <: OrdinaryLLGAlgorithm end
struct SIBStochastic <: StochasticLLGAlgorithm end
isfsal(alg::SIBOrdinary) = false
alg_order(alg::SIBOrdinary) = 1


# @cache struct StochasticSIB   Cache{uType} <: OrdinaryDiffEqMutableCache
#     up::uType
#     Acache::uType
#     N::Int64
#     αkbT::Float64
#     detM::Vector{Float64}
#     detMx::Vector{Float64}
#     detMy::Vector{Float64}
#     detMz::Vector{Float64}
# end

@cache struct SIBCache{uType} <: OrdinaryDiffEqMutableCache
    up::uType
    Acache::uType
    N::Int64
    αkbT::Float64
    detM::Vector{Float64}
    detMx::Vector{Float64}
    detMy::Vector{Float64}
    detMz::Vector{Float64}
end

function OrdinaryDiffEq.alg_cache(alg::SIBOrdinary, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt, reltol_internal, p, calck, iip)
    Acache = similar(u)
    up = similar(u)
    N = size(u)[1]
    detM = zeros(Float64, N)
    detMx = zeros(Float64, N)
    detMy = zeros(Float64, N)
    detMz = zeros(Float64, N)
    SIBCache(up, Acache, N, 0.,  detM, detMx, detMy, detMz)
end

function StochasticDiffEq.alg_cache(alg::SIBStochastic, prob::StochasticLLGProblem, u, dW, dZ, p, rate_prototype, noise_rate_prototype, jump_prototype, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, uprev, f, t, dt, iip)
    Acache = similar(u)
    up = similar(u)
    N = size(u)[1]
    detM = zeros(Float64, N)
    detMx = zeros(Float64, N)
    detMy = zeros(Float64, N)
    detMz = zeros(Float64, N)
    SIBCache(up, Acache, N, prob.αkbT, detM, detMx, detMy, detMz)
end
function initialize!(integrator, cache::SIBCache)
    # integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    # integrator.fsallast = zero(integrator.fsalfirst)

end


@fastmath @muladd function OrdinaryDiffEq.perform_step!(integrator, cache::SIBCache, repeat_step = false)
    @unpack t, dt, u, f, p = integrator
    # @unpack M, Mx, My, Mz, a, Acache, N, up  = cache
    @unpack Acache, N, up, detM, detMx, detMy, detMz = cache

    # do the predictor step:
    f(Acache, u, p, t)
    @. Acache = .5dt * (Acache) # + W.dW * √(2α * kbT))
    
    fdet(x11, x21, x31, x12, x22, x32, x13, x23, x33) = x11 * (x22 * x33 - x23 * x32) - x12 * (x21 * x33 - x23 * x31) + x13 * (x21 * x32 - x22 * x31)
    
    combined_cross(uloc, Ax, Ay, Az) = (uloc[1] + uloc[2] * Az - Ay * uloc[3], uloc[2] + uloc[3] * Ax - uloc[1] * Az, uloc[3] + uloc[1] * Ay - uloc[2] * Ax)

    @inbounds for i in 1:N
        uloc = @view u[i,:]
        Ax, Ay, Az = @view Acache[i,:]
        a1, a2, a3 = combined_cross(uloc, Ax, Ay, Az)
        detMx[i] =  fdet(a1, a2, a3, -Az, 1, Ax, Ay, -Ax, 1)
        detMy[i] = fdet(1, Az, -Ay, a1, a2, a3, Ay, -Ax, 1)
        detMz[i] = fdet(1, Az, -Ay, -Az, 1, Ax, a1, a2, a3)
        detM[i] = 1 / fdet(1, Az, -Ay, -Az, 1, Ax, Ay, -Ax, 1)
    end
    @. up[:,1] = detMx * detM
    @. up[:,2] = detMy * detM
    @. up[:,3] = detMz * detM

    # # and the final step:
    @. up = .5 * ( up + u)
    f(Acache,  up, p, t + .5dt) 
    @. Acache = .5dt * (Acache) # + W.dW * √(2α * kbT))

    @inbounds for i in 1:N
        uloc = @view u[i,:]
        Ax, Ay, Az = @view Acache[i,:]
        a1, a2, a3 = combined_cross(uloc, Ax, Ay, Az)
        detMx[i] =  fdet(a1, a2, a3, -Az, 1, Ax, Ay, -Ax, 1)
        detMy[i] = fdet(1, Az, -Ay, a1, a2, a3, Ay, -Ax, 1)
        detMz[i] = fdet(1, Az, -Ay, -Az, 1, Ax, a1, a2, a3)
        detM[i] = 1 / fdet(1, Az, -Ay, -Az, 1, Ax, Ay, -Ax, 1)
    end
    @. u[:,1] = detMx * detM
    @. u[:,2] = detMy * detM
    @. u[:,3] = detMz * detM
end

@fastmath @muladd function StochasticDiffEq.perform_step!(integrator, cache::SIBCache)
    @unpack t, dt, u, f, W, p = integrator
    # @unpack M, Mx, My, Mz, a, Acache, N, up  = cache
    @unpack Acache, N, up, αkbT, detM, detMx, detMy, detMz = cache

    # do the predictor step:
    f(Acache, u, p, t)
    @. Acache = .5dt * (Acache + W.dW * √(2αkbT))
    
    fdet(x11, x21, x31, x12, x22, x32, x13, x23, x33) = x11 * (x22 * x33 - x23 * x32) - x12 * (x21 * x33 - x23 * x31) + x13 * (x21 * x32 - x22 * x31)
    
    combined_cross(uloc, Ax, Ay, Az) = (uloc[1] + uloc[2] * Az - Ay * uloc[3], uloc[2] + uloc[3] * Ax - uloc[1] * Az, uloc[3] + uloc[1] * Ay - uloc[2] * Ax)

    @inbounds for i in 1:N
        uloc = @view u[i,:]
        Ax, Ay, Az = @view Acache[i,:]
        a1, a2, a3 = combined_cross(uloc, Ax, Ay, Az)
        detMx[i] =  fdet(a1, a2, a3, -Az, 1, Ax, Ay, -Ax, 1)
        detMy[i] = fdet(1, Az, -Ay, a1, a2, a3, Ay, -Ax, 1)
        detMz[i] = fdet(1, Az, -Ay, -Az, 1, Ax, a1, a2, a3)
        detM[i] = 1 / fdet(1, Az, -Ay, -Az, 1, Ax, Ay, -Ax, 1)
    end
    @. up[:,1] = detMx * detM
    @. up[:,2] = detMy * detM
    @. up[:,3] = detMz * detM

    # # and the final step:
    @. up = .5 * ( up + u)
    f(Acache,  up, p, t + .5dt) 
    @. Acache = .5dt * (Acache + W.dW * √(2αkbT))

    @inbounds for i in 1:N
        uloc = @view u[i,:]
        Ax, Ay, Az = @view Acache[i,:]
        a1, a2, a3 = combined_cross(uloc, Ax, Ay, Az)
        detMx[i] =  fdet(a1, a2, a3, -Az, 1, Ax, Ay, -Ax, 1)
        detMy[i] = fdet(1, Az, -Ay, a1, a2, a3, Ay, -Ax, 1)
        detMz[i] = fdet(1, Az, -Ay, -Az, 1, Ax, a1, a2, a3)
        detM[i] = 1 / fdet(1, Az, -Ay, -Az, 1, Ax, Ay, -Ax, 1)
    end
    @. u[:,1] = detMx * detM
    @. u[:,2] = detMy * detM
    @. u[:,3] = detMz * detM
end