struct OrdinaryLLGProblem{F,uType,tType,P,K} <: AbstractODEProblem{uType,tType,true}
    f::F
    u0::uType
    tspan::tType
    p::P
    kwargs::K


    #     function OrdinaryLLGProblem(f::AbstractODEFunction,
    #         u0, tspan, p=NullParameters();
    #         kwargs...)
    #         _tspan = promote_tspan(tspan)
    #         new{typeof(f),
    #             typeof(u0),
    #             typeof(_tspan),
    #             typeof(p),
    #             typeof(kwargs)
    #         }(f, u0, _tspan, p, kwargs)
    #     end
end
@add_kwonly function OrdinaryLLGProblem{true}(f::AbstractODEFunction,
    u0, tspan, p=NullParameters();
    kwargs...)
    _tspan = promote_tspan(tspan)
    OrdinaryLLGProblem{typeof(f),
        typeof(u0),
        typeof(_tspan),
        typeof(p),
        typeof(kwargs)
    }(f, u0, _tspan, p, kwargs)
end

struct StochasticLLGProblem{uType,tType,P,NP,K,ND} <: AbstractSDEProblem{uType,tType,true,ND}
    f::Function
    g::Function
    αkbT::Number
    u0::uType
    tspan::tType
    p::P
    noise_rate_prototype::ND
    noise::NP
    kwargs::K
    seed::UInt64
end

# these two functions create the appropriate LLGproblem, either an OrdinaryLLGProblem or a StochasticLLGProblem
# function LLGProblem(A, u0, tspan, p = NullParameters() ; 
#     kwargs...)
#     _tspan = promote_tspan(tspan)
#     OrdinaryLLGProblem{typeof(u0),typeof(_tspan),typeof(p),typeof(noise),typeof(kwargs),typeof(nothing)}(A, u0, tspan, p, kwargs)
# end

function LLGProblem(A, u0, tspan, p=NullParameters(); αkbT=0.0,
    noise_rate_prototype=nothing,
    seed=UInt64(0),
    kwargs...)
    _tspan = promote_tspan(tspan)

    if αkbT != 0.0
        g(u, p, t) = 1
        f = SciMLBase.SDEFunction{true}(A, g)

        # f = convert(SciMLBase.SDEFunction{true}, A, g)
        noise = StochasticDiffEq.WienerProcess(0.0, zeros(size(u0)), 0.0, save_everystep=false)
        StochasticLLGProblem{typeof(u0),typeof(_tspan),typeof(p),typeof(noise),typeof(kwargs),typeof(nothing)}(f, g, αkbT, u0, tspan, p, noise_rate_prototype, noise, kwargs, seed)
    else
        f = SciMLBase.ODEFunction{true,SciMLBase.DEFAULT_SPECIALIZATION}(A)

        # f = convert(SciMLBase.ODEFunction{true}, A)
        OrdinaryLLGProblem{typeof(f),typeof(u0),typeof(_tspan),typeof(p),typeof(kwargs)}(f, u0, tspan, p, kwargs)
    end

end
# SciMLBase.remaker_of(::OrdinaryLLGProblem) = OrdinaryLLGProblem
abstract type OrdinaryLLGAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type StochasticLLGAlgorithm <: StochasticDiffEqAlgorithm end

StochasticDiffEq.alg_compatible(prob::StochasticLLGProblem, alg::StochasticLLGAlgorithm) = true




function DiffEqBase.__solve(prob::StochasticLLGProblem,
    alg::StochasticLLGAlgorithm,
    timeseries=[], ts=[], ks=nothing, # needed for variable rate
    recompile::Type{Val{recompile_flag}}=Val{true};
    kwargs...) where {recompile_flag}
    integrator = DiffEqBase.__init(prob, alg, timeseries, ts, recompile; kwargs...)
    solve!(integrator)
    integrator.sol
end


function DiffEqBase.__solve(prob::OrdinaryLLGProblem,
    alg::OrdinaryLLGAlgorithm,
    timeseries=(), ts=(), ks=(), # needed for variable rate
    recompile::Type{Val{recompile_flag}}=Val{true};
    kwargs...) where {recompile_flag}
    integrator = DiffEqBase.__init(prob, alg, timeseries, ts, ks, recompile; kwargs...)
    solve!(integrator)
    integrator.sol
end

function llg_f(A!, Acache)
    # remake the function to solve the LLG. A is defined as:
    # ∂m = m × A
    # returned is an iip function f(m) = m × A
    # note that we need to introduce an Acache to prevent allocation
    function f(du, u, p, t)
        A!(Acache, u, p, t)
        @. du[:, 1] = u[:, 2] * Acache[:, 3] - u[:, 3] * Acache[:, 2]
        @. du[:, 2] = u[:, 3] * Acache[:, 1] - u[:, 1] * Acache[:, 3]
        @. du[:, 3] = u[:, 1] * Acache[:, 2] - u[:, 2] * Acache[:, 1]
    end
    f
end

function normalize_llg(norm_cache)
    integrator -> (
        @. norm_cache = sqrt(integrator.u[:, 1]^2 + integrator.u[:, 2]^2 + integrator.u[:, 3]^2);
        @. integrator.u[:, 1] /= norm_cache;
        @. integrator.u[:, 2] /= norm_cache;
        @. integrator.u[:, 3] /= norm_cache)
end

function DiffEqBase.__solve(prob::OrdinaryLLGProblem,
    alg::OrdinaryDiffEqAlgorithm,
    timeseries=(), ts=(), ks=(), # needed for variable rate
    recompile::Type{Val{recompile_flag}}=Val{true};
    kwargs...) where {recompile_flag}
    Acache = similar(prob.u0)

    fA = ODEFunction(llg_f(prob.f, Acache), analytic=prob.f.analytic) # make sure to copy over the analytic part
    prob2 = ODEProblem(fA, prob.u0, prob.tspan, prob.p; prob.kwargs)

    norm_cache = zeros(Float64, size(prob.u0)[1])
    affect! = normalize_llg(norm_cache)
    norm_times = prob.tspan[1]:kwargs[:normalize_dt]:prob.tspan[2]
    if norm_times[end] < prob.tspan[2]
        # the last step will not be normalized, so add it
        append!(norm_times, prob.tspan[2])
    end
    cb = PresetTimeCallback(norm_times, affect!)
    if haskey(kwargs, :callback)
        cb = CallBackSet(kwargs[:callback], cb)
    end
    integrator = DiffEqBase.__init(prob2, alg, timeseries, ts, ks, recompile; callback=cb, kwargs...) # callback is handled by __init, so this should not give problems!


    solve!(integrator)
    integrator.sol
end

