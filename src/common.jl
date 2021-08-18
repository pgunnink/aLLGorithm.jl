struct OrdinaryLLGProblem{uType,tType,P,K} <: AbstractODEProblem{uType,tType,true}
    f::Function
    u0::uType
    tspan::tType
    p::P
    kwargs::K
end

struct StochasticLLGProblem{uType,tType,P,NP,K,ND} <: AbstractSDEProblem{uType,tType,true,ND} 
    f::Function
    αkbT::Number
    u0::uType
    tspan::tType
    p::P
    noise_rate_prototype::ND
    noise::NP
    kwargs::K
end

# these two functions create the appropriate LLGproblem, either an OrdinaryLLGProblem or a StochasticLLGProblem
# function LLGProblem(A, u0, tspan, p = NullParameters() ; 
#     kwargs...)
#     _tspan = promote_tspan(tspan)
#     OrdinaryLLGProblem{typeof(u0),typeof(_tspan),typeof(p),typeof(noise),typeof(kwargs),typeof(nothing)}(A, u0, tspan, p, kwargs)
# end

function LLGProblem(A, u0, tspan, p = NullParameters() ; αkbT = 0.,
        noise_rate_prototype = nothing,
        noise = nothing,
        kwargs...)
    _tspan = promote_tspan(tspan)

    if αkbT != 0.
        f = convert(SciMLBase.SDEFunction{true}, A,  x -> 1)

        StochasticLLGProblem{typeof(u0),typeof(_tspan),typeof(p),typeof(noise),typeof(kwargs),typeof(nothing)}(f, αkbT, u0, tspan, p, noise_rate_prototype, noise, kwargs)
    else
        f = convert(SciMLBase.ODEFunction{true}, A)

        OrdinaryLLGProblem{typeof(u0),typeof(_tspan),typeof(p),typeof(kwargs)}(f, u0, tspan, p, kwargs)
    end

end

abstract type OrdinaryLLGAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type StochasticLLGAlgorithm <: StochasticDiffEqAlgorithm end

alg_compatible(prob::OrdinaryLLGProblem,alg::OrdinaryLLGAlgorithm) = true
alg_compatible(prob::StochasticLLGProblem,alg::StochasticLLGAlgorithm) = true




function DiffEqBase.__solve(prob::StochasticLLGProblem,
    alg::StochasticLLGAlgorithm,
    timeseries = (),ts = (),ks = (), # needed for variable rate
recompile::Type{Val{recompile_flag}} = Val{true};
    kwargs...) where recompile_flag
    integrator = DiffEqBase.__init(prob, alg, timeseries, ts, recompile;kwargs...)
    solve!(integrator)
    integrator.sol
end


function DiffEqBase.__solve(prob::OrdinaryLLGProblem,
    alg::OrdinaryLLGAlgorithm,
    timeseries = (),ts = (),ks = (), # needed for variable rate
recompile::Type{Val{recompile_flag}} = Val{true};
    kwargs...) where recompile_flag
    integrator = DiffEqBase.__init(prob, alg, timeseries, ts, ks, recompile;kwargs...)
    solve!(integrator)
    integrator.sol
end

# function llg_f(A, N)
#     # remake the function to solve the LLG. A is defined as:
#     # ∂m = m × A
#     # returned is an iip function f(m) = m × A
#     function f(du, u, p, t)

#     end
#     f
# end

# function DiffEqBase.__solve(prob::OrdinaryLLGProblem,
#     alg::OrdinaryDiffEqAlgorithm,
#     timeseries = (),ts = (),ks = (), # needed for variable rate
# recompile::Type{Val{recompile_flag}} = Val{true};
#     kwargs...) where recompile_flag
#     prob2 = remake(prob; f = llg_f(prob.f))

#     integrator = DiffEqBase.__init(prob, alg, timeseries, ts, ks, recompile;kwargs...)
#     solve!(integrator)
#     integrator.sol
# end

