using aLLGorithm
using Test
using Symbolics
include("llg_problems.jl")


@variables x, y, z, t


# @inferred aLLGorithm.OrdinaryLLGProblem single_spin_problem()
function runtimefunction()
    γ = 0.2
    α = 0.01

    expr = [-y * γ * α, x * γ * α, -γ]
    Symbolics.build_function(expr, [x y z], t, expression=Val{false})[2]
end

# @inferred aLLGorithm.OrdinaryLLGProblem LLGProblem(runtimefunction(), [1.0 2.0 3.0], (0.0, 10.0))

@testset "Solve" begin
    @inferred SciMLBase.ODESolution solve(LLGProblem(runtimefunction(), [1.0 2.0 3.0], (0.0, 10.0)), SIBOrdinary(), dt=0.1)
end


@testset "ODE solve" begin
    solve(LLGProblem(runtimefunction(), [1.0 2.0 3.0], (0.0, 10.0), αkbT=0.001), SIBStochastic(), dt=0.1)
end


