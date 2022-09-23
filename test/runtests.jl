using aLLGorithm
using Test
using Symbolics
include("llg_problems.jl")




@inferred aLLGorithm.OrdinaryLLGProblem single_spin_problem()
function runtimefunction()
    @variables x, y, z, t
    γ = 0.2
    α = 0.01

    expr = [-y * γ * α, x * γ * α, -γ]
    Symbolics.build_function(expr, [], t, expression=Val{false})[2]
end

@inferred aLLGorithm.OrdinaryLLGProblem LLGProblem(runtimefunction(), [1.0 2.0 3.0], (0.0, 10.0))
