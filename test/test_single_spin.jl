using aLLGorithm
using Test

γ = .2
α = .01
function A!(Acache, u, t)
    Acache[1] = - u[2] * γ * α
    Acache[2] = u[1] * γ * α
    Acache[3] = -γ
end

function analytical_solution(u, p, t)
    ϕ = γ * t
    nz = tanh(α * ϕ)
    nx = cos(ϕ) * √(1 - nz^2)
    ny = sin(ϕ) * √(1 - nz^2)
    [nx ny nz]
end

f = ODEFunction(A!, analytic = analytical_solution)

prob =  LLGProblem(f, [1. 0. 0.], (0., 100.))
# integrator = init(prob, SIBOrdinary(), dt = 0.01)
# step!(integrator)
sol = solve(prob, SIBOrdinary(), dt = 0.01)



##

using Plots
plot(sol, denseplot = false, plot_analytic = true)
