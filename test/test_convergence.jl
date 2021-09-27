using aLLGorithm, Test, DiffEqDevTools
import DifferentialEquations: Euler

include("llg_problems.jl")

dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts, single_spin_problem(), Euler(), normalize_dt = dts[end])

sim.𝒪est[:final]