import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DrWatson
using DifferentialEquations
using Plots
using DataFrames
using JLD2

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

ω0 = 1.5    # собственная частота
γ  = 0.3    # коэффициент затухания
F0 = 1.0    # амплитуда внешней силы
Ω  = 2.0    # частота внешней силы

x0 = 1.0
y0 = 0.0
tspan = (0.0, 40.0)

function oscillator_free!(du, u, p, t)
    ω0 = p
    x, y = u
    du[1] = y
    du[2] = -ω0^2 * x
end

prob1 = ODEProblem(oscillator_free!, [x0, y0], tspan, ω0)
sol1 = solve(prob1, Tsit5(), saveat=0.05)

function oscillator_damped!(du, u, p, t)
    ω0, γ = p
    x, y = u
    du[1] = y
    du[2] = -ω0^2 * x - 2γ * y
end

prob2 = ODEProblem(oscillator_damped!, [x0, y0], tspan, (ω0, γ))
sol2 = solve(prob2, Tsit5(), saveat=0.05)

function oscillator_forced!(du, u, p, t)
    ω0, γ, F0, Ω = p
    x, y = u
    du[1] = y
    du[2] = -ω0^2 * x - 2γ * y + F0 * sin(Ω * t)
end

prob3 = ODEProblem(oscillator_forced!, [x0, y0], tspan, (ω0, γ, F0, Ω))
sol3 = solve(prob3, Tsit5(), saveat=0.05)

pt = plot(layout=(3,1), size=(900, 900),
    xlabel="t", ylabel="x(t)",
    legend=:topright)

plot!(pt[1], sol1.t, [u[1] for u in sol1.u], label="x(t)", lw=2, color=:blue)
plot!(pt[1], title="Случай 1: без затухания, без силы")

plot!(pt[2], sol2.t, [u[1] for u in sol2.u], label="x(t)", lw=2, color=:red)
plot!(pt[2], title="Случай 2: с затуханием")

plot!(pt[3], sol3.t, [u[1] for u in sol3.u], label="x(t)", lw=2, color=:green)
plot!(pt[3], title="Случай 3: с затуханием и внешней силой")

savefig(plotsdir(script_name, "three_cases_time.png"))

pf = plot(layout=(1,3), size=(1200, 400),
    xlabel="x", ylabel="ẋ", legend=false)

plot!(pf[1], [u[1] for u in sol1.u], [u[2] for u in sol1.u], lw=2, color=:blue)
plot!(pf[1], title="Случай 1")

plot!(pf[2], [u[1] for u in sol2.u], [u[2] for u in sol2.u], lw=2, color=:red)
plot!(pf[2], title="Случай 2")

plot!(pf[3], [u[1] for u in sol3.u], [u[2] for u in sol3.u], lw=2, color=:green)
plot!(pf[3], title="Случай 3")

savefig(plotsdir(script_name, "three_cases_phase.png"))

df = DataFrame(
    случай = ["1: без затухания", "2: с затуханием", "3: с силой"],
    x_final = [last(sol1.u)[1], last(sol2.u)[1], last(sol3.u)[1]],
    y_final = [last(sol1.u)[2], last(sol2.u)[2], last(sol3.u)[2]]
)
println("Итоги:")
println(df)

@save datadir(script_name, "three_cases.jld2") df sol1 sol2 sol3
