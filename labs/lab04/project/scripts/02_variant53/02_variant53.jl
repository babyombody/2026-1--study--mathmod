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

x0 = 1.9
y0 = 0.9
tspan = (0.0, 49.0)

ω0_sq_1 = 1.9

function osc_free!(du, u, p, t)
    ω0_sq = p
    x, y = u
    du[1] = y
    du[2] = -ω0_sq * x
end

prob1 = ODEProblem(osc_free!, [x0, y0], tspan, ω0_sq_1)
sol1 = solve(prob1, Tsit5(), saveat=0.05)

γ2 = 2.9 / 2
ω0_sq_2 = 3.9

function osc_damped!(du, u, p, t)
    ω0_sq, γ = p
    x, y = u
    du[1] = y
    du[2] = -ω0_sq * x - 2γ * y
end

prob2 = ODEProblem(osc_damped!, [x0, y0], tspan, (ω0_sq_2, γ2))
sol2 = solve(prob2, Tsit5(), saveat=0.05)

γ3 = 4.9 / 2
ω0_sq_3 = 5.9
F0_3 = 6.9
Ω3 = 7.9

function osc_forced!(du, u, p, t)
    ω0_sq, γ, F0, Ω = p
    x, y = u
    du[1] = y
    du[2] = -ω0_sq * x - 2γ * y + F0 * sin(Ω * t)
end

prob3 = ODEProblem(osc_forced!, [x0, y0], tspan, (ω0_sq_3, γ3, F0_3, Ω3))
sol3 = solve(prob3, Tsit5(), saveat=0.05)

pt = plot(layout=(3,1), size=(900, 900),
    xlabel="t", ylabel="x(t)",
    legend=:topright)

plot!(pt[1], sol1.t, [u[1] for u in sol1.u], label="x(t)", lw=2, color=:blue)
plot!(pt[1], title="Вариант 53: ẍ + 1.9x = 0")

plot!(pt[2], sol2.t, [u[1] for u in sol2.u], label="x(t)", lw=2, color=:red)
plot!(pt[2], title="Вариант 53: ẍ + 2.9ẋ + 3.9x = 0")

plot!(pt[3], sol3.t, [u[1] for u in sol3.u], label="x(t)", lw=2, color=:green)
plot!(pt[3], title="Вариант 53: ẍ + 4.9ẋ + 5.9x = 6.9 sin(7.9t)")

savefig(plotsdir(script_name, "variant53_time.png"))

pf = plot(layout=(1,3), size=(1200, 400),
    xlabel="x", ylabel="ẋ", legend=false)

plot!(pf[1], [u[1] for u in sol1.u], [u[2] for u in sol1.u], lw=2, color=:blue)
plot!(pf[1], title="Случай 1")

plot!(pf[2], [u[1] for u in sol2.u], [u[2] for u in sol2.u], lw=2, color=:red)
plot!(pf[2], title="Случай 2")

plot!(pf[3], [u[1] for u in sol3.u], [u[2] for u in sol3.u], lw=2, color=:green)
plot!(pf[3], title="Случай 3")

savefig(plotsdir(script_name, "variant53_phase.png"))

df = DataFrame(
    случай = ["1: без затухания", "2: с затуханием", "3: с силой"],
    x_final = [last(sol1.u)[1], last(sol2.u)[1], last(sol3.u)[1]],
    y_final = [last(sol1.u)[2], last(sol2.u)[2], last(sol3.u)[2]]
)
println("Итоги варианта 53:")
println(df)

@save datadir(script_name, "variant53.jld2") df sol1 sol2 sol3
