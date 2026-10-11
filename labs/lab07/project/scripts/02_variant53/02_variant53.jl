import Pkg
Pkg.activate("/home/dabakharev/work/study/2026-1/2026-1==study--mathmod/labs/lab07/project")
using DrWatson
using DifferentialEquations
using Plots
Plots.default(fmt = :png)
using DataFrames
using JLD2

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

N  = 1304.0
n0 = 10.0
tspan = (0.0, 20.0)

println("Параметры варианта 53:")
println("  N = $N, n0 = $n0")

function model1!(du, u, p, t)
    N = p
    n = u[1]
    du[1] = (0.76 + 0.000016 * n) * (N - n)
end

prob1 = ODEProblem(model1!, [n0], tspan, N)
sol1 = solve(prob1, Tsit5(), saveat=0.05)

function model2!(du, u, p, t)
    N = p
    n = u[1]
    du[1] = (0.000016 + 0.6 * n) * (N - n)
end

prob2 = ODEProblem(model2!, [n0], tspan, N)
sol2 = solve(prob2, Tsit5(), saveat=0.001, abstol=1e-8, reltol=1e-8)

function model3!(du, u, p, t)
    N = p
    n = u[1]
    α1 = 0.7 * abs(sin(7t))
    α2 = 0.7 * abs(sin(3t))
    du[1] = (α1 + α2 * n) * (N - n)
end

prob3 = ODEProblem(model3!, [n0], tspan, N)
sol3 = solve(prob3, Tsit5(), saveat=0.05)

N_half = N / 2
idx = findfirst(u -> u[1] >= N_half, sol2.u)

if idx !== nothing && idx > 1
    n_prev = sol2.u[idx-1][1]
    n_curr = sol2.u[idx][1]
    t_prev = sol2.t[idx-1]
    t_curr = sol2.t[idx]
    t_max = t_prev + (N_half - n_prev) / (n_curr - n_prev) * (t_curr - t_prev)
else
    t_max = NaN
end

println("\nСлучай 2: n(t) достигает N/2 = $N_half")
println("  при t ≈ $t_max")

p = plot(layout=(3,1), size=(900, 900),
    xlabel="t", ylabel="n(t)",
    legend=:bottomright)

plot!(p[1], sol1.t, [u[1] for u in sol1.u], label="n(t)", lw=2, color=:blue)
plot!(p[1], title="Случай 1: преобладает платная реклама")

plot!(p[2], sol2.t, [u[1] for u in sol2.u], label="n(t)", lw=2, color=:red)
vline!(p[2], [t_max], label="t_max ≈ $t_max", color=:black, linestyle=:dash)
plot!(p[2], title="Случай 2: преобладает сарафанное радио")

plot!(p[3], sol3.t, [u[1] for u in sol3.u], label="n(t)", lw=2, color=:green)
plot!(p[3], title="Случай 3: коэффициенты зависят от времени")

savefig(plotsdir(script_name, "variant53.png"))

v2 = [(0.000016 + 0.6 * sol2.u[i][1]) * (N - sol2.u[i][1]) for i in 1:length(sol2.u)]

pv = plot(sol2.t, v2, lw=2, color=:red,
    xlabel="t", ylabel="dn/dt",
    title="Скорость распространения (случай 2)",
    legend=false)
vline!(pv, [t_max], color=:black, linestyle=:dash, label="t_max")
savefig(plotsdir(script_name, "variant53_speed.png"))

df = DataFrame(
    случай = ["1: платная", "2: сарафан", "3: с sin"],
    n_final = [last(sol1.u)[1], last(sol2.u)[1], last(sol3.u)[1]]
)
println(df)

@save datadir(script_name, "variant53.jld2") df sol1 sol2 sol3 t_max
