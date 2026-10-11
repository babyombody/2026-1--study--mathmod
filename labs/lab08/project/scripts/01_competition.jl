# # Модель конкуренции двух фирм
#
# **Задание 1 (методичка).** Исследовать две модели конкуренции:
#
# 1. Чисто рыночная конкуренция (симметричный член M₁M₂);
# 2. Конкуренция с социально-психологическим фактором (асимметричный M₁M₂).

import Pkg
Pkg.activate("/home/dabakharev/work/study/2026-1/2026-1==study--mathmod/labs/lab08/project")
using DrWatson
using DifferentialEquations
using Plots
Plots.default(fmt = :png)
using DataFrames
using JLD2

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# ## Параметры (из примера 1 методички)

p_cr = 20.0
N    = 10.0
q    = 1.0
τ1   = 10.0
τ2   = 16.0
p1   = 9.0
p2   = 7.0

# ## Коэффициенты системы

a1 = p_cr / (τ1^2 * p1^2 * N * q)
a2 = p_cr / (τ2^2 * p2^2 * N * q)
b  = p_cr / (τ1^2 * p1^2 * τ2^2 * p2^2 * N * q)
c1 = (p_cr - p1) / (τ1 * p1)
c2 = (p_cr - p2) / (τ2 * p2)

println("Коэффициенты:")
println("  a₁ = $a1, a₂ = $a2")
println("  b  = $b")
println("  c₁ = $c1, c₂ = $c2")

# ## Начальные условия

M0 = [2.0, 1.0]
θ_span = (0.0, 30.0)

# ## Случай 1: симметричная конкуренция

function competition_sym!(du, u, p, t)
    a1, a2, b, c1, c2 = p
    M1, M2 = u
    du[1] = M1 - (b/c1) * M1 * M2 - (a1/c1) * M1^2
    du[2] = (c2/c1) * M2 - (b/c1) * M1 * M2 - (a2/c1) * M2^2
end

prob1 = ODEProblem(competition_sym!, M0, θ_span, (a1, a2, b, c1, c2))
sol1 = solve(prob1, Tsit5(), saveat=0.05)

# ## Случай 2: асимметричная конкуренция (соц-психологический фактор)

δ = 0.002

function competition_asym!(du, u, p, t)
    a1, a2, b, c1, c2, δ = p
    M1, M2 = u
    du[1] = M1 - (b/c1 + δ) * M1 * M2 - (a1/c1) * M1^2
    du[2] = (c2/c1) * M2 - (b/c1) * M1 * M2 - (a2/c1) * M2^2
end

prob2 = ODEProblem(competition_asym!, M0, θ_span, (a1, a2, b, c1, c2, δ))
sol2 = solve(prob2, Tsit5(), saveat=0.05)

# ## Визуализация

p = plot(layout=(2,1), size=(900, 800),
    xlabel="θ (безразмерное время)", ylabel="M",
    legend=:topright)

plot!(p[1], sol1.t, [u[1] for u in sol1.u], label="M₁", lw=2, color=:blue)
plot!(p[1], sol1.t, [u[2] for u in sol1.u], label="M₂", lw=2, color=:green)
plot!(p[1], title="Случай 1: чисто рыночная конкуренция")

plot!(p[2], sol2.t, [u[1] for u in sol2.u], label="M₁", lw=2, color=:blue)
plot!(p[2], sol2.t, [u[2] for u in sol2.u], label="M₂", lw=2, color=:green)
plot!(p[2], title="Случай 2: с соц-психологическим фактором")

savefig(plotsdir(script_name, "competition.png"))

# ## Итог

df = DataFrame(
    случай = ["1: симметричный", "2: асимметричный"],
    M1_final = [last(sol1.u)[1], last(sol2.u)[1]],
    M2_final = [last(sol1.u)[2], last(sol2.u)[2]]
)
println(df)

@save datadir(script_name, "competition.jld2") df sol1 sol2
