# # Модель эпидемии: вариант 53
#
# **Условие:**
# - N = 6 159
# - I(0) = 173
# - R(0) = 61
# - S(0) = N − I(0) − R(0) = 5 925
#
# Рассмотреть два случая: I(0) ≤ I* и I(0) > I*.
# α и β задаём самостоятельно: α = 0.01, β = 0.02.
# I* = 150, тогда I(0) = 173 > I* (случай Б).
# Для случая А возьмём I* = 200.

import Pkg
Pkg.activate("/home/dabakharev/work/study/2026-1/2026-1==study--mathmod/labs/lab06/project")
using DrWatson
using DifferentialEquations
using Plots
Plots.default(fmt = :png)
using DataFrames
using JLD2

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# ## Параметры варианта 53

N  = 6159.0
I0 = 173.0
R0 = 61.0
S0 = N - I0 - R0

α = 0.01
β = 0.02
tspan = (0.0, 300.0)

println("Параметры варианта 53:")
println("  N  = $N")
println("  S0 = $S0, I0 = $I0, R0 = $R0")
println("  α  = $α, β = $β")

# ## Функция правой части SIR с порогом I*

function sir_with_threshold!(du, u, p, t)
    α, β, I_star = p
    S, I, R = u
    if I > I_star
        du[1] = -α * S
        du[2] = α * S - β * I
    else
        du[1] = 0.0
        du[2] = -β * I
    end
    du[3] = β * I
end

# ## Случай A: I(0) ≤ I*
#
# Задаём I* = 200 (больше I(0)), эпидемия не развивается.

I_star_A = 200.0
probA = ODEProblem(sir_with_threshold!, [S0, I0, R0], tspan, (α, β, I_star_A))
solA = solve(probA, Tsit5(), saveat=0.1)

# ## Случай B: I(0) > I*
#
# Задаём I* = 150 (меньше I(0)), эпидемия развивается.

I_star_B = 150.0
probB = ODEProblem(sir_with_threshold!, [S0, I0, R0], tspan, (α, β, I_star_B))
solB = solve(probB, Tsit5(), saveat=0.1)

# ## Визуализация

p = plot(layout=(2,1), size=(900, 800),
    xlabel="t", ylabel="Численность",
    legend=:topright)

plot!(p[1], solA.t, [u[1] for u in solA.u], label="S(t)", lw=2, color=:blue)
plot!(p[1], solA.t, [u[2] for u in solA.u], label="I(t)", lw=2, color=:red)
plot!(p[1], solA.t, [u[3] for u in solA.u], label="R(t)", lw=2, color=:green)
plot!(p[1], title="Вариант 53, случай A: I(0)=$I0 ≤ I*=$I_star_A")

plot!(p[2], solB.t, [u[1] for u in solB.u], label="S(t)", lw=2, color=:blue)
plot!(p[2], solB.t, [u[2] for u in solB.u], label="I(t)", lw=2, color=:red)
plot!(p[2], solB.t, [u[3] for u in solB.u], label="R(t)", lw=2, color=:green)
plot!(p[2], title="Вариант 53, случай B: I(0)=$I0 > I*=$I_star_B")

savefig(plotsdir(script_name, "variant53.png"))

# ## Сводка

df = DataFrame(
    случай = ["A (I₀ ≤ I*)", "B (I₀ > I*)"],
    S_final = [last(solA.u)[1], last(solB.u)[1]],
    I_final = [last(solA.u)[2], last(solB.u)[2]],
    R_final = [last(solA.u)[3], last(solB.u)[3]]
)
println("\nИтоги варианта 53:")
println(df)

@save datadir(script_name, "variant53.jld2") df solA solB
