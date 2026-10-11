# # Модель эпидемии (SIR)
#
# **Задание 1 (методичка).** Исследовать простейшую модель эпидемии:
#
# $$\frac{dS}{dt}, \frac{dI}{dt}, \frac{dR}{dt}$$
#
# с двумя случаями: $I(0) \le I^*$ и $I(0) > I^*$.

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

# ## Параметры

N   = 2000.0    # общая численность популяции
α   = 0.01      # коэффициент заболеваемости
β   = 0.02      # коэффициент выздоровления
I_star = 150.0  # критическое значение числа инфицированных

I0 = 300.0
R0 = 0.0
S0 = N - I0 - R0
tspan = (0.0, 200.0)

println("Начальные условия:")
println("  N  = $N")
println("  S0 = $S0, I0 = $I0, R0 = $R0")
println("  α  = $α, β = $β, I* = $I_star")

# ## Случай A: I(0) ≤ I*
#
# Если число инфицированных ниже критического, распространения нет:
# восприимчивые не заражаются. Инфицированные только выздоравливают.

function sir_below!(du, u, p, t)
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

# Случай A: искусственно делаем I* больше I0
probA = ODEProblem(sir_below!, [S0, I0, R0], tspan, (α, β, I0 + 100.0))
solA = solve(probA, Tsit5(), saveat=0.1)

# ## Случай B: I(0) > I*
#
# Эпидемия развивается: восприимчивые заражаются, потом выздоравливают.

function sir_above!(du, u, p, t)
    α, β, I_star = p
    S, I, R = u
    du[1] = -α * S
    du[2] = α * S - β * I
    du[3] = β * I
end

probB = ODEProblem(sir_above!, [S0, I0, R0], tspan, (α, β, I0 - 100.0))
solB = solve(probB, Tsit5(), saveat=0.1)

# ## Визуализация

p = plot(layout=(2,1), size=(900, 800),
    xlabel="t", ylabel="Численность",
    legend=:topright)

plot!(p[1], solA.t, [u[1] for u in solA.u], label="S(t)", lw=2, color=:blue)
plot!(p[1], solA.t, [u[2] for u in solA.u], label="I(t)", lw=2, color=:red)
plot!(p[1], solA.t, [u[3] for u in solA.u], label="R(t)", lw=2, color=:green)
plot!(p[1], title="Случай A: I(0) ≤ I*")

plot!(p[2], solB.t, [u[1] for u in solB.u], label="S(t)", lw=2, color=:blue)
plot!(p[2], solB.t, [u[2] for u in solB.u], label="I(t)", lw=2, color=:red)
plot!(p[2], solB.t, [u[3] for u in solB.u], label="R(t)", lw=2, color=:green)
plot!(p[2], title="Случай B: I(0) > I*")

savefig(plotsdir(script_name, "epidemic_both_cases.png"))

# ## Сводка

df = DataFrame(
    случай = ["A (I₀ ≤ I*)", "B (I₀ > I*)"],
    S_final = [last(solA.u)[1], last(solB.u)[1]],
    I_final = [last(solA.u)[2], last(solB.u)[2]],
    R_final = [last(solA.u)[3], last(solB.u)[3]]
)
println("\nИтоги:")
println(df)

@save datadir(script_name, "epidemic.jld2") df solA solB
