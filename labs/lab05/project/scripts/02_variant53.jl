# # Модель «хищник-жертва»: вариант 53
#
# **Условие:**
# $$\frac{dx}{dt} = -0.34x(t) + 0.051x(t)y(t)$$
# $$\frac{dy}{dt} = 0.33y(t) - 0.031x(t)y(t)$$
#
# Начальные условия: **x₀ = 9, y₀ = 30**.

import Pkg
Pkg.activate("/home/dabakharev/work/study/2026-1/2026-1==study--mathmod/labs/lab05/project")
using DrWatson
using DifferentialEquations
using Plots
using DataFrames
using JLD2

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# ## Параметры варианта 53

x0 = 9.0
y0 = 30.0
tspan = (0.0, 200.0)

# Система из условия:
#   dx/dt = -0.34x + 0.051xy
#   dy/dt =  0.33y - 0.031xy
#
# В «жертвы-хищники» x — это хищники, y — жертвы.
# Но формально система та же.

function predator_prey!(du, u, p, t)
    x, y = u
    du[1] = -0.34*x + 0.051*x*y
    du[2] =  0.33*y - 0.031*x*y
end

prob = ODEProblem(predator_prey!, [x0, y0], tspan)
sol = solve(prob, Tsit5(), saveat=0.1)

# ## Стационарное состояние
#
# Из dx/dt = 0: x(-0.34 + 0.051y) = 0 → y* = 0.34 / 0.051
# Из dy/dt = 0: y(0.33 - 0.031x) = 0 → x* = 0.33 / 0.031

x_star = 0.33 / 0.031
y_star = 0.34 / 0.051
println("Стационарное состояние (вариант 53):")
println("  x* = $x_star")
println("  y* = $y_star")

# ## Временные ряды x(t), y(t)

pt = plot(layout=(2,1), size=(900, 700),
    xlabel="t", ylabel="Численность",
    legend=:topright)

plot!(pt[1], sol.t, [u[1] for u in sol.u], label="x(t)", lw=2, color=:blue)
plot!(pt[1], sol.t, [u[2] for u in sol.u], label="y(t)", lw=2, color=:red)
plot!(pt[1], title="Вариант 53: динамика популяций")

plot!(pt[2], [u[1] for u in sol.u], [u[2] for u in sol.u], lw=2, color=:purple, label="Фазовая траектория")
scatter!(pt[2], [x_star], [y_star], markersize=8, color=:black, label="Стационарная точка")
plot!(pt[2], xlabel="x", ylabel="y", title="Вариант 53: фазовый портрет")

savefig(plotsdir(script_name, "variant53.png"))

# ## Сводка

df = DataFrame(
    параметр = ["x0", "y0", "x*", "y*"],
    значение = [x0, y0, x_star, y_star]
)
println("\nИтоги варианта 53:")
println(df)

@save datadir(script_name, "variant53.jld2") df sol x_star y_star
