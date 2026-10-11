# # Модель Лотки-Вольтерры: базовый анализ
#
# **Задание 1 (методичка).** Исследовать модель «хищник-жертва»:
#
# $$\frac{dx}{dt} = ax - bxy, \quad \frac{dy}{dt} = -cy + dxy$$
#
# где $x$ — число жертв, $y$ — число хищников.
#
# Коэффициенты задаём самостоятельно (близко к примеру 1 из методички).

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

# ## Параметры модели

a = 0.5   # скорость естественного прироста жертв
b = 0.05  # коэффициент взаимодействия (гибель жертв)
c = 0.2   # коэффициент естественной смертности хищников
d = 0.02  # коэффициент увеличения числа хищников

x0 = 5.0
y0 = 10.0
tspan = (0.0, 400.0)

# ## Система уравнений Лотки-Вольтерры
#
# В методичке обозначения чуть другие (a,b,c,d — смертность/прирост),
# но структура та же: жертвы растут, хищники убывают без жертв.

function lotka_volterra!(du, u, p, t)
    a, b, c, d = p
    x, y = u
    du[1] = a*x - b*x*y
    du[2] = -c*y + d*x*y
end

prob = ODEProblem(lotka_volterra!, [x0, y0], tspan, (a, b, c, d))
sol = solve(prob, Tsit5(), saveat=0.1)

# ## Стационарное состояние
#
# Ненулевое стационарное состояние:
# $$x^* = \frac{c}{d}, \quad y^* = \frac{a}{b}$$

x_star = c / d
y_star = a / b
println("Стационарное состояние:")
println("  x* = $x_star (жертвы)")
println("  y* = $y_star (хищники)")

# ## Временные ряды x(t), y(t)

pt = plot(layout=(2,1), size=(900, 700),
    xlabel="t", ylabel="Численность",
    legend=:topright)

plot!(pt[1], sol.t, [u[1] for u in sol.u], label="x(t) — жертвы", lw=2, color=:blue)
plot!(pt[1], sol.t, [u[2] for u in sol.u], label="y(t) — хищники", lw=2, color=:red)
plot!(pt[1], title="Динамика популяций (Лотка-Вольтерра)")

# Фазовый портрет
plot!(pt[2], [u[1] for u in sol.u], [u[2] for u in sol.u], lw=2, color=:purple, label="Фазовая траектория")
scatter!(pt[2], [x_star], [y_star], markersize=8, color=:black, label="Стационарная точка (x*, y*)")
plot!(pt[2], xlabel="x (жертвы)", ylabel="y (хищники)", title="Фазовый портрет")

savefig(plotsdir(script_name, "lotka_volterra.png"))

# ## Сводка

df = DataFrame(
    параметр = ["a", "b", "c", "d", "x*", "y*"],
    значение = [a, b, c, d, x_star, y_star]
)
println("\nПараметры и стационарное состояние:")
println(df)

@save datadir(script_name, "lotka_volterra.jld2") df sol x_star y_star
