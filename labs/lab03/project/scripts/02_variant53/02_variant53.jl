using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots
using DataFrames
using JLD2

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

x0 = 321000.0
y0 = 123000.0
tspan = (0.0, 30.0)

function model1!(du, u, p, t)
    x, y = u
    du[1] = -0.336*x - 0.877*y + sin(t + 1) + 1
    du[2] = -0.441*x - 0.232*y + cos(t + 2) + 1
end

prob1 = ODEProblem(model1!, [x0, y0], tspan)
sol1 = solve(prob1, Tsit5(), saveat=0.05)

function model2!(du, u, p, t)
    x, y = u
    du[1] = -0.432*x - 0.815*y + sin(2t) + 2
    du[2] = -0.336*x*y - 0.245*y + cos(t) + 2
end

prob2 = ODEProblem(model2!, [x0, y0], tspan)
sol2 = solve(prob2, Tsit5(), saveat=0.05)

p = plot(layout=(2,1), size=(900, 700),
    xlabel="t (дни)", ylabel="Численность",
    legend=:topright)

plot!(p[1], sol1.t, [u[1] for u in sol1.u], label="X (регулярные)", lw=2, color=:blue)
plot!(p[1], sol1.t, [u[2] for u in sol1.u], label="Y (регулярные)", lw=2, color=:red)
plot!(p[1], title="Модель 1: регулярные vs регулярные (вариант 53)")

plot!(p[2], sol2.t, [u[1] for u in sol2.u], label="X (регулярные)", lw=2, color=:blue)
plot!(p[2], sol2.t, [u[2] for u in sol2.u], label="Y (партизаны)", lw=2, color=:red)
plot!(p[2], title="Модель 2: регулярные vs партизаны (вариант 53)")

savefig(plotsdir(script_name, "variant53_models.png"))

function winner(sol)
    xf = last(sol.u)[1]
    yf = last(sol.u)[2]
    if xf > yf
        return "X (x=$xf)"
    elseif yf > xf
        return "Y (y=$yf)"
    else
        return "Ничья"
    end
end

println("Модель 1 (регулярные): ", winner(sol1))
println("Модель 2 (партизаны):  ", winner(sol2))

println("\nАнализ условий победы:")
println("  Модель 1: регулярные vs регулярные")
println("    X побеждает, если b·y₀² > c·x₀²")
println("    b·y₀² = ", 0.877 * y0^2)
println("    c·x₀² = ", 0.441 * x0^2)
println("    Победитель (теория): ", 0.877*y0^2 > 0.441*x0^2 ? "X" : "Y")

println("\n  Модель 2: регулярные vs партизаны")
println("    X побеждает, если (b/2)·x₀² > c·y₀")
println("    (b/2)·x₀² = ", (0.815/2) * x0^2)
println("    c·y₀ = ", 0.336 * y0)
println("    Победитель (теория): ", (0.815/2)*x0^2 > 0.336*y0 ? "X" : "Y")

df = DataFrame(
    модель = ["Регулярные vs регулярные", "Регулярные vs партизаны"],
    x_финал = [last(sol1.u)[1], last(sol2.u)[1]],
    y_финал = [last(sol1.u)[2], last(sol2.u)[2]],
    победитель = [winner(sol1), winner(sol2)]
)
println("\nИтоги:")
println(df)

@save datadir(script_name, "variant53.jld2") df sol1 sol2
