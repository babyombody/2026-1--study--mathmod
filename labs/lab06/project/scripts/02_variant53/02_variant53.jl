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

I_star_A = 200.0
probA = ODEProblem(sir_with_threshold!, [S0, I0, R0], tspan, (α, β, I_star_A))
solA = solve(probA, Tsit5(), saveat=0.1)

I_star_B = 150.0
probB = ODEProblem(sir_with_threshold!, [S0, I0, R0], tspan, (α, β, I_star_B))
solB = solve(probB, Tsit5(), saveat=0.1)

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

df = DataFrame(
    случай = ["A (I₀ ≤ I*)", "B (I₀ > I*)"],
    S_final = [last(solA.u)[1], last(solB.u)[1]],
    I_final = [last(solA.u)[2], last(solB.u)[2]],
    R_final = [last(solA.u)[3], last(solB.u)[3]]
)
println("\nИтоги варианта 53:")
println(df)

@save datadir(script_name, "variant53.jld2") df solA solB
