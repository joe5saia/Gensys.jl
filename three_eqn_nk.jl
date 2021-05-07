include("src/Gensys.jl")
using Plots

θ = 1
β = 0.97
κ = 0.172
Ψπ = 1.5
Ψy = 0.15
ρy = 0.90

Y = [:y, :π, :r, :Ey, :Eπ, :uis, :uis_grwth]
z = [:uis, :uis_grwth]
η = [:y, :π]

function set_values!(M::Matrix, Y::Vector{Symbol}, eqn::Integer, variable::Symbol, value)
  M[eqn, findfirst(Y .== variable)] = value
end


function set_values!(M::Matrix, Y::Vector{Symbol}, eqn::Integer, variable::Vector{Symbol}, value::Vector)
  for (var, val) in zip(variable, value)
    set_values!(M, Y, eqn, var, val)
  end
end

# Γ₀ Y_t = Γ₁ Y_t-1 + C + Ψ z_t + Π η_t
Γ₀ = zeros(Float64, length(Y), length(Y))
set_values!(Γ₀, Y, 1, [:y, :Ey, :r, :uis], [1, -1, 1/θ, -1])
set_values!(Γ₀, Y, 2, [:π, :Eπ, :y], [1, -β, -κ])
set_values!(Γ₀, Y, 3, [:r, :y, :π], [1, -ϕy, -ϕπ])
set_values!(Γ₀, Y, 4, [:y], [1])
set_values!(Γ₀, Y, 5, [:π], [1])
set_values!(Γ₀, Y, 6, [:uis, :uis_grwth], [1, -0.95])
set_values!(Γ₀, Y, 7, [:uis_grwth], [1])

Γ₁ = zeros(Float64, length(Y), length(Y))
set_values!(Γ₁, Y, 4, :Ey, 1)
set_values!(Γ₁, Y, 5, :Eπ, 1)
set_values!(Γ₁, Y, 6, :uis, ρy)
set_values!(Γ₁, Y, 7, :uis_grwth, 0.75)

Ψ = zeros(Float64, length(Y), length(z))
set_values!(Ψ, z, 6, :uis, 1)
set_values!(Ψ, z, 7, :uis_grwth, 1)

Π = zeros(Float64, length(Y), length(η))
set_values!(Π, η, 4, [:y], [1])
set_values!(Π, η, 5, [:π], [1])

C = zeros(Float64, length(Y), 1)

# Solve
G1, Const, impact, fmat, fwt, ywt, gev, eu, loose = Gensys.gensysdt(Γ₀, Γ₁, C, Ψ, Π);

# Simulation
T = 60
shock_t = 25

function plot_sim(x, Y, vars)
  p = plot(x[1,:], label="Output", legend=:topleft);
  plot!(x[2,:], label="Inflation");
  plot!([0; x[4,1:end-1]], label="Expected Output");
  plot!(x[end-1,:], label="Persistient IS_shock");
  plot!(x[end,:], label="Transient IS_shock");
  vline!([shock_t], color="black", linestyle=:dot, label="Shock_date");
  vline!([findmax(abs.(x[1,:]))[2] + findfirst(abs.(x[1,findmax(abs.(x[1,:]))[2]:end]) .< 0.1)], color="black", label="Normal")
  return p
end

function plot_sim(x, Y, vars)
  p = plot(bg=:black, legend=:topleft);
  for v in vars
    y = x[findfirst(Y .== v), :]
    plot!(p, y, label="$(v)")
  end
  vline!([shock_t], color="white", linestyle=:dot, label="Shock_date");
  return p
end

function shock_sim(u, Y, G1, impact)
  T = size(u,2)
  x = zeros(length(Y), T)
  for t in 2:T
    x[:, t] = G1*x[:, t-1] + impact*u[:, t]
  end
  return x
end

function stock_price(y, i, γy=1.5, γi=1, β=β)
  p = 0
  for t in 1:length(y)
    p += β^t*(γy * y[t] - γi * i[t])
  end
  return p
end


# Level shock
u = zeros(size(Ψ,2), T+1);
u[1, shock_t] = 1;
x = shock_sim(u, Y, G1, impact)
pss = plot_sim(x, Y, [:y, :r, :uis, :uis_grwth])
title!(pss, "Response at Steady State");


u = zeros(size(Ψ,2), T+1);
u[1, 1:shock_t-2] .= -0.15;
x0 = shock_sim(u, Y, G1, impact)
pr0 = plot_sim(x0, Y, [:y, :r, :uis, :uis_grwth])
title!(pr0, "Base Recovery from Recession");


u = zeros(size(Ψ,2), T+1);
u[1, 1:shock_t-2] .= -0.15;
u[1, shock_t] = 1;
x1 = shock_sim(u, Y, G1, impact)
pr1 = plot_sim(x1, Y, [:y, :r, :uis, :uis_grwth])
title!(pr1, "Response in Recession");


prdiff = plot_sim(x1 .- x0, Y, [:y, :r, :uis, :uis_grwth])
title!(prdiff, "Differential Response in Recession");

l = @layout [a b; c d]
plot(pss, prdiff, pr0, pr1, layout=l)
plot!(size=(2400,1600))



## Growth Shock
u = zeros(size(Ψ,2), T+1);
u[2, shock_t] = 0.25;
x = shock_sim(u, Y, G1, impact)
pss = plot_sim(x, Y, [:y, :r, :uis, :uis_grwth])
title!(pss, "Response at Steady State");


u = zeros(size(Ψ,2), T+1);
u[1, 1:shock_t-2] .= -0.15;
x0 = shock_sim(u, Y, G1, impact)
pr0 = plot_sim(x0, Y, [:y, :r, :uis, :uis_grwth])
title!(pr0, "Base Recovery from Recession");


u = zeros(size(Ψ,2), T+1);
u[1, 1:shock_t-2] .= -0.15;
u[2, shock_t] = 0.25;
x1 = shock_sim(u, Y, G1, impact)
pr1 = plot_sim(x1, Y, [:y, :r, :uis, :uis_grwth])
title!(pr1, "Response in Recession");


prdiff = plot_sim(x1 .- x0, Y, [:y, :r, :uis, :uis_grwth])
title!(prdiff, "Differential Response in Recession");

l = @layout [a b; c d]
plot(pss, prdiff, pr0, pr1, layout=l)
plot!(size=(2400,1600))
