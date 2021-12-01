using InfiniteOpt, Gurobi, GaussianRandomFields, Interpolations, Random, PyPlot, Statistics

Random.seed!(42)

# Setup the random fields
ts = 0:0.05:24
num_samples = 100
cov = CovarianceFunction(1, SquaredExponential(1.5, σ = 0.5))
grf = GaussianRandomField(cov, Spectral(), ts)
raw_samples = [sample(grf) for i in 1:num_samples]
shifted_samples = [[s[i] + sin(ts[i] / 24 * 2pi) + 3 for i in eachindex(ts)] for s in raw_samples]
interps = map(s -> LinearInterpolation(ts, s), shifted_samples)

# Plot the random fields 
figure()
plot(ts, interps[1].(ts), label = L"$\hat{\xi}_1(t)$")
plot(ts, interps[2].(ts), label = L"$\hat{\xi}_2(t)$")
plot(ts, interps[3].(ts), label = L"$\hat{\xi}_3(t)$")
plot(ts, sin.(ts / 24 * 2pi) .+ 3, linestyle = "--", label = L"$\mu(t)$")
xlabel(L"Time $t$")
ylabel("Unit Power Cost")
legend(loc = "best")
xlim((0, 24))
# savefig("power_fields.png", dpi = 300)

# Setup the simulation and simulate
sim = InfiniteModel(Gurobi.Optimizer)
set_silent(sim)
@infinite_parameter(sim, t in [0, 24], supports = collect(ts))
@parameter_function(sim, yg[i = 1:num_samples] == t -> t >= 12 ? 1 : 0)
@parameter_function(sim, c[i = 1:num_samples] == t -> interps[i](t))
@variable(sim, yb[1:num_samples], Infinite(t))
@variable(sim, yc[1:num_samples], Infinite(t))
@constraint(sim, [i in 1:num_samples], yc[i] == yg[i] * c[i])
@constraint(sim, [i in 1:num_samples], deriv(yb[i], t) == yg[i] - 1)
@constraint(sim, [i in 1:num_samples], yb[i](0) == 50 / 2)
optimize!(sim)

# Plot the simulation
opt_yc = value.(yc)
figure()
plot(ts, mean(opt_yc), color = "C0")
plot(ts, mean(opt_yc) + std(opt_yc), color = "C0", linestyle = "--", alpha = 0.4)
plot(ts, max.(mean(opt_yc) - std(opt_yc), 0), color = "C0", linestyle = "--", alpha = 0.4)
xlabel(L"Time $t$")
ylabel(L"$f(t, \xi(t))$")
xlim((0, 24))
# savefig("stochastic_power_sim.png", dpi = 300)

# Setup the optimization problem and solve
ts = 0:0.5:24
model = InfiniteModel(Gurobi.Optimizer)
set_silent(model)
@infinite_parameter(model, t in [0, 24], supports = collect(ts))
@variable(model, yg[1:num_samples] >= 0, Infinite(t))
@variable(model, yb[1:num_samples], Infinite(t))
@variable(model, 0 <= z <= 100)
# fix(z, 78.3333, force = true) # uncomment to simulate the deterministic reponse
@parameter_function(model, c[i = 1:num_samples] == t -> interps[i](t))
@objective(model, Min, 1 / num_samples * sum(∫(yg[i] * c[i], t) for i in 1:num_samples) + 0.3z)
@constraint(model, [i in 1:num_samples], deriv(yb[i], t) == yg[i] - 1)
@constraint(model, [i in 1:num_samples], yb[i](0) == z / 2)
@constraint(model, [i in 1:num_samples], yb[i](24) == z / 2)
@constraint(model, [i in 1:num_samples], 0.2z <= yb[i])
@constraint(model, [i in 1:num_samples], yb[i] <= z)
optimize!(model)

# Retrieve the results
opt_yg = value.(yg)
opt_yb = value.(yb)
opt_z = value(z)
opt_cost = objective_value(model)

println("Expected Cost: ", opt_cost)
println("Optimal Battery Size: ", opt_z)

figure()
plot(ts, mean(opt_yg), color = "C0", label = L"$yg(t, \xi(t))$")
plot(ts, mean(opt_yg) + std(opt_yg), color = "C0", linestyle = "--", alpha = 0.4, label = "")
plot(ts, max.(mean(opt_yg) - std(opt_yg), 0), color = "C0", linestyle = "--", alpha = 0.4, label = "")
plot(ts, mean(opt_yb), color = "C1", label = L"$yb(t, \xi(t))$")
plot(ts, mean(opt_yb) + std(opt_yb), color = "C1", linestyle = "--", alpha = 0.4, label = "")
plot(ts, mean(opt_yb) - std(opt_yb), color = "C1", linestyle = "--", alpha = 0.4, label = "")
legend(loc = "best")
xlabel(L"Time $t$")
ylabel("Power")
xlim((0, 24))
# savefig("determ_power_opt.png", dpi = 300)


# Deterministic Model
model2 = InfiniteModel(Gurobi.Optimizer)
set_silent(model2)
@infinite_parameter(model2, t in [0, 24], supports = collect(ts))
@variable(model2, yg >= 0, Infinite(t))
@variable(model2, yb, Infinite(t))
@variable(model2, 0 <= z <= 100)
@parameter_function(model2, c == t -> sin(t / 24 * 2pi) + 3)
@objective(model2, Min, ∫(yg * c, t) + 0.1z)
@constraint(model2, deriv(yb, t) == yg - 1)
@constraint(model2, yb(0) == z / 2)
@constraint(model2, yb(24) == z / 2)
@constraint(model2, 0.2z <= yb)
@constraint(model2, yb <= z)
optimize!(model2)
opt_z2 = value(z)
println("Deterministic Battery Size: ", opt_z2)
