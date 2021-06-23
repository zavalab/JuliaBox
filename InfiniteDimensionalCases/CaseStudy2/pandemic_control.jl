using InfiniteOpt, Ipopt, Distributions, PyPlot, Random

Random.seed!(13)

# Set the SEIR parameters
γ = 0.303
β = 0.727
N = 1e5
ξ_min = 0.1 
ξ_max = 0.6

# Set the domain information
i_max = 0.02
ϵ = 0.005
t0 = 0
tf = 200
extra_ts = [0.001, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08, 0.2, 0.4, 0.8]
num_samples = 20
α = 0.5

# Set the objective mode 
objective_mode = 3 # 1 -> integral, 2 -> max, 3 -> CVaR

# Set the intial condition values
y0 = Dict(:e => 1 / N, :i => 0, :r => 0, :s => 1 - 1 / N)

# Initialize the model
model = InfiniteModel(Ipopt.Optimizer)

# Set the infinite parameters 
@infinite_parameter(model, t ∈ [t0, tf], num_supports = 101, 
                    derivative_method = OrthogonalCollocation(2))
@infinite_parameter(model, ξ ~ Uniform(ξ_min, ξ_max), num_supports = num_samples)
add_supports(t, extra_ts)

# Set the infinite variables 
var_inds = [:s, :e, :i, :r]
@variable(model, y[var_inds] ≥ 0, Infinite(t, ξ))
@variable(model, ysi, Infinite(t, ξ))
@variable(model, 0 ≤ yu ≤ 0.8, Infinite(t), start = 0.2)

if objective_mode == 1
    # Set the integral objective 
    @objective(model, Min, 1 / (tf - t0) * ∫(yu, t))
elseif objective_mode == 2
    # Set the max objective 
    @variable(model, z)
    @objective(model, Min, z)
    @constraint(model, z ≥ yu)
else
    # Set the CVaR objective 
    @variable(model, z)
    @variable(model, ym ≥ 0, Infinite(t))
    @objective(model, Min, z + 1 / (1 - α) * 𝔼(ym, t))
    @constraint(model, ym ≥ yu - z)
end

# Define the initial conditions
@constraint(model, [v ∈ var_inds], y[v](0, ξ) == y0[v])

# Define the SEIR equations
@constraints(model, begin 
    ∂(y[:s], t) == -(1 - yu) * β * ysi
    ∂(y[:e], t) == (1 - yu) * β * ysi - ξ * y[:e]
    ∂(y[:i], t) == ξ * y[:e] - γ * y[:i]
    ∂(y[:r], t) == γ * y[:i]
    ysi == y[:s] * y[:i]
end)

# Define the infection limit
@constraint(model, y[:i] ≤ i_max)

# Optimize and get the results
optimize!(model)
state_opt = value.(y, ndarray = true) * 100 # make the population fractions into percentages
control_opt = value(yu) * 100
obj_opt = objective_value(model)
ts = value(t)
ξs = value(ξ)

# Plot the results
fig, ax = plt.subplots(3,1, sharex = true)

r_mean = mean(state_opt[:r], dims = 2)
r_std = std(state_opt[:r], dims = 2)
ax[1].plot(ts, r_mean, "C0", label = L"$y_r(t, \xi)$")
ax[1].plot(ts, r_mean + r_std, "--C0", alpha = 0.4)
ax[1].plot(ts, r_mean - r_std, "--C0", alpha = 0.4)

s_mean = mean(state_opt[:s], dims = 2)
s_std = std(state_opt[:s], dims = 2)
ax[1].plot(ts, s_mean, "C1", label = L"$y_s(t, \xi)$")
ax[1].plot(ts, s_mean + s_std, "--C1", alpha = 0.4)
ax[1].plot(ts, s_mean - s_std, "--C1", alpha = 0.4)
ax[1].set_ylabel("Pop. (%)")
ax[1].legend(loc = "best")

i_mean = mean(state_opt[:i], dims = 2)
i_std = std(state_opt[:i], dims = 2)
ax[2].plot(ts, i_mean, "C2", label = L"$y_i(t, \xi)$")
ax[2].plot(ts, i_mean + i_std, "--C2", alpha = 0.4)
ax[2].plot(ts, i_mean - i_std, "--C2", alpha = 0.4)

e_mean = mean(state_opt[:e], dims = 2)
e_std = std(state_opt[:e], dims = 2)
ax[2].plot(ts, e_mean, "C3", label = L"$y_e(t, \xi)$")
ax[2].plot(ts, e_mean + e_std, "--C3", alpha = 0.4)
ax[2].plot(ts, e_mean - e_std, "--C3", alpha = 0.4)
ax[2].set_ylabel("Pop.  (%)")
ax[2].legend(loc = "best")
# ax[2].set_ylim([-2, 12])

ax[3].plot(ts, control_opt, "C4", label = L"$y_u(t)$")
ax[3].set_ylim([-2, 102])
ax[3].set_xlabel("Time (Days)")
ax[3].set_ylabel("Isolation (%)")
ax[3].legend(loc = "best")

xlim([0, 200])
# fig.savefig("covid_integral.png", dpi = 300, transparent = true)
