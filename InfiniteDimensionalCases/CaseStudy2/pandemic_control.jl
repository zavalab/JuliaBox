using JuMP, InfiniteOpt, Ipopt, Distributions, PyPlot, Random

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

# Set the intial condition values
e0 = 1 / N
i0 = 0
r0 = 0
s0 = 1 - 1 / N

# Initialize the model
model = InfiniteModel(Ipopt.Optimizer)

# Set the infinite parameters 
@infinite_parameter(model, t ∈ [t0, tf], num_supports = 101, 
                    derivative_method = OrthogonalCollocation(2))
@infinite_parameter(model, ξ ∈ Uniform(ξ_min, ξ_max), num_supports = num_samples)
add_supports(t, extra_ts)

# Set the infinite variables 
@infinite_variable(model, s(t, ξ) ≥ 0)
@infinite_variable(model, e(t, ξ) ≥ 0)
@infinite_variable(model, i(t, ξ) ≥ 0)
@infinite_variable(model, r(t, ξ) ≥ 0)
@infinite_variable(model, si(t, ξ))
@infinite_variable(model, 0 ≤ u(t) ≤ 0.8, start = 0.2)

# Set the CVaR objective 
@hold_variable(model, λ)
@infinite_variable(model, q(t) >= 0)
@objective(model, Min, λ + 1 / α * expect(q, t))
@constraint(model, q >= u - λ)

# Set the integral objective 
# @objective(model, Min, ∫(u, t))

# Set the max objective 
# @hold_variable(model, λ)
# @objective(model, Min, λ)
# @constraint(model, λ >= u)

# Define the initial conditions
@BDconstraint(model, t == 0, s == s0)
@BDconstraint(model, t == 0, e == e0)
@BDconstraint(model, t == 0, i == i0)
@BDconstraint(model, t == 0, r == r0)

# Define the SEIR equations
@constraint(model, s_constr, deriv(s, t) == -(1 - u) * β * si)
@constraint(model, e_constr, deriv(e, t) == (1 - u) * β * si - ξ * e)
@constraint(model, i_constr, deriv(i, t) == ξ * e - γ * i)
@constraint(model, r_constr, deriv(r, t) == γ * i)
@constraint(model, si == s * i)

# Define the infection limit
@constraint(model, i <= i_max)

# Optimize and get the results
optimize!(model)
r_opt = value(r, ndarray = true) * 100 # make the population fractions into percentages
s_opt = value(s, ndarray = true) * 100
i_opt = value(i, ndarray = true) * 100
e_opt = value(e, ndarray = true) * 100
u_opt = value(u) * 100
obj_opt = objective_value(model)
ts = value(t)
ξs = value(ξ)

# Plot the results
fig, ax = plt.subplots(3,1, sharex = true)

r_mean = mean(r_opt, dims = 2)
r_std = std(r_opt, dims = 2)
ax[1].plot(ts, r_mean, "C0", label = L"$y_r(t, \xi)$")
ax[1].plot(ts, r_mean + r_std, "--C0", alpha = 0.4)
ax[1].plot(ts, r_mean - r_std, "--C0", alpha = 0.4)

s_mean = mean(s_opt, dims = 2)
s_std = std(s_opt, dims = 2)
ax[1].plot(ts, s_mean, "C1", label = L"$y_s(t, \xi)$")
ax[1].plot(ts, s_mean + s_std, "--C1", alpha = 0.4)
ax[1].plot(ts, s_mean - s_std, "--C1", alpha = 0.4)
ax[1].set_ylabel("Pop. (%)")
ax[1].legend(loc = "best")

i_mean = mean(i_opt, dims = 2)
i_std = std(i_opt, dims = 2)
ax[2].plot(ts, i_mean, "C2", label = L"$y_i(t, \xi)$")
ax[2].plot(ts, i_mean + i_std, "--C2", alpha = 0.4)
ax[2].plot(ts, i_mean - i_std, "--C2", alpha = 0.4)

e_mean = mean(e_opt, dims = 2)
e_std = std(e_opt, dims = 2)
ax[2].plot(ts, e_mean, "C3", label = L"$y_e(t, \xi)$")
ax[2].plot(ts, e_mean + e_std, "--C3", alpha = 0.4)
ax[2].plot(ts, e_mean - e_std, "--C3", alpha = 0.4)
ax[2].set_ylabel("Pop.  (%)")
ax[2].legend(loc = "best")
# ax[2].set_ylim([-2, 12])

ax[3].plot(ts, u_opt, "C4", label = L"$y_u(t)$")
ax[3].set_ylim([-2, 102])
ax[3].set_xlabel("Time (Days)")
ax[3].set_ylabel("Isolation (%)")
ax[3].legend(loc = "best")

xlim([0, 200])
fig.savefig("covid_integral.png", dpi = 300)
