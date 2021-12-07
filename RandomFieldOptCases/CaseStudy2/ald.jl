using InfiniteOpt, Plots, Interpolations, GaussianRandomFields, Random, Statistics, KNITRO

# Set the parameters
L = 500 # μm
ξ_mean = 2*10^-4
Tp = 200 + 273.15 #K
mp = 323.63*1.66*10^-27 #kg
σp = 2.86 #μm^-2
ϵ = .33
τ = 8.28
r = 6*10^-1 #μm
Ao = 1.1 *10^3 #m^2/kg
Vp = 7*10^-5 #m^3/kg
sbar = Ao/Vp * 10^-6 #1/μm
kb = 1.380649*10^-11 #boltzmann constant (with μm)
so = 1/σp #μm^2
vth = (2*kb*Tp/pi/mp)^.5 #μm/s
D = ϵ/τ * 2/3 * r * (8*kb*Tp/pi/mp)^(.5) #μm^2/s
γ = sbar*1/4*vth #s^-1 #precursor  species  reaction  constant
η = so*1/4*vth #μm^3 s^-1 #site reaction constant

# Setpoint parameters
"""
a adjusts amplitude
b adjusts y offset
c adjusts x offset
d adjusts slope 
"""
a, b, c, d = 0, 1, 500, .125 
coverage_setpoint(x) = ((a - b) / (1 + exp(x - (c / 2))^d)) + b

# Set the sampling parameters
Random.seed!(69)
xs = collect(range(0, L, length = 100))
num_ts = 10
num_samples = 10

# Setup the random fields
grf = GaussianRandomField(CovarianceFunction(1, SquaredExponential(300, σ=0.0008)), Spectral(), xs)
shifted_samples = [max.(sample(grf) .+ ξ_mean, 0.00001) for i in 1:num_samples]
interps = map(s -> LinearInterpolation(xs, s), shifted_samples)

# Setup the model
m = InfiniteModel(KNITRO.Optimizer)
@infinite_parameter(m, t in [0, 10], num_supports = num_ts)
@infinite_parameter(m, x in [0, L], supports = xs)

@variable(m, 0 <= yθ[i = 1:num_samples] <= 1, Infinite(t, x), start = 1)
@variable(m, 0 <= zp <= .2, start = 0.039)
@variable(m, 0 <= yp[i = 1:num_samples], Infinite(t, x), start = 0.039)

@parameter_function(m, ξ[i = 1:num_samples] == x -> interps[i](x))
@parameter_function(m, yθ_setpoint == coverage_setpoint(x))

@objective(m, Min, 1 / num_samples * sum(∫((yθ[i](10, z) - yθ_setpoint)^2, x) for i in 1:num_samples))

@constraint(m, [i = 1:num_samples], yθ[i](0, x) == 1)
@constraint(m, [i = 1:num_samples], yp[i](t, 0) == zp, DomainRestrictions(t => [supports(t)[2], 10]))
@constraint(m, [i = 1:num_samples], yp[i](0, x) == 0)
@constraint(m, [i = 1:num_samples], ∂(yp[i], t)(t, L) == 0)
@constraint(m, c1[i = 1:num_samples], ∂(yp[i], t) == D * @∂(yp[i], x^2) - γ * ξ[i] * yp[i] * yθ[i])
@constraint(m, c2[i = 1:num_samples], ∂(yθ[i], t) == -η * ξ[i] * yp[i] * yθ[i])

# Solve and extract the results
optimize!(m)
yθs = value.(yθ, ndarray = true)
yps = value.(yp, ndarray = true)
ts = value(t)
xs = value(x)

# Plot the final coverage profile
mean_θ = mean(yθs)[end, :]
std_θ = std(yθs)[end, :]
plot(xs, 1 .- mean_θ, color = :green, label = "Solution")
plot!(xs, 1 .- min.(mean_θ + std_θ, 1), color = :red, linestyle = :dash, label = "+σ")
plot!(xs, 1 .- max.(mean_θ - std_θ, 0), color = :blue, linestyle = :dash, label = "-σ")
plot!(xs, 1 .- coverage_setpoint.(xs), color = :magenta, label = "Setpoint")
ylabel!("Coverage (1 - yθ)")
xlabel!("x (μm)")
