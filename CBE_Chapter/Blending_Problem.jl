using JuMP, Ipopt

# Define sets
W = ["w1", "w2", "w3"]; # Waste streams
T = ["tA", "tB"]; # Waste treatment facilities

# Define data
fin = [10; 20; 10]; # Waste flows (gpm)
fin = Dict(zip(W,fin));
Cin = [0.1; 0.5; 0.25]; # Waste concentration (kg/gpm)
Cin = Dict(zip(W,Cin));
fmax = [20; 30]; # Max waste flow capacity (gpm)
fmax = Dict(zip(T,fmax));
Cmax = [0.4; 0.3]; # Max waste concentration capacity (kg/gpm)
Cmax = Dict(zip(T,Cmax));
α = [50; 10] # Processing cost (usd/gpm)
α = Dict(zip(T, α));

# Define optimization problem
m = Model(Ipopt.Optimizer)
set_silent(m)
@variables m begin
    f[W,T]>=0 # Waste to treatment flow
    ftot[T]>=0 # Treatment inlet flow
    Ctot[T]>=0 # Treatment inlet concentration
end
@constraint(m, [w in W], fin[w] == sum(f[w,t] for t in T))
@constraint(m, λ[t in T], ftot[t] == sum(f[w,t] for w in W))
@constraint(m, [t in T], Ctot[t] * ftot[t] == sum(f[w,t] * Cin[w] for w in W))
@constraint(m, [t in T], Ctot[t] <= Cmax[t])
@constraint(m, [t in T], ftot[t] <= fmax[t])
@objective(m, Min, sum(α[t] * ftot[t] for t in T))

# Optimize model
optimize!(m)

# Get solution
println("ftot = ", value.(ftot), "\n")
println("Ctot = ", value.(Ctot), "\n")
println("f = ", value.(f), "\n")
println("λ = ", dual.(λ))
