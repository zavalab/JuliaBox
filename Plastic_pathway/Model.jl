using JuMP
using Gurobi
using DelimitedFiles
using CSV
using DataFrames

# import, define and pretreat data
node_matrix = CSV.read("node_matrix.csv", DataFrame)
tech_matrix = CSV.read("tech_matrix.csv", DataFrame)
supply_matrix = CSV.read("supply_matrix.csv", DataFrame)
demand_matrix = CSV.read("demand_matrix.csv", DataFrame)

RawMaterials = collect(skipmissing(node_matrix[:,1]))
Intermediates = collect(skipmissing(node_matrix[:,2]))
Products = collect(skipmissing(node_matrix[:,3]))
Technology = node_matrix[:,4]
Materials = vcat(RawMaterials, Intermediates, Products)

tech = tech_matrix[:,1]
input = tech_matrix[:,2]
output = tech_matrix[:,3]
gamma_input = tech_matrix[:,4]
cost_tech = unique(tech_matrix[:,5])
CO2_tech = unique(tech_matrix[:,6])

cost_raw = supply_matrix[:,2]
CO2_raw  = supply_matrix[:,3]
capacity_raw = supply_matrix[:,4]

demand_lower = demand_matrix[:,2]

arc1 = []
for j in 1:size(input)[1]
    append!(arc1, [(input[j], tech[j])])
end
arc2 = []
for j in 1:size(input)[1]
    append!(arc2, [(tech[j], output[j])])
end
arc2 = unique(arc2)
arcs1 = tuple(arc1...)
arcs2 = tuple(arc2...)
arcs = (arcs1... , arcs2...)

cost_tech_dict = Dict(Technology .=> cost_tech)
CO2_tech_dict = Dict(Technology .=> CO2_tech)

gamma_output = ones(length(arcs2))
gamma_input_dict = Dict(arcs1 .=> gamma_input)
gamma_output_dict = Dict(arcs2 .=> gamma_output)
gamma_dict = merge(gamma_input_dict, gamma_output_dict)

cost_raw_dict = Dict(RawMaterials .=> cost_raw)
CO2_raw_dict = Dict(RawMaterials .=> CO2_raw)
capacity_raw_dict = Dict(RawMaterials .=> capacity_raw)

demand_lower_dict = Dict(Products .=> demand_lower)

# carbon price (USD/tonne CO2eq)
Ctax = 0

# optimization 
MOD = JuMP.Model(optimizer_with_attributes(Gurobi.Optimizer))
@variable(MOD, 0<= f[arc in arcs])
@variable(MOD, 0<= s[R in RawMaterials])
@variable(MOD, 0<= d[P in Products])
@variable(MOD, 0<= Tech[T in Technology])

@objective(MOD, Min,
                  sum(cost_tech_dict[T] * Tech[T]
                 + CO2_tech_dict[T] * Tech[T] * Ctax/1000 for T in Technology)
                 + sum(cost_raw_dict[R] * s[R]
                 + CO2_raw_dict[R] * s[R] * Ctax/1000 for R in RawMaterials))

@constraint(MOD, ConsSupply[i=RawMaterials], s[i] - sum(f[(ii,t)] for (ii,t) in arcs if ii==i) == 0)
@constraint(MOD, CapRawMaterial[i=RawMaterials], s[i] <= capacity_raw_dict[i])

@constraint(MOD, ConsIntermediate[j=Intermediates], sum(f[(t,jj)] for (t,jj) in arcs if jj==j) - sum(f[(jj,t)] for (jj,t) in arcs if jj==j) == 0)

@constraint(MOD, ConsDemand[k=Products], sum(f[(t,kk)] for (t,kk) in arcs if kk==k) - d[k] == 0)
@constraint(MOD, CapDemand[k=Products], d[k] >= demand_lower_dict[k])
    
@constraint(MOD, ConsTechnology1[t=Technology, m=Materials; (m,t) in arcs], sum(f[(mm,tt)] for (mm,tt) in arcs if tt==t && mm==m) + gamma_dict[m,t] * Tech[t]  == 0)
@constraint(MOD, ConsTechnology2[t=Technology, m=Materials; (t,m) in arcs], -sum(f[(tt,mm)] for (tt,mm) in arcs if mm==m && tt==t) + gamma_dict[t,m] * Tech[t] == 0)

JuMP.optimize!(MOD)
