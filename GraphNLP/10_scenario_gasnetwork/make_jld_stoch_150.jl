using JLD2
using DelimitedFiles
using Plots

z = 0.8                                     # gas compressibility  - []
rhon = 0.72         		                # density of natural gas at normal conditions - [kg/m3]
R = 8314.0       			                # universal gas constant [J/kgmol-K]
M = 18.0    			                    # gas molar mass of natural gas [kg/kgmol]
Tgas = 293.15      		                    # reference temperature [K]
Cp = 2.34        		                    # heat capacity @ constant pressure [kJ/kg-K]
Cv = 1.85        		                    # heat capacity @ constant volume [kJ/kg-K]

# define scaling factors for optimization
ffac=(1e+6*rhon)/(24*3600)                  # from MMSCM/day to kg/s
ffac2=(3600)/(1e+4*rhon)                    # from kg/s to scmx10^-4/hr
pfac=1e+5                                   # from bar to Pa
pfac2=1e-5                                  # from Pa to bar
dfac=1e-3                                   # from mm to m
lfac=1e+3                                   # from km to m

gam = Cp/Cv       		     	            # expansion coefficient [-]
n_poly = gam
nu2 = gam*z*R*Tgas/M  			            # gas speed of sound squared
om = (gam-1.0)/gam 		     	            # aux constant [-]
c4 = (1/ffac2)*(Cp*Tgas)                    #[kW/(scmx10-4/hr)]


# create dictionaries for storing data
junction_data1  = Dict()
pipeline_data = Dict()
compressor_data = Dict()

# add data to supply junction
j1 = 1
supply_data = Dict("smin" => 0,"smax" => 200)
data = Dict(:pmax => 54,:pmin => 50,:supplies => [supply_data],:demands => [],:demand_values => [])
junction_data1[j1]  = data



# add data for a chain of pipelines and compressors
p1 = 1
pipe_length = 300000
diameter = 0.92
epsilon = 0.025
area = (1/4)*pi*diameter*diameter
lam = (2*log10(3.7*diameter/(epsilon*dfac)))^(-2)      #pipe friction coefficient ignoring reynolds number
c1 = (pfac2/ffac2)*(nu2/area)
c2 = area*(ffac2/pfac2)
c3 = area*(pfac2/ffac2)*(8*lam*nu2)/(pi*pi*diameter^5)

junction_from = 1
junction_to = 2
pipeline_data[p1] = Dict(:pipe_length => 300000,:diameter => 0.92, :area => area, :lam => lam, :c1 => c1, :c2 => c2, :c3 => c3, :epsilon => epsilon,
:from_node => junction_from,:to_node => junction_to)
junction_data1[junction_to] =  Dict(:pmax => 70,:pmin => 50,:supplies => [],:demands => [],:demand_values => [])


junction_index = 2
pipeline_index = 2
compressor_index = 1
for j = 1:11   #11 pipelines with compressors
    junction_from = junction_index

    #compressor
    junction_to = junction_index + 1
    junction_data1[junction_to] = Dict(:pmax => 70,:pmin => 50,:supplies => [],:demands => [],:demand_values => [])
    compressor_data[compressor_index] = Dict(:from_node => junction_from,:to_node => junction_to)

    #pipeline
    junction_from = junction_to
    junction_to = junction_from + 1
    junction_data1[junction_to] = Dict(:pmax => 70,:pmin => 50,:supplies => [],:demands => [],:demand_values => [])


    pipe_length = 100000
    diameter = 0.92
    area = (1/4)*pi*diameter*diameter
    lam = (2*log10(3.7*diameter/(epsilon*dfac)))^(-2)   #pipe friction coefficient ignoring reynolds number
    c1 = (pfac2/ffac2)*(nu2/area)
    c2 = area*(ffac2/pfac2)
    c3 = area*(pfac2/ffac2)*(8*lam*nu2)/(pi*pi*diameter^5)

    pipeline_data[pipeline_index] = Dict(:pipe_length => pipe_length,:diameter => diameter, :area => area, :lam => lam, :c1 => c1, :c2 => c2, :c3 => c3, :epsilon => epsilon,
    :from_node => junction_from,:to_node => junction_to)

    global junction_index = junction_to 
    global compressor_index += 1
    global pipeline_index += 1

end

# create last pipeline and demand
pend = pipeline_index
pipe_length = 300000
diameter = 0.92
area = (1/4)*pi*diameter*diameter
lam = (2*log10(3.7*diameter/(epsilon*dfac)))^(-2)   #pipe friction coefficient ignoring reynolds number
c1 = (pfac2/ffac2)*(nu2/area)
c2 = area*(ffac2/pfac2)
c3 = area*(pfac2/ffac2)*(8*lam*nu2)/(pi*pi*diameter^5)
junction_from = junction_index
junction_to = junction_index + 1
pipeline_data[pend] = Dict(:pipe_length => 300000,:diameter => 0.92, :area => area, :lam => lam, :c1 => c1, :c2 => c2, :c3 => c3, :epsilon => epsilon,
:from_node => junction_from,:to_node => junction_to)


# define junction models for 150 scenarios. 
junctions = Array{Any,1}(undef,150)
for i in 1:150
    junctions[i] = copy(junction_data1)
end

# create an array to store demand values for each scenario
dvals = Array{Any,2}(undef,(150,24))

# define baselevels for all scenarios
# create baselevel for first 60 scenarios
for i in 1:60
    dvals[i,:] = fill(.6*7*ffac*ffac2,1,24)
end


# create baselevel for next 40 scenarios
for i in 61:100
    dvals[i,:] = fill(.7*7*ffac*ffac2,1,24)
end

# create baselevels for next 50 scenarios
for i in 101:150
    dvals[i,:] = fill(.63*7*ffac*ffac2,1,24)
end


# change demands across the scenarios
for i in 1:40
    dvals[i,10:20] .= (.6+.0025*i)*7*ffac*ffac2
end

for i in 41:60
    dvals[i,8:22] .= (.6+.0025*(i-40))*7*ffac*ffac2
end

for i in 61:90
    mult = i-60
    dvals[i,10:20] .= (.7+.0025*mult)*7*ffac*ffac2
end

for i in 91:100
    mult = i-90
    dvals[i,8:22] .=(.7 + .0025*mult)*7*ffac*ffac2
end


for i in 101:140
    mult = i-100
    dvals[i,10:20] .= (.63+.001*mult)*7*ffac*ffac2
end

for i in 141:150
    mult = i - 140
    dvals[i,8:22] .= (.63 + .001*mult)*7*ffac*ffac2
end


jend = junction_to


for i in 1:150
    junctions[i][jend] = Dict(:pmax => 70,:pmin => 34,:supplies => [],:demands => [],:demand_values => [dvals[i,:]])
end

r = 1:24

plt = plot(r,dvals[1,:],legend=false)
for i in 2:150
    plot!(r,dvals[i,:])
end
xlabel!("Time (hr)")
ylabel!("Demand")
display(plt)


# To create the JLD2 file, uncomment the code below
#=
jld_file = jldopen("13_pipelines_150.jld2","w")
for i in 1:150
    jld_file["junction_data$i"] = junctions[i]
end


jld_file["pipeline_data"] = pipeline_data
jld_file["compressor_data"] = compressor_data
close(jld_file)
=#

# plot the stochastic demand scenarios
