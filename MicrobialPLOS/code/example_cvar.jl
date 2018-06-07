# Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (zavalatejeda@wisc.edu)

# Call Libraries
using Plasmo, JuMP, Ipopt, JLD
include("param.jl")

# Load data file
datasets = load("data.jld")
data = datasets["P1"]
outputpath = "../output/cvar/"

# Options --------------------
opts = Dict(:save_y => true,
            :save_e => true,
            :outputpath => outputpath,
            :solver => IpoptSolver(output_file=outputpath*"io.out",linear_solver="mumps"))

# Arguments ------------------
args = Dict(
    :min_norm=>:CVaR,       # Prediction error formulation (L1 or CVaR)
    :reg_norm=>:L2,         # Prior formulation (none, L1, or L2)
    :model =>:LV,           # Dynamic model (LV or Saturable)
    :dilution_rate => 1/20,     # Diluting medium by 1/20
    :y_sig_rel=>.05,            # Relative measurement error
    :y_sig_min=>.1,             # Absolute measurement error
    :n_total_species=>12,       # Total number of species
    :n_disc=>Dict(:M=>5,:P=>120), # Order of discretization
    :lambda=>5e1,                   # Prior coefficient
    :eps=>1e-2,                     # Small number
    :maxr=>5,                    # Maximum absolute value of parameters
    :beta=>0.9)                 # Beta value for CVaR formulation

# Add discretization information to the data dictionary
add_discretize_info!(data,args)   

# Prior means, starting values, scaling factors
args[:r_prior]=Dict((i,j)=>0.0 for i=1:args[:n_total_species] for j=0:args[:n_total_species])
args[:r_start]=Dict((i,j)=>0.0 for i=1:args[:n_total_species] for j=0:args[:n_total_species]) 
args[:y_start]=Dict((exp,i,ii,j,k)=>0.0
                    for exp in keys(data) for i in 1:length(data[exp]) for ii in 1:length(data[exp][i]) for j in 1:data[exp][i][ii][:n_species] for k in 1:data[exp][i][ii][:n_disc_time])
args[:r_scale]=1
args[:y_scale]=Dict(exp=>[Dict((ii,j,k)=>1
                               for ii =1:length(data[exp][i]) for j in 1:data[exp][i][ii][:n_species] for k in 1:data[exp][i][ii][:n_disc_time])
                          for i=1:length(data[exp])]  for exp in keys(data))
# Absolute measurement error
args[:y_sig_abs]=Dict((exp,i,ii,j,k) => max(args[:y_sig_min],data[exp][i][ii][:y][j][k])
                      for exp in keys(data) for i in 1:length(data[exp]) for ii in 1:length(data[exp][i]) for j in 1:data[exp][i][ii][:n_species] for k in 1:data[exp][i][ii][:n_time])

# Solve the problem --------------
outs = Dict()
mkpath(opts[:outputpath])
(outs[:graph],outs[:m_parent],outs[:m_children],
 outs[:node_parent],outs[:node_children]) = param(data,args)
outs[:graph].solver=opts[:solver]
t0=time_ns()
outs[:status]=Plasmo.solve(outs[:graph])
outs[:time_elapsed]=time_ns()-t0

# Save the output ----------------
output(data,args,opts,outs)

# Inference analysis -------------
# inference_analysis(data,args,opts,500)
