include("electrolyzer_struct.jl")
include("save_optimization.jl")
include("formulation.jl")
using CPUTime

# Purpose of this is to just run the sensitivity analysis of the MS-elec work

# Import the data from the ERCOT ISO
DAM_10yr = DataFrame(XLSX.readtable((@__DIR__)*"/../data/ERCOT_DAM_AVG_2014-2024.xlsx", "Sheet1"))
spp_arr_DAM = Float64.(DAM_10yr[:, "SettlementPointPrice"])

DAM_10yr_pan = DataFrame(XLSX.readtable((@__DIR__)*"/../data/ERCOT_PAN_2014-2024.xlsx", "Sheet1"))
spp_arr_DAM_pan = Float64.(DAM_10yr_pan[:, "SettlementPointPrice"])

DAM_22yr_avg = vcat(spp_arr_DAM, spp_arr_DAM)
DAM_22yr_pan = vcat(spp_arr_DAM_pan, spp_arr_DAM_pan)

yearly_avg = [6.91, 6.91, 6.76, 6.88, 6.92, 6.81, 6.67, 7.18, 8.32, 8.04, 8.13].*10 # $/MWh
distrib_10 = repeat(yearly_avg, inner=8760)
distrib_10 = distrib_10[1:length(spp_arr_DAM)]

distrib_22yr = vcat(distrib_10, distrib_10)

# CASE STUDY 1: ERCOT PARTICIPATION
################################################################################
# Set up basecase electrolyzer 
θ_basecase = Electrolyzer(2.2, #MW
    0.65, # % LHV
    3.2/1000000, # V/hr
    1.9, # V
    10, # lifetime
    0.0204/500,  # start up deg (% LHV)
    3, # $/kg H2
    .05, # %
    1816*1000, # plant cost 
    250*1000, # stack cost
    0.05) #standby load 

# Run basecase of electrolyzer
# Flexible
t_base_avg = @CPUelapsed basecase_avg = run_3st_opt(θ_basecase, DAM_22yr_avg)
t_base_pan = @CPUelapsed basecase_pan = run_3st_opt(θ_basecase, DAM_22yr_pan)

# Non-Flexible
# t_base_avg_constant = @CPUelapsed basecase_avg_constant = run_3st_opt_constant(θ_basecase, DAM_22yr_avg)
# t_base_pan_constant = @CPUelapsed basecase_pan_constant = run_3st_opt_constant(θ_basecase, DAM_22yr_pan)
# t_distrib = @CPUelapsed basecase_distr_constant = run_3st_opt_constant(θ_basecase, distrib_22yr)


# Save the results 
# Degradation
save_path = (@__DIR__)*"/data/test/"
# Flexible
pop_and_save(basecase_avg, t_base_avg, θ_basecase, "../data/test/basecase_avg.JSON")
pop_and_save(basecase_pan, t_base_pan, θ_basecase, "../data/test/basecase_pan.JSON")

# Non-Flexible
# pop_and_save(basecase_avg_constant, t_base_avg_constant, θ_basecase, "../data/test/basecase_avg_constant.JSON")
# pop_and_save(basecase_pan_constant, t_base_pan_constant, θ_basecase, "../data/test/basecase_pan_constant.JSON")
# pop_and_save(basecase_distr_constant, t_distrib, θ_basecase, "../data/test/basecase_distrib_constant.JSON")
