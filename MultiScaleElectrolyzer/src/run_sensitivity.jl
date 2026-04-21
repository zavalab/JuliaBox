include("electrolyzer_struct.jl")
include("save_optimization.jl")
include("formulation.jl")
using CPUTime

# Purpose of this is to just run the sensitivity analysis of the MS-elec work

# Import the data from the ERCOT ISO
DAM_10yr = DataFrame(XLSX.readtable((@__DIR__)*"/../data/ERCOT_DAM_AVG_2014-2024.xlsx", "Sheet1"))
spp_arr_DAM = Float64.(DAM_10yr[:, "SettlementPointPrice"])

# DAM_10yr_pan = DataFrame(XLSX.readtable((@__DIR__)*"/ERCOT_west_test_sorted.xlsx", "Sheet1"))
# spp_arr_DAM_pan = Float64.(DAM_10yr_pan[:, "SettlementPointPrice"])

DAM_20yr_avg = vcat(spp_arr_DAM, spp_arr_DAM)
# DAM_20yr_pan = vcat(spp_arr_DAM_pan, spp_arr_DAM_pan)


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
    250*1000,# stack cost
    0.05) #sb %


# t_basecase = @CPUelapsed basecase = run_3st_opt(θ_basecase, DAM_20yr_avg)
# pop_and_save(basecase, t_basecase ,θ_basecase, "../data/casestudy2/basecase.JSON")

# CASE STUDY 2: SENSITIVITY ANALYSIS
################################################################################

# Tornado Plot:
# Varying these parameters: 
# % LHV - 0.60, 0.70
# δ_on - 1.0 μV/hr, 5.0 μV/hr
# δ_start - 
# ψ^stack

# Semi arbitrary - DOE
# θ_l_eff = Electrolyzer(2.2,
#     0.60, # % LHV
#     3.2/1000000, 1.9, 10, 0.0204/500,3,.05, 1816*1000, 250*1000, 0.05)
# θ_h_eff = Electrolyzer(2.2,
#     0.70, # % LHV
#     3.2/1000000, 1.9, 10, 0.0204/500,3,.05, 1816*1000, 250*1000, 0.05)

# D. Lim et al 2021
θ_l_don = Electrolyzer(2.2, 0.65, 
    1.0/1000000, 
    1.9, 10, 0.0204/500,3,.05, 1816*1000, 250*1000, 0.05)
θ_h_don = Electrolyzer(2.2, 0.65, 
    5.4/1000000, 
    1.9, 10, 0.0204/500,3,.05, 1816*1000, 250*1000, 0.05)

# Arbitrary
# θ_l_dst = Electrolyzer(2.2, 0.65, 3.2/1000000, 1.9, 10, 
#     0.0,
#     3,.05, 1816*1000, 250*1000, 0.05)
# θ_h_dst = Electrolyzer(2.2, 0.65, 3.2/1000000, 1.9, 10, 
#     0.0204*3/500,
#     3,.05, 1816*1000, 250*1000, 0.05)

# Semi arbitrary, also Krishnan et al. 2023
# θ_l_stk = Electrolyzer(2.2, 0.65, 3.2/1000000, 1.9, 10, 0.0204/500,3,.05, 1816*1000, 
#     150*1000, 0.05) #
# θ_h_stk = Electrolyzer(2.2, 0.65, 3.2/1000000, 1.9, 10, 0.0204/500,3,.05, 1816*1000, 
#     350*1000, 0.05)

θ_l_sb = Electrolyzer(2.2, 0.65,3.2/1000000, 1.9, 10, 0.0204/500,3,.05, 1816*1000, 250*1000,
    0.01) #

θ_h_sb = Electrolyzer(2.2, 0.65, 3.2/1000000, 1.9, 10, 0.0204/500,3,.05, 1816*1000, 250*1000,
    0.1) # 


# Run and save the results 
# (save in between to reduce chance of memory limitations ruining results)
# Efficiency
# t_l_eff = @CPUelapsed l_eff = run_3st_opt(θ_l_eff, DAM_20yr_avg)
# pop_and_save(l_eff,t_l_eff, θ_l_eff,  "../data/casestudy2/low_eff.JSON")

# t_h_eff = @CPUelapsed h_eff = run_3st_opt(θ_h_eff, DAM_20yr_avg)
# pop_and_save(h_eff, t_h_eff, θ_h_eff, "../data/casestudy2/high_eff.JSON")

# On Degradation
# t_l_don = @CPUelapsed l_don = run_3st_opt(θ_l_don, DAM_20yr_avg)
# pop_and_save(l_don, t_l_don, θ_l_don , "../data/casestudy2/low_don.JSON")
# t_h_don = @CPUelapsed h_don = run_3st_opt(θ_h_don, DAM_20yr_avg)
# pop_and_save(h_don, t_h_don, θ_h_don , "../data/casestudy2/high_don.JSON")

# Start Degradation
# t_l_dst = @CPUelapsed l_dst = run_3st_opt(θ_l_dst, DAM_20yr_avg)
# pop_and_save(l_dst, t_l_dst, θ_l_dst , "/../data/casestudy2/low_dst.JSON")
# t_h_dst = @CPUelapsed h_dst = run_3st_opt(θ_h_dst, DAM_20yr_avg)
# pop_and_save(h_dst, t_h_dst ,θ_h_dst, "../data/casestudy2/high_dst.JSON")

# Stack Cost
# t_l_stk = @CPUelapsed l_stk = run_3st_opt(θ_l_stk, DAM_20yr_avg)
# pop_and_save(l_stk, t_l_stk, θ_l_stk , "../data/casestudy2/low_stk.JSON")
# t_h_dst = @CPUelapsed h_stk = run_3st_opt(θ_h_stk, DAM_20yr_avg)
# pop_and_save(h_stk, t_h_dst ,θ_h_stk, "../data/casestudy2/high_stk.JSON")

# SB %
# t_l_sb = @CPUelapsed l_sb = run_3st_opt(θ_l_sb, DAM_20yr_avg, gap=0.03)
# pop_and_save(l_sb, t_l_sb, θ_l_sb , "../data/tornado/low_sb.JSON")
# t_h_sb = @CPUelapsed h_sb = run_3st_opt(θ_h_sb, DAM_20yr_avg, gap=0.03)
# pop_and_save(h_sb, t_h_sb ,θ_h_sb, "../data/tornado/high_sb.JSON")


t_no_sb = @CPUelapsed no_sb = run_2st_opt(θ_basecase, DAM_20yr_avg, gap=0.03)
pop_and_save(no_sb, t_no_sb, θ_basecase, "../data/tornado/no_sb.JSON")

################################################################################
