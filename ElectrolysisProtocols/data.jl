# Power-to-X Modelling and Optimization
#************************************************************************
# Input Data
#************************************************************************
using DataFrames
using XLSX
using PyCall
using Plots
digs = 2 # rounding digits
#************************************************************************

# ERCOT Real-Time Market (RTM) Data
Scenario_rt_15min = DataFrame(XLSX.readtable("ercot_rtm_July2024.xlsx","July $date"))
T = collect(1:nrow(Scenario_rt_15min)) # ERCOT RTM one-day data, length(T) = 96 = 4 * 24 
time_res_rt = 0.25 # time resolution of RTM in hour
lambda_M = Scenario_rt_15min[!,:"Bus average LMP"] # RTM market price in $/MWh

# Electrolyzer
max_load = 1.0 # maximum power load rate in %
min_load = 0.15 # minimum power load rate in %
start_cost = 15 # cold start-up cost in $/MW

A = 18.031 # slope of H2 production curve in kg/MWh
B = 2.628 # intercept of H2 production curve in kg/h
eta_full_load = A # constant production efficieny in kg/MWh

# Storage
soc_0 = 0 # initial storage state of charge (SOC) in kg

# Protocols for structured scheduling (slot length = L, period = p, duty cycle = d, capacity = c)

prot_diff_L = [
    [1   , 0.15, 1   , 0.15], # L = 4, p = 2, d = 0.5 , c = 1.0
    [0.5 , 0.15, 0.5 , 0.15], # L = 4, p = 2, d = 0.5 , c = 0.5
    [0.15, 0.15, 1   , 1   ], # L = 4, p = 4, d = 0.5 , c = 1.0    
    [0.15, 0.15, 0.5 , 0.5 ], # L = 4, p = 4, d = 0.5 , c = 0.5    
    [1   , 1   , 1   , 0.15], # L = 4, p = 4, d = 0.75, c = 1.0
    [1   , 0.15, 0.15, 0.15], # L = 4, p = 4, d = 0.25, c = 1.0
    [1   , 1   , 0.15, 0.15, 0.15, 0.15, 1   , 1   ], # L = 8, p = 8, d = 0.5 , c = 1.0
    [0.15, 0.15, 0.5 , 0.5 , 0.5 , 0.5 , 0.15, 0.15], # L = 8, p = 8, d = 0.5 , c = 0.5
    [0.15, 0.15, 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ], # L = 8, p = 8, d = 0.75, c = 0.5
    [0.5 , 0.5 , 0.15, 0.15, 0.15, 0.15, 0.15, 0.15]  # L = 8, p = 8, d = 0.25, c = 0.5
    ]

prot_L_8 = [
    [1   , 0.15, 1   , 0.15, 1   , 0.15, 1   , 0.15], # L = 8, p = 2, d = 0.5 , c = 1.0
    [0.5 , 0.15, 0.5 , 0.15, 0.5 , 0.15, 0.5 , 0.15], # L = 8, p = 2, d = 0.5 , c = 0.5
    [0.15, 0.15, 1   , 1   , 0.15, 0.15, 1   , 1   ], # L = 8, p = 4, d = 0.5 , c = 1.0    
    [0.15, 0.15, 0.5 , 0.5 , 0.15, 0.15, 0.5 , 0.5 ], # L = 8, p = 4, d = 0.5 , c = 0.5    
    [1   , 1   , 1   , 0.15, 1   , 1   , 1   , 0.15], # L = 8, p = 4, d = 0.75, c = 1.0
    [1   , 0.15, 0.15, 0.15, 1   , 0.15, 0.15, 0.15], # L = 8, p = 4, d = 0.25, c = 1.0
    [1   , 1   , 0.15, 0.15, 0.15, 0.15, 1   , 1   ], # L = 8, p = 8, d = 0.5 , c = 1.0
    [0.15, 0.15, 0.5 , 0.5 , 0.5 , 0.5 , 0.15, 0.15], # L = 8, p = 8, d = 0.5 , c = 0.5
    [0.15, 0.15, 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ], # L = 8, p = 8, d = 0.75, c = 0.5
    [0.5 , 0.5 , 0.15, 0.15, 0.15, 0.15, 0.15, 0.15]  # L = 8, p = 8, d = 0.25, c = 0.5
    ]

prot_L_4 = [
    [1   , 0.15, 1   , 0.15], # L = 4, p = 2, d = 0.5 , c = 1.0
    [0.5 , 0.15, 0.5 , 0.15], # L = 4, p = 2, d = 0.5 , c = 0.5
    [0.15, 0.15, 1   , 1   ], # L = 4, p = 4, d = 0.5 , c = 1.0    
    [0.15, 0.15, 0.5 , 0.5 ], # L = 4, p = 4, d = 0.5 , c = 0.5    
    [1   , 1   , 1   , 0.15], # L = 4, p = 4, d = 0.75, c = 1.0
    [0.5 , 0.5 , 0.5 , 0.15], # L = 4, p = 4, d = 0.75 ,c = 0.5
    [1   , 0.15, 0.15, 0.15], # L = 4, p = 4, d = 0.25, c = 1.0
    [0.5 , 0.15, 0.15, 0.15], # L = 4, p = 4, d = 0.25, c = 0.5
    [0.15, 1   , 0.15 ,1   ], # L = 4, p = 2, d = 0.5 , c = 1.0 # repeated
    [0.5 , 0.15, 0.15 ,0.5 ], # L = 4, p = 4, d = 0.5 , c = 0.5 # repeated
    ]






