using DifferentialEquations
using PyPlot
using Interpolations
using Distributions
using JuMP
using Ipopt
using Gurobi
using JLD
include("setup.jl")


nHours_Horizon = 24
expectedrevenue = 409723*P_nominal/1000  #452171, 409723
expectedBand = 8273*P_nominal
costBand = expectedrevenue/expectedBand


function OptimalControl(u0, t_start, soc_min, soc_max)
    Totaltime_Horizon = nHours_Horizon*60*60                                     # seconds
    Nt_FR_Horizon = round(Int, Totaltime_Horizon/dt_FR)+1                        # number of FR time step
    Nt_FR_start = round(Int, t_start/dt_FR)
    nHours_start = round(Int, t_start/3600)

    TIME_FR = 0:dt_FR:Totaltime_Horizon
    expectedrevenue = 241407*P_nominal/1000

    csp_avg0 = u0[1]
    csn_avg0 = u0[Ncp+1]
    soc0 = csn_avg0/csnmax
    delta_sei0 = u0[Ncp+Ncn+7]
    cf0 = u0[Ncp+Ncn+Nsei+8]
    fade = cf0/Qmax
    capacity_remain = 1 - fade

    dt = 2                                         # seconds
    Nt = round(Int, Totaltime_Horizon/dt)+1        # number of temporal grid points
    TIMEG = 1:(Nt)                                 # set of temporal grid points
    RealTime = 0:dt:Totaltime_Horizon              # set of real time (s) at temporal grid points
    TIMEGm = 1:Nt-1                                # set of temporal grid points minus 1
    mTIMEG = 2:Nt                                  # set of temporal grid points except 1


    m = Model(solver=GurobiSolver(OutputFlag=0))
    @variable(m, 0<=FR_band[h in 1:nHours_Horizon]<=P_nominal*maxC)                      #kw
    @variable(m, -P_nominal*maxC<=buy_from_grid[h in 1:nHours_Horizon]<=P_nominal*maxC)  #kw
    @variable(m, power[i in 1:Nt_FR_Horizon], start= signal[Nt_FR_start + i]*P_nominal*3 )
    @constraint(m, [i in 1:Nt_FR_Horizon], power[i] ==  signal[Nt_FR_start + i]*FR_band[max(1, Int(ceil(TIME_FR[i]/3600)))] + buy_from_grid[max(1, Int(ceil(TIME_FR[i]/3600)))])


    @variable(m, capacity_remain*soc_min_stop+0.01<=soc[i in 1:Nt_FR_Horizon]<=soc_max_stop*capacity_remain-0.01)
    @constraint(m, [i in 1:Nt_FR_Horizon-1], soc[i+1]==soc[i] + power[i]*2/3600/P_nominal)
    @constraint(m, soc[1] == soc0)
    @constraint(m, soc[Nt_FR_Horizon] >= soc_min*capacity_remain)
    @constraint(m, soc[Nt_FR_Horizon] <= soc_max*capacity_remain)

    @variable(m, 0<=buy_from_grid_plus[h in 1:nHours_Horizon]<=P_nominal*maxC)
    @constraint(m, [h in 1:nHours_Horizon], buy_from_grid_plus[h] >= buy_from_grid[h])
    @objective(m, :Min, -sum(FR_band[h]* (FR_price[nHours_start+h] - costBand) for h in 1:nHours_Horizon) + sum(buy_from_grid_plus[h]*grid_price[nHours_start+h] for h in 1:nHours_Horizon))

    status = JuMP.solve(m)

    println("FR_price:  ", FR_price[nHours_start+1])
    println("FR_band:  ",getvalue(FR_band[1]))
    println("buy_from_grid:  ",getvalue(buy_from_grid[1]))
    println("revenue   ", FR_price[nHours_start+1]*getvalue(FR_band[1]))
    println("cost      ", getvalue(buy_from_grid_plus[1])*grid_price[nHours_start+1])

    FR_band = getvalue(FR_band)[1]
    grid_band = getvalue(buy_from_grid)[1]
    if status != :Optimal
        FR_band = 3*P_nominal
        power_next = P_nominal*soc0 + FR_band*mean(signal[Nt_FR_start+1:Nt_FR_start+Nt_FR_Horizon])
        if power_next <= P_nominal*capacity_remain*soc_min
            grid_band = (P_nominal*capacity_remain*soc_min - power_next)/nHours_Horizon
        elseif power_next >= P_nominal*capacity_remain*soc_max
            grid_band = (P_nominal*capacity_remain*soc_max - power_next)/nHours_Horizon
        else
            grid_band = 0
        end
    end

    return FR_band, grid_band
end


MPC(u0, "modified_sim_MPC")