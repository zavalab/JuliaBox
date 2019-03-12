using DifferentialEquations
using PyPlot
using Interpolations
using Distributions
using JuMP
using Ipopt
using JLD
include("setup.jl")


nHours_Horizon = 1
function OptimalControl(u0, t_start, soc_min, soc_max)
    Totaltime_Horizon = nHours_Horizon*60*60                                     # seconds
    Nt_FR_Horizon = round(Int, Totaltime_Horizon/dt_FR)+1           		 # number of FR time step
    Nt_FR_start = round(Int, t_start/dt_FR)
    nHours_start = round(Int, t_start/3600)

    TIME_FR = 0:dt_FR:Totaltime_Horizon
    expectedrevenue = 241407*P_nominal/1000        # 241407, 409723

    csp_avg0 = u0[1]
    csn_avg0 = u0[Ncp+1]
    soc0 = csn_avg0/csnmax
    delta_sei0 = u0[Ncp+Ncn+7]
    cf0 = u0[Ncp+Ncn+Nsei+8]
    fade = cf0/Qmax
    capacity_remain = 1 - fade

    dt = 2					   # seconds
    Nt = round(Int, Totaltime_Horizon/dt)+1        # number of temporal grid points
    TIMEG = 1:(Nt)                         	   # set of temporal grid points
    RealTime = 0:dt:Totaltime_Horizon              # set of real time (s) at temporal grid points
    TIMEGm = 1:Nt-1                        	   # set of temporal grid points minus 1
    mTIMEG = 2:Nt                          	   # set of temporal grid points except 1


    m = Model(solver=IpoptSolver(print_level = 0, linear_solver = "ma27", max_cpu_time = 600.0))         #default linear solver mumps
    @variable(m, 0<=FR_band[h in 1:nHours_Horizon]<=P_nominal*maxC)        	         #kw
    @variable(m, 0<=buy_from_grid[h in 1:nHours_Horizon]<=P_nominal*maxC)  #kw
    @variable(m, 0<=waste[i in 1:Nt_FR_Horizon]<=P_nominal*maxC)
    @variable(m, power[i in 1:Nt_FR_Horizon], start= signal[Nt_FR_start + i]*P_nominal*3 )
    @constraint(m, [i in 1:Nt_FR_Horizon], power[i] ==  signal[Nt_FR_start + i]*FR_band[max(1, Int(ceil(TIME_FR[i]/3600)))] + buy_from_grid[max(1, Int(ceil(TIME_FR[i]/3600)))]
    		      	 		   	    	- waste[i] )

    it0 = zeros(Nt)
    for i = 1:Nt
    	it0[i] = signal[Int(ceil((t_start+i*dt)/dt_FR))]*P_nominal*3/3.3
    end

    theta_p_guess = min(0.9, csp_avg0/cspmax)
    theta_n_guess = min(0.9, csn_avg0/csnmax)


    Un_guess = 9.99877 -9.99961*theta_n_guess.^0.5 -9.98836*theta_n_guess + 8.2024*theta_n_guess.^1.5  + 0.23584./theta_n_guess -2.03569*theta_n_guess.^(-0.5) -
         1.47266*exp.(-1.14872*theta_n_guess+2.13185)-
         9.9989*tanh.(0.60345*theta_n_guess-1.58171)
    Up_guess = 7.49983 - 13.7758*theta_p_guess.^0.5 + 21.7683*theta_p_guess - 12.6985*theta_p_guess.^1.5 + 0.0174967./theta_p_guess -0.41649*theta_p_guess.^(-0.5) -
        0.0161404*exp.(100*theta_p_guess-97.1069) +
        0.363031 * tanh.(5.89493 *theta_p_guess -4.21921)


    pot0 = max(min(Up_guess-Un_guess, V_max), V_min)


    #Continuous states
    @variable(m, 1<=csp_avg[t in 2:Nt]<=cspmax-1,start = csp_avg0)
    @variable(m, 1<=csp_s[t in TIMEG]<=cspmax-1,start = csp_avg0)
    @variable(m, 1<=csn_s[t in TIMEG]<=csnmax-1,start = csn_avg0)


    back = 0.01
    if soc_min != soc_max
        @variable(m, capacity_remain*(back+soc_min_stop)*csnmax<=csn_avg[t in 2:Nt]<=capacity_remain*(soc_max_stop-back)*csnmax,start = csn_avg0)  # be a little more conservative
        JuMP.setlowerbound(csn_avg[Nt], capacity_remain*soc_min*csnmax)
    	JuMP.setupperbound(csn_avg[Nt], capacity_remain*soc_max*csnmax)
    else
	@variable(m, capacity_remain*(back+soc_min_stop)*csnmax<=csn_avg[t in 2:Nt-1]<=capacity_remain*(soc_max_stop-back)*csnmax,start = csn_avg0)
    end


    @variable(m, iint[t in TIMEG], start= it0[t])
    @variable(m, phi_p[t in TIMEG], start= Up_guess)
    @variable(m, phi_n[t in TIMEG], start= Un_guess)
    @variable(m, pot[t in TIMEG], start = pot0)
    @variable(m, it[t in TIMEG], start= it0[t])
    @variable(m, delta_sei[t in 2:Nt], start = delta_sei0)
    @variable(m, cf0<=cf[t in 2:Nt]<=cf0+0.01*Qmax, start=0)

    # Dynamics
    #C1. Governing Equations
    #Positive electrode
    @constraint(m, csp_avg_Euler[t = 2:Nt-1], (csp_avg[t+1] - csp_avg[t])/dt == - 15 * Dp/Rpp/Rpp*(csp_avg[t+1] - csp_s[t+1]))
    @constraint(m, csp_s_Euler[t = 2:Nt], 0 == 1e3*(5*(csp_s[t] - csp_avg[t]) + Rpp*it[t]/F/Dp/ap/lp))
    t = 1
    @constraint(m, (csp_avg[t+1] - csp_avg0)/dt == - 15 * Dp/Rpp/Rpp*(csp_avg[t+1] - csp_s[t+1]))
    @constraint(m, 0 == 1e3*(5*(csp_s[t] - csp_avg0) + Rpp*it[t]/F/Dp/ap/lp))

    #Negative electrode
    if soc_min != soc_max
        @constraint(m, csn_avg_Euler[t = 2:Nt-1], (csn_avg[t+1] - csn_avg[t])/dt == - 15 * Dn/Rpn/Rpn*(csn_avg[t+1] - csn_s[t+1]))
    	@constraint(m, csn_s_Euler[t = 2:Nt], 0 == 1e3*(5*(csn_s[t] - csn_avg[t]) - Rpn*iint[t]/F/Dn/an/lnn))
	
        t = 1
        @constraint(m, (csn_avg[t+1] - csn_avg0)/dt == - 15 * Dn/Rpn/Rpn*(csn_avg[t+1] - csn_s[t+1]))
        @constraint(m, 0 == 1e3*(5*(csn_s[t] - csn_avg0) - Rpn*iint[t]/F/Dn/an/lnn))
    else
        @constraint(m, csn_avg_Euler[t = 2:Nt-2], (csn_avg[t+1] - csn_avg[t])/dt == - 15 * Dn/Rpn/Rpn*(csn_avg[t+1] - csn_s[t+1]))
        @constraint(m, csn_s_Euler[t = 2:Nt-1], 0 == 1e3*(5*(csn_s[t] - csn_avg[t]) - Rpn*iint[t]/F/Dn/an/lnn))

        t = 1
        @constraint(m, (csn_avg[t+1] - csn_avg0)/dt == - 15 * Dn/Rpn/Rpn*(csn_avg[t+1] - csn_s[t+1]))
        @constraint(m, 0 == 1e3*(5*(csn_s[t] - csn_avg0) - Rpn*iint[t]/F/Dn/an/lnn))
        t = Nt-1
        @constraint(m, (soc_min*csnmax*capacity_remain - csn_avg[t])/dt == - 15 * Dn/Rpn/Rpn*(soc_min*csnmax*capacity_remain - csn_s[t+1]))
	t = Nt
        @constraint(m,  0 == 1e3*(5*(csn_s[t] - soc_min*csnmax*capacity_remain) - Rpn*iint[t]/F/Dn/an/lnn))
    end


    @variable(m, 0<=theta_p[t in TIMEG]<=1, start=theta_p_guess)
    @variable(m, 0<=theta_n[t in TIMEG]<=1, start= theta_n_guess)
    @constraint(m, theta_p_def[t = TIMEG], theta_p[t] == csp_s[t]/cspmax)
    @constraint(m, theta_n_def[t = TIMEG], theta_n[t] == csn_s[t]/csnmax)

    @variable(m, Up[t in TIMEG], start= Up_guess)
    @variable(m, Un[t in TIMEG], start = Un_guess)
    
    @variable(m, jp[t in TIMEG])
    @variable(m, jn[t in TIMEG])
    @NLconstraint(m, [t = TIMEG], jp[t] - (it[t]/ap/F/lp) == 0)
    @NLconstraint(m, jp_def[t = TIMEG], jp[t] == 2*kp*ce^(0.5)*(cspmax-csp_s[t])^(0.5)*csp_s[t]^(0.5)*sinh(0.5*F/R/T*(phi_p[t]-Up[t])))

    t = 1
    @NLconstraint(m, jn[t] == 2*kn*ce^(0.5)*(csnmax-csn_s[t])^(0.5)*csn_s[t]^(0.5)*sinh(0.5*F/R/T*(phi_n[t]-Un[t]+delta_sei0/Kappa_sei*it[t]/an/lnn)))
    @NLconstraint(m, jn_def[t = 2:Nt], jn[t] == 2*kn*ce^(0.5)*(csnmax-csn_s[t])^(0.5)*csn_s[t]^(0.5)*sinh(0.5*F/R/T*(phi_n[t]-Un[t]+delta_sei[t]/Kappa_sei*it[t]/an/lnn)))
    @NLconstraint(m, [t = TIMEG], 0 == jn[t] + (iint[t]/an/F/lnn))
    

    @NLconstraint(m, Up_def[t = TIMEG], Up[t] == 7.49983 - 13.7758*theta_p[t]^0.5 + 21.7683*theta_p[t] - 12.6985*theta_p[t]^1.5 + 0.0174967/theta_p[t] 
	-0.41649*theta_p[t]^(-0.5) - 0.0161404*exp(100*theta_p[t]-97.1069) + 0.363031 * tanh(5.89493 *theta_p[t] -4.21921))
    @NLconstraint(m, Un_def[t = TIMEG], Un[t] == 9.99877 -9.99961*theta_n[t]^0.5 -9.98836*theta_n[t] + 8.2024*theta_n[t]^1.5  + 0.23584/theta_n[t] 
	-2.03569*theta_n[t]^(-0.5) - 1.47266*exp(-1.14872*theta_n[t]+2.13185)- 9.9989*tanh(0.60345*theta_n[t]-1.58171))
   
    @constraint(m, pot_def[t = TIMEG], 0 == pot[t] - phi_p[t] + phi_n[t])
    @constraint(m, [t in mTIMEG], power[max(1, Int(ceil(RealTime[t]/dt_FR)))] == pot[t]*it[t])

    #C3. SEI layer Equations
    @variable(m, isei[t in TIMEG], start =  an*lnn*ksei*exp(-1*F/R/T*(Un_guess - Urefs + delta_sei0/Kappa_sei*it0[t]/an/lnn)))
    @constraint(m, sei1[t = TIMEG], 0 == -iint[t] + it[t] -isei[t])

    t=1
    @NLconstraint(m, 0 == -isei[t] + an*lnn*ksei*(exp(-1*F/R/T*(phi_n[t] - Urefs + delta_sei0/Kappa_sei*it[t]/an/lnn))))
    @NLconstraint(m, sei2[t = 2:Nt], 0 == -isei[t] + an*lnn*ksei*(exp(-1*F/R/T*(phi_n[t] - Urefs + delta_sei[t]/Kappa_sei*it[t]/an/lnn))))

    t = 1
    @constraint(m, (delta_sei[t+1] - delta_sei0)/dt == isei[t+1]*M_sei/F/rho_sei/an/lnn)
    @constraint(m, delta_sei_Euler[t = 2:Nt-1], (delta_sei[t+1] - delta_sei[t])/dt == isei[t+1]*M_sei/F/rho_sei/an/lnn)


    #C4. Charge stored
    @constraint(m, (cf[2]-cf0)/dt == isei[1]/3600)
    @constraint(m, cf_Euler[t = 2:Nt-1], (cf[t+1]-cf[t])/dt == isei[t]/3600)



    @objective(m, :Min, -sum(FR_band[h]* FR_price[nHours_start+h] for h in 1:nHours_Horizon) + sum(buy_from_grid[h]*grid_price[nHours_start+h] for h in 1:nHours_Horizon) + 
    		expectedrevenue/(1-soc_retire) * (cf[Nt]-cf0)/Qmax)

    status = JuMP.solve(m)

    println("FR_price:  ", FR_price[nHours_start+1])
    println("FR_band:  ",getvalue(FR_band[1]))
    println("buy_from_grid:  ",getvalue(buy_from_grid[1]))
 

    println("revenue   ", FR_price[nHours_start+1]*getvalue(FR_band[1]))
    println("cost      ", getvalue(buy_from_grid)*grid_price[nHours_start+1])
    println("penalty   ", expectedrevenue/(1-soc_retire) * (getvalue(cf[Nt])-cf0)/Qmax)

    println("capacity fade:  ",(getvalue(cf[Nt])-cf0)/Qmax)

    FR_band = getvalue(FR_band)[1]
    grid_band = getvalue(buy_from_grid)[1]
    waste = getvalue(waste)[1:Nt_FR_hour+1]
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
	if grid_band >= 0
	    waste[1:end] = 0
	else
	    waste[1:end] = -grid_band	    
	    grid_band = 0
	end
    end

    return FR_band, grid_band, waste
end


MPC(u0, "MPC_flexible")
