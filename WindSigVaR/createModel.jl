# wind Turbine model for a single scenario without pitch and torque controller
# find optimal theta and Tgen profile
# load constraint is not added
# Yankai Cao, Victor Zavala
# UW-Madison, 2016

c = Array{Float64}(3)
c[1] = 0.1550510257216822
c[2] = 0.6449489742783178
c[3] = 1.0
a = Array{Float64}(3, 3)
a[1,1] = 0.19681547722366
a[1,2] = 0.39442431473909
a[1,3] = 0.37640306270047
a[2,1] = -0.06553542585020
a[2,2] = 0.29207341166523
a[2,3] = 0.51248582618842
a[3,1] = 0.02377097434822
a[3,2] = -0.04154875212600
a[3,3] = 0.11111111111111
ncp = 3
CP=collect(1:3)

function createModel(SCEN, prob=1)

m = Model(solver=IpoptSolver(linear_solver = "ma57")) 

#Continuous states
@variable(m, wr[s in SCEN, t in TIMEG]<=1.5*wr_rated0, start=wr_initial[s,t])                  #rotor speed [rad/s]
@variable(m, xfa[s in SCEN, t in TIMEG, j in CP], start=xfa_initial[s,t])                      #tower top displacement fore-aft [m]
@variable(m, vel_xfa[s in SCEN, t in TIMEG, j in CP], start=vel_xfa_initial[s,t])              #tower top velocity (fore-aft) [m/s]

#Continuous inputs
@variable(m, Tgen[s in SCEN,t in TIMEG]<=1.2*P_rated0/wg_ref0/1e4,start =Tgen_initial[s,t]/1e4)#Generator torque (at high speed shaft) [N m]
@variable(m, 0<=theta[s in SCEN, t in TIMEG]<=theta_max, start = theta_initial[s,t])           #collective pitch angle [rad]


#Outputs
@variable(m, Fz[s in SCEN, t in mTIMEG], start = Fz_initial[s,t]/1e5)			       #Aerodynamic thrust force [1e5 N]
@variable(m, Mz[s in SCEN, t in TIMEG], start = Mz_initial[s,t]/1e5)                           #Aerodynamic torque [1e5 N m]
@variable(m, Power[s in SCEN, t in TIMEG],   start = Power_initial[s,t]/1e4)                   #Electrical power obtained from generator [1e4 W]
@variable(m, 0<=lambda_eff[s in SCEN, t in TIMEG] <= 30, start = lambda_eff_initial[s,t])      #Effective tip to speed ratio accounting for tower movement
@NLconstraint(m, lambda_eff_definition[s = SCEN, t = TIMEG], lambda_eff[s,t]*( Vm[s,t] - vel_xfa[s,t,ncp]) ==  wr[s,t] * Rr )
@variable(m, -10<=Ct[s in SCEN, t in TIMEG]<=10, start = Ct_initial[s,t])      	       	       #Thrust coefficient [-]
@variable(m, -10<=Cm[s in SCEN, t in TIMEG]<=10, start = Cm_initial[s,t])		       #Power coefficient [-]
@constraint(m, Ct_definition[s = SCEN, t = TIMEG], Ct[s,t] ==
                 0.1831 + 0.009459*theta[s,t] + 0.09179*lambda_eff[s,t] -7.645e-6*theta[s,t]^2 -
                 0.008423*theta[s,t]*lambda_eff[s,t] - 0.001975*lambda_eff[s,t]^2 )
@constraint(m, Cm_definition[s = SCEN, t = TIMEG],    Cm[s,t] ==
                 0.1218 + 0.005477*theta[s,t] + 0.0944*lambda_eff[s,t] -7.576e-5*theta[s,t]^2 -
                 0.005538*theta[s,t]*lambda_eff[s,t] - 0.005738*lambda_eff[s,t]^2)
@NLconstraint(m, Fz_definition[s = SCEN, t = mTIMEG], Fz[s,t] == 0.5*rho*( Vm[s,t]  - vel_xfa[s,t,ncp])^2*Area*Ct[s,t]/1e5); 
@NLconstraint(m, Mz_definition[s = SCEN, t = TIMEG], Mz[s,t]*lambda_eff[s,t] == 0.5*rho*( Vm[s,t] - vel_xfa[s,t, ncp])*( Vm[s,t] - vel_xfa[s,t, ncp])*Area*Rr*Cm[s,t]/1e5); 
@NLconstraint(m, Power_definition[s = SCEN, t = TIMEG], Power[s,t] == (1-P1)*Tgen[s,t]*wr[s,t]*Ngear)
@variable(m, MyTB[s in SCEN, t in TIMEG], start = MyTB_initial[s,t]/1e7)       	      	       #Tower base moment in y-direction [1e7 N m]
@constraint(m, MyTB_definition[s = SCEN, t = TIMEG], MyTB[s,t] == H/Kfafz * (wfa * wfa * xfa[s,t,ncp] + 2 * zfa * wfa * vel_xfa[s,t,ncp])/1e7)

# Dynamics
@variable(m, vel_xfa_dot[s in SCEN, t in mTIMEG, j in CP])
@constraint(m, vel_xfa_dot_cal[s = SCEN, t = mTIMEG, j in CP], vel_xfa_dot[s,t,j]  == -2 * zfa * wfa * vel_xfa[s,t,j] - wfa * wfa *(xfa[s,t,j]) + 1e5 * Kfafz * Fz[s,t])
@constraint(m, wr_Euler[s = SCEN, t = TIMEGm], (wr[s,t+1] - wr[s,t])/dt == (1e5*Mz[s,t+1] - 1e4*Ngear * Tgen[s,t+1])/Jr)
# RK collocation
@constraint(m, xfa_RK[s = SCEN, t = TIMEGm, j in CP],  xfa[s,t+1,j] - xfa[s,t,ncp] - dt *sum{a[k,j]*vel_xfa[s,t+1,k], k in CP} == 0)
@constraint(m, vel_xfa_RK[s = SCEN, t = TIMEGm, j in CP],  vel_xfa[s,t+1,j] - vel_xfa[s,t,ncp] - dt *sum{a[k,j]*vel_xfa_dot[s,t+1,k], k in CP} == 0)


## initial state
#@constraint(m, wr_start[s = SCEN, t = [1]],  	  	wr_initial[s,t]     == wr[s,t])
#@constraint(m, xfa_start[s = SCEN, t = [1], j in CP],           xfa_initial[s,t]    == xfa[s,t, j])
#@constraint(m, vel_xfa_start[s = SCEN, t = [1], j in CP],       vel_xfa_initial[s,t]== vel_xfa[s,t, j])




return m
end

