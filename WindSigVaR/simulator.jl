# defined functions to do simulations and analysis
# Yankai Cao, Victor Zavala
# UW-Madison, 2016

using Distributions
using Interpolations
using ODE
using MAT
using JuMP
using Ipopt

# time
Totaltime = 10*60                      # seconds
dt = 0.5                               # time step per finite element[s]
Nt = round(Int, Totaltime/dt)+1        # number of temporal grid points
TIMEG = 1:(Nt)                         # set of temporal grid points
RealTime = 0:dt:Totaltime              # set of real time (s) at temporal grid points
TIMEGm = 1:Nt-1                        # set of temporal grid points minus 1
mTIMEG = 2:Nt                          # set of temporal grid points except 1
SimTime = 0:0.1:Totaltime              # real time set for simulation
SNt = length(SimTime)

#constant
pi = 3.14159265359
rpstorpm = 9.549297

# Parameters
Ngear = 97                             #Gear ratio generato shaft [1]
Rr = 63                                #rotor radius [m]
Area = pi * Rr * Rr                    #rotor disk area [m2]
Jr = 11776047                          #Inertia of rotor [kg m2]
C0 = 0                                 #Torque loss coefficient (constant) [N m]
C1 = 0                                 #Torque loss coefficient (linear) [1]
C2 = 0                                 #Torque loss coefficient (quadratic) [1/N/m]
P0 = 0                                 #Power loss coefficient (constant) [Watt]
P1 = 0.056                             #Power loss coefficient (linear) [1]
P2 = 0                                 #Power loss coefficient (quadratic) [1/Watt]
rho = 1.225                            #Air density At sea level and at 15 °C [kg/m3]
H =90                                  #Tower height [m]
mTe = 436750                           #Tower equivalent modal mass [kg]
cTe = 17782                            #Tower structural damping
kTe = 1810000                          #Bending stiffness
Kfafz = 1.0/mTe                        #Tower thrust force effectiveness in fore-aft movement []
wfa = sqrt(kTe/mTe)                    #tower natural frequency first mode, fore-aft movement [rad/s]
zfa = cTe/mTe/wfa/2.0                  #tower damping coefficient first mode fore-aft movement [1]

P_rated0 = 5.0*10^6                    #Rated power [W]
PM_rated0 = P_rated0/(1-P1)            #Rated mechanical power [W]
wr_rated0 = 12.1/rpstorpm              #Nominal rated rotor speed [rad/s]
wg_rated0 = wr_rated0*Ngear            #Nominal rated generator speed [rad/s] 
V_rated = 11.2                         #Rated wind speed [m/s]
V_in = 3                               #Cut-in wind speed [m/s]
V_out = 25                             #Cut-out wind speed [m/s]
xfa0 = -0.014                          #Static tower-top displacement in absence of thrust forces
lambda_opt = 8.22586354673             #Optimal tip speed ratio
cp_max =  0.5100607180223101           #Peak power coefficient
theta_min = 0                          #Minimum blade pitch angle [degree]
theta_max = 30                         #Maximum blade pitch angle [degree]
vel_theta_max = 8                      #Maximum blade pitch rate [rad/s]
vel_Tgen_max = 15000                   #maximum generator torque rate N.m/s
vel_P_rated_max =  100
vel_wg_rated_max = 1


# Parameters for pitch control 
KP = 0.01882681/5	
KI = 0.008068634/5     
KD = 0
# Parameters for torque control
slip_percentage = 10
wg_tr0 = wr_rated0*Ngear/(1+0.01*slip_percentage)
wg_250 = wr_rated0*Ngear*0.99
slope250 = PM_rated0/(wg_250)/(wg_250 - wg_tr0)
Rgn2K = 0.5*rho*Area*Rr^3*cp_max/(lambda_opt^3)/(Ngear)^3			
wg_b0 = (slope250 - sqrt(slope250*( slope250 - 4.0*Rgn2K*wg_tr0)))/(2.0*Rgn2K)   
wg_ref0 = wg_rated0


#controller law of Torque in region 2
function f1_cal( x, Rgn2K)
    Rgn2K.*x.*x
end
#controller law	of Torque in region 3
function f3_cal(x, PM_rated, wg_25)
    PM_rated./(x) 
end
#controller law	of Torque in region 2.5
function f2_cal( x, slope25, wg_b)
    slope25.*( x - wg_b)+Rgn2K*wg_b*wg_b
end


function system_w_local_control(Vm, wg_ref, wg_b, wg_25, P_rated, x0, VmTime=SimTime)
    function f(t, x)
       	PM_rated = P_rated/(1-P1)
     	(wr, xfa, vel_xfa, delta_wg_integral) = x

   	# local controller
   	wg = wr*Ngear
    	delta_wg = wg - wg_ref
    	theta = KP*delta_wg + KI*delta_wg_integral
    	theta = theta*180/pi
    	theta = max(min(theta, theta_max), theta_min)
    	itp = interpolate((VmTime,), Vm, Gridded(Linear()))
    	V = itp[t]
	slope25 = (PM_rated/wg_25 - Rgn2K*wg_b*wg_b)/(wg_25 - wg_b)	 
	if (wg<=wg_b)
            Tgen = f1_cal(wg, Rgn2K)
        elseif (wg_b<wg)&&(wg<=wg_25)
            Tgen = f2_cal(wg, slope25, wg_b)
        else
            Tgen = f3_cal(wg, PM_rated, wg_25)
        end

        Veff = V - vel_xfa
    	lambda_eff =  wr * Rr / Veff	 
        Ct = 0.1831 + 0.009459*theta + 0.09179*lambda_eff -7.645e-6*theta.^2 -
                0.008423*theta.*lambda_eff - 0.001975*lambda_eff.^2

        Cm = 0.1218 + 0.005477*theta + 0.0944*lambda_eff -7.576e-5*theta.^2 -
                0.005538*theta.*lambda_eff - 0.005738*lambda_eff.^2	      
        Fz = 0.5*rho*Veff^2*Area*Ct
    	Mz = 0.5*rho*Veff^2*Area*Rr*Cm/lambda_eff
    	Mr = Tgen + C0 + C1 * Tgen + C2 * Tgen * Tgen

    	# Dynamics
    	wr_dot = (Mz - Ngear * Mr)/Jr
    	xfa_dot = vel_xfa
    	vel_xfa_dot  = -2 * zfa * wfa * vel_xfa - wfa * wfa *xfa + Kfafz * Fz

    	delta_wg_integral_dot = delta_wg
    	if delta_wg_integral + delta_wg < 0
       	    delta_wg_integral_dot = 0
        end	 
    	[wr_dot; xfa_dot; vel_xfa_dot; delta_wg_integral_dot]
    end
    (t2, pos) = ode23(f, x0, VmTime)
    (t2, pos)
end



function get_accurate_initial_ss(V, wg_ref, wg_b, wg_25, P_rated)
    function get_initial_ss(V)
        vel_xfa = 0
        Veff = V
        if V<= V_rated
            wr = lambda_opt*V/Rr;
        else
            wr = wg_ref/Ngear;
        end
        lambda_eff =  wr * Rr / Veff;
	theta = 0
	Ct = 0.1831 + 0.009459*theta + 0.09179*lambda_eff -7.645e-6*theta.^2 -
      	 0.008423*theta.*lambda_eff - 0.001975*lambda_eff.^2
        Fz = 0.5*rho*Veff^2*Area*Ct
        xfa = Kfafz * Fz /wfa / wfa
        [wr; xfa; vel_xfa; 0]
    end
    x0 = get_initial_ss(V)
    (t, pos) = system_w_local_control(V*ones(length(SimTime)), wg_ref, wg_b, wg_25, P_rated, x0, SimTime)
    wr = map(v -> v[1], pos)
    xfa = map(v -> v[2], pos)
    vel_xfa = map(v ->v[3], pos)
    delta_wg_integral = map(v -> v[4], pos)
    [wr[end]; xfa[end]; vel_xfa[end]; delta_wg_integral[end]]
end



# this function is to do simulation with local controller
# if initial condition is not provided, we compute the steady state assuming the wind profile is constantly Vm[1]
function simulation(Vm, P_rated::Float64, wg_ref, wg_b, wg_25, OutTime=SimTime, VmTime=SimTime, x0 = nothing)
    if(x0 == nothing)
        x0 = get_accurate_initial_ss(Vm[1], wg_ref0, wg_b0, wg_250, P_rated0)
    end
    (t, pos) = system_w_local_control(Vm, wg_ref, wg_b, wg_25, P_rated, x0, VmTime)
	 
    wr = map(v -> v[1], pos)
    xfa = map(v -> v[2], pos)
    vel_xfa = map(v ->v[3], pos)
    delta_wg_integral = map(v -> v[4], pos)
    PM_rated = P_rated/(1-P1)       
    slope25 = (PM_rated/wg_25 - Rgn2K*wg_b*wg_b)/(wg_25 - wg_b)
    wg = wr*Ngear
    delta_wg = wr*Ngear - wg_ref
    theta = (KP*delta_wg + KI*delta_wg_integral)*180/pi
    theta = max(min(theta, theta_max), theta_min)
    Tgen  = Array(Float64, length(wr))	 	 
    for i = 1:length(wg)
        if (wg[i]<=wg_b)
            Tgen[i] = f1_cal(wg[i], Rgn2K)
        elseif (wg_b<wg[i])&&(wg[i]<=wg_25)
            Tgen[i] = f2_cal(wg[i], slope25, wg_b)
        else
            Tgen[i] = f3_cal(wg[i],PM_rated, wg_25)
        end
    end	
    V_itp = interpolate((VmTime,), Vm, Gridded(Linear()))
    Veff = V_itp[t] - vel_xfa
    lambda_eff =  wr * Rr ./ Veff
    Ct = 0.1831 + 0.009459*theta + 0.09179*lambda_eff -7.645e-6*theta.^2 -
      	 0.008423*theta.*lambda_eff - 0.001975*lambda_eff.^2
    Cm = 0.1218 + 0.005477*theta + 0.0944*lambda_eff -7.576e-5*theta.^2 -
           0.005538*theta.*lambda_eff - 0.005738*lambda_eff.^2
	 Fz = 0.5*rho*Veff.^2*Area.*Ct
    Mz = 0.5*rho*Veff.^2*Area*Rr.*Cm./lambda_eff
	 Power = Tgen.*wr*Ngear*(1 - P1)
    Mr = Tgen + C0 + C1 * Tgen + C2 * Tgen .* Tgen
    MyTB = H/Kfafz * (wfa .* wfa .* xfa + 2 * zfa .* wfa .* vel_xfa)	 
    f1 = f1_cal(wg, Rgn2K)
    f2 = f2_cal(wg, slope25,wg_b) 
    f3 = f3_cal(wg,PM_rated,wg_25)

    wr_itp = interpolate((t,), wr, Gridded(Linear()))
    wr = wr_itp[OutTime]
    theta_itp = interpolate((t,), theta, Gridded(Linear()))
    theta = theta_itp[OutTime]
    Tgen_itp = interpolate((t,), Tgen, Gridded(Linear()))
    Tgen = Tgen_itp[OutTime]
    Power_itp = interpolate((t,), Power, Gridded(Linear()))
    Power = Power_itp[OutTime]
    Mz_itp = interpolate((t,), Mz, Gridded(Linear()))
    Mz = Mz_itp[OutTime]
    Fz_itp = interpolate((t,), Fz, Gridded(Linear()))
    Fz = Fz_itp[OutTime]
    Mr_itp = interpolate((t,), Mr, Gridded(Linear()))
    Mr = Mr_itp[OutTime]
    xfa_itp = interpolate((t,), xfa, Gridded(Linear()))
    xfa = xfa_itp[OutTime]
    vel_xfa_itp = interpolate((t,), vel_xfa, Gridded(Linear()))
    vel_xfa = vel_xfa_itp[OutTime]
    MyTB_itp = interpolate((t,), MyTB, Gridded(Linear()))
    MyTB = MyTB_itp[OutTime]
    lambda_eff_itp = interpolate((t,), lambda_eff, Gridded(Linear()))
    lambda_eff = lambda_eff_itp[OutTime]
    Ct_itp = interpolate((t,), Ct, Gridded(Linear()))
    Ct = Ct_itp[OutTime]
    Cm_itp = interpolate((t,), Cm, Gridded(Linear()))
    Cm = Cm_itp[OutTime]
    wg_itp = interpolate((t,), wg, Gridded(Linear()))
    wg = wg_itp[OutTime]
    f1_itp = interpolate((t,), f1, Gridded(Linear()))
    f1 = f1_itp[OutTime]
    f2_itp = interpolate((t,), f2, Gridded(Linear()))
    f2 = f2_itp[OutTime]
    f3_itp = interpolate((t,), f3, Gridded(Linear()))
    f3 = f3_itp[OutTime]
    delta_wg_integral_itp = interpolate((t,), delta_wg_integral, Gridded(Linear()))
    delta_wg_integral = delta_wg_integral_itp[OutTime]	 
    wr, theta, Tgen, Power, Mz, Fz, Mr, xfa, vel_xfa,  MyTB, lambda_eff, Ct, Cm, wg, slope25, f1, f2, f3, delta_wg_integral
end


# this function is to do simulation with the profile of theta and Tgen provided
function simulation_w_input(Vm, theta_input, Tgen_input, OutTime=SimTime, inputTime=SimTime, VmTime=SimTime, x0 = nothing)
    function f(t, x)
        (wr, xfa, vel_xfa) = x
        # local controller
        theta_itp = interpolate((inputTime,), theta_input, Gridded(Linear()))
        theta = theta_itp[t]
        Tgen_itp = interpolate((inputTime,), Tgen_input, Gridded(Linear()))
        Tgen = Tgen_itp[t]
        itp = interpolate((VmTime,), Vm, Gridded(Linear()))
        V = itp[t]
        Veff = V - vel_xfa
        lambda_eff =  wr * Rr / Veff
	Ct = 0.1831 + 0.009459*theta + 0.09179*lambda_eff -7.645e-6*theta.^2 -
               0.008423*theta.*lambda_eff - 0.001975*lambda_eff.^2
        Cm = 0.1218 + 0.005477*theta + 0.0944*lambda_eff -7.576e-5*theta.^2 -
               0.005538*theta.*lambda_eff - 0.005738*lambda_eff.^2
        Fz = 0.5*rho*Veff^2*Area*Ct
        Mz = 0.5*rho*Veff^2*Area*Rr*Cm/lambda_eff
        Mr = Tgen + C0 + C1 * Tgen + C2 * Tgen * Tgen
        # Dynamics
        wr_dot = (Mz - Ngear * Mr)/Jr
        xfa_dot = vel_xfa
        vel_xfa_dot  = -2 * zfa * wfa * vel_xfa - wfa * wfa *xfa + Kfafz * Fz
        [wr_dot; xfa_dot; vel_xfa_dot]
    end
    x0 = get_accurate_initial_ss(Vm[1], wg_ref0, wg_b0, wg_250, P_rated0)[1:3]
    (t, pos) = ode23(f, x0, VmTime)
    wr = map(v -> v[1], pos)
    xfa = map(v -> v[2], pos)
    vel_xfa = map(v ->v[3], pos)
    theta_itp = interpolate((inputTime,), theta_input, Gridded(Linear()))
    Tgen_itp = interpolate((inputTime,), Tgen_input, Gridded(Linear()))	
    theta = theta_itp[t]
    Tgen  = Tgen_itp[t]
    wg = wr*Ngear
    V_itp = interpolate((VmTime,), Vm, Gridded(Linear()))
    Veff = V_itp[t] - vel_xfa
    lambda_eff =  wr * Rr ./ Veff
    Ct = 0.1831 + 0.009459*theta + 0.09179*lambda_eff -7.645e-6*theta.^2 -
         0.008423*theta.*lambda_eff - 0.001975*lambda_eff.^2
    Cm = 0.1218 + 0.005477*theta + 0.0944*lambda_eff -7.576e-5*theta.^2 -
         0.005538*theta.*lambda_eff - 0.005738*lambda_eff.^2
    Fz = 0.5*rho*Veff.^2*Area.*Ct
    Mz = 0.5*rho*Veff.^2*Area*Rr.*Cm./lambda_eff
    Power = Tgen.*wr*Ngear*(1 - P1)
    Mr = Tgen + C0 + C1 * Tgen + C2 * Tgen .* Tgen
    MyTB = H/Kfafz * (wfa .* wfa .* xfa + 2 * zfa .* wfa .* vel_xfa)

    wr_itp = interpolate((t,), wr, Gridded(Linear()))
    wr = wr_itp[OutTime]
    Power_itp = interpolate((t,), Power, Gridded(Linear()))
    Power = Power_itp[OutTime]
    Mz_itp = interpolate((t,), Mz, Gridded(Linear()))
    Mz = Mz_itp[OutTime]
    Fz_itp = interpolate((t,), Fz, Gridded(Linear()))
    Fz = Fz_itp[OutTime]
    Mr_itp = interpolate((t,), Mr, Gridded(Linear()))
    Mr = Mr_itp[OutTime]
    xfa_itp = interpolate((t,), xfa, Gridded(Linear()))
    xfa = xfa_itp[OutTime]
    vel_xfa_itp = interpolate((t,), vel_xfa, Gridded(Linear()))
    vel_xfa = vel_xfa_itp[OutTime]
    MyTB_itp = interpolate((t,), MyTB, Gridded(Linear()))
    MyTB = MyTB_itp[OutTime]
    lambda_eff_itp = interpolate((t,), lambda_eff, Gridded(Linear()))
    lambda_eff = lambda_eff_itp[OutTime]
    Ct_itp = interpolate((t,), Ct, Gridded(Linear()))
    Ct = Ct_itp[OutTime]
    Cm_itp = interpolate((t,), Cm, Gridded(Linear()))
    Cm = Cm_itp[OutTime]
    wg_itp = interpolate((t,), wg, Gridded(Linear()))
    wg = wg_itp[OutTime]
    wr, Power, Mz, Fz, Mr, xfa, vel_xfa,  MyTB, lambda_eff, Ct, Cm, wg
end


# get wind profile with a certain mean wind speed and seed
function wind_profile(speed::Int, seed::Int, OutTime=SimTime)
    if speed < 10
        speed = string("0", speed)
    else
        speed = string(speed)
    end
    if seed < 10
        seed = string("00", seed)
    else
        seed = string("0", seed)
    end
    filename = string("wind2/V_lngHC/V_lngHC_", speed, seed, ".mat")
    vars = matread(filename)
    time_sample = []
    speed_sample = []
    time_sample = vars["Time"][:,1]
    speed_sample = vars["V_lngHC"][:,1]
    itp = interpolate((collect(time_sample),), speed_sample, Gridded(Linear()))
    Vm_single = itp[OutTime+50]
    Vm_single = max(min(Vm_single, V_out), V_in) #max(Vm_single, 0)
    return Vm_single
end

# extract maximum values of each segment/stage, according to block maximum method
function select_segment_maxima(signal::AbstractArray{Float64,1}, nsegment=10, dt=collect(1:length(signal)))
    result = []
    size=round(Int, length(signal)/nsegment)
    for i in 1:nsegment
        signal_piece = signal[size*(i-1)+1: size*i]
        if i == 1
            result = [maximum(signal_piece)]
        else
            push!(result, maximum(signal_piece))
        end
    end
    return result
end

# extract local maxima according to peak over threshold method
function select_local_maxima(signal::AbstractArray{Float64,1}, dt=collect(1:length(signal)))
    threshold = mean(signal) + 1.4*std(signal)
    smallnumber = 1e-6
    slope = diff(signal)
    is_maxima = [false; (slope[1:end-1] .>=smallnumber)  & (slope[2:end].<= -smallnumber) & (signal[2:end-1] .>=threshold); false]
    return signal[is_maxima] 
end


function simulation_r_MyTB_mP(input)
    P_rated = input[1]
    wg_ref = input[2]	 
    speed = round(Int, input[3])
    seed = round(Int, input[4])	 	 
    wg_b = input[5]
    wg_25 = input[6]
    Vm_single = wind_profile(speed, seed)
    wr, theta, Tgen, Power, Mz, Fz, Mr, xfa, vel_xfa,  MyTB, lambda_eff, Ct, Cm, wg, slope25, f1, f2, f3, delta_wg_integral = simulation(Vm_single, P_rated, wg_ref, wg_b, wg_25, SimTime, SimTime)
    return [MyTB; mean(Power)]
end

#for each mean wind bin, use peak over	threhold method to compute key parameters of Gumbel distribution assuing theta and Tgen profiles are given
function local_maxMyTB_w_input(input)
    i = round(Int, input[1])
    seedStart = round(Int, input[2])
    seedEnd = round(Int, input[3])
    MyTB_list = []
    Pmean = []
    t = 0
    for j in seedStart:seedEnd
        theta_input = input[4+t*Nt: (t+1)*Nt+3]
        Tgen_input = input[4+t*Nt+length(seedSet)*Nt: (t+1)*Nt+3+length(seedSet)*Nt]
        Vm_single =	 wind_profile(i,j)
        wr, Power, Mz, Fz, Mr, xfa, vel_xfa,  MyTB, lambda_eff, Ct, Cm, wg = simulation_w_input(Vm_single, theta_input, Tgen_input, SimTime, RealTime, SimTime)
        push!(Pmean, mean(Power))
        if j == 0
            MyTB_list = MyTB
        else
            MyTB_list = vcat(MyTB_list, MyTB)
        end
	t = t + 1
    end
    local_maxima = select_local_maxima(MyTB_list)
    c = std(local_maxima)*sqrt(6)/pi
    Euler = 0.577216
    x0 = mean(local_maxima)-c*Euler
    n = length(local_maxima)/length(seedSet)
    return  [c, x0, n, mean(Pmean)]
end

#for each mean wind bin, use peak over  threhold method to compute key parameters of Gumbel distribution if local controller is used
function local_max_MyTB(input)
    i = round(Int, input[1])
    P_rated = input[2]
    wg_ref = input[3]
    seedStart = round(Int, input[4])
    seedEnd = round(Int, input[5])
    wg_b = input[6]
    wg_25 = input[7]
    MyTB_bin = Float64[]
    Pmean = []
    for j in seedStart:seedEnd
	Vm_single =	wind_profile(i, j)
        wr, theta, Tgen, Power, Mz, Fz, Mr, xfa, vel_xfa,  MyTB, lambda_eff, Ct, Cm, wg, slope25, f1, f2, f3, delta_wg_integral = simulation(Vm_single, P_rated, wg_ref, wg_b, wg_25)
        push!(Pmean, mean(Power))
	MyTB_bin = [MyTB_bin; MyTB]
    end
    local_maxima = select_local_maxima(MyTB_bin)
    c = std(local_maxima)*sqrt(6)/pi
    Euler = 0.577216
    x0 = mean(local_maxima)-c*Euler
    n = length(local_maxima)/length(seedSet)
    return  [c; x0; n; mean(Pmean)]
end


#for each mean wind bin, use block maximum method to compute key parameters of Gumbel distribution if local controller is used
function segment_max_MyTB(input)
    i = round(Int, input[1])
    P_rated = input[2]
    wg_ref = input[3]
    seedStart = round(Int, input[4])
    seedEnd = round(Int, input[5])
    nsegment = round(Int, input[6])
    wg_b = input[7]
    wg_25 = input[8]
    Pmean = []
    segment_maxima = Float64[]
    println("speed: ", i)
    for j in seedStart:seedEnd
	Vm_single = wind_profile(i, j)
 	wr, theta, Tgen, Power, Mz, Fz, Mr, xfa, vel_xfa,  MyTB, lambda_eff, Ct, Cm, wg, slope25, f1, f2, f3, delta_wg_integral = simulation(Vm_single, P_rated, wg_ref, wg_b, wg_25)
        push!(Pmean, mean(Power))
        segment_maxima = [segment_maxima; select_segment_maxima(MyTB, nsegment)]
    end
    println("segment_maxima: ", segment_maxima)
    c = std(segment_maxima)*sqrt(6)/pi
    Euler = 0.577216
    x0 = mean(segment_maxima)-c*Euler
    n = length(segment_maxima)/length(seedSet)
    return  [c; x0; n; mean(Pmean)]
end


# compute characteristic_load with year as input 
function characteristic_load(Vbin, c, x0, n, deltaP, year = 50)
    m = Model(solver=IpoptSolver(print_level = 0))
    @variable(m, F, start = 100)
    @NLconstraint(m,   sum{  (1-   exp(-exp(-(1e6*F-x0[h])/c[h]))^n[h])*deltaP[h],  h in 1:length(Vbin)}*1e7 == 3.8*50/year)
    solve(m)
    return getvalue(F)*1e6     #N.m
end
