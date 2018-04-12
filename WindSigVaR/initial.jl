push!(LOAD_PATH, pwd())
include("simulator.jl")
using JLD

include("setup.jl")


Vm =Array(Float64, S, Nt)              #Mean wind speed [m/s]
wr_initial = Array(Float64, S, Nt)
theta_initial = Array(Float64, S, Nt)
Tgen_initial = Array(Float64, S, Nt)
Power_initial = Array(Float64, S, Nt)
Fz_initial = Array(Float64, S, Nt)
Mz_initial = Array(Float64, S, Nt)
Mr_initial = Array(Float64, S, Nt)
xfa_initial = Array(Float64, S, Nt)
vel_xfa_initial = Array(Float64, S, Nt)
MyTB_initial = Array(Float64, S, Nt)
lambda_eff_initial = Array(Float64, S, Nt)
Ct_initial = Array(Float64, S, Nt)
Cm_initial = Array(Float64, S, Nt)
wg_initial = Array(Float64, S, Nt)
slope25_initial = Array(Float64, S)
f1_initial   = Array(Float64, S, Nt)
f2_initial   = Array(Float64, S, Nt)
f3_initial   = Array(Float64, S, Nt)
f4_initial   = Array(Float64, S, Nt)
delta_wg_integral_initial   = Array(Float64, S, Nt)
bin_max_mean_initial = Array(Float64, length(Vbin))
bin_max_square_mean_initial = Array(Float64, length(Vbin))
max_square_mean_initial = Array(Float64, length(Vbin), length(seedSet))

s=1
for i in Vbin
    bm = 0
    sold = s
    for j in seedSet
    	println("speed   ", i, "seed   ", j)
    	Vm_single = wind_profile(i, j)
        Vm_single_opt=wind_profile(i, j, RealTime)
        wr, theta, Tgen, Power, Mz, Fz, Mr, xfa, vel_xfa,  MyTB, lambda_eff, Ct, Cm, wg, slope25, f1, f2, f3, delta_wg_integral = simulation(Vm_single, P_rated0, wg_ref0, wg_b0,wg_250, RealTime, SimTime)
 
	Vm[s, :] = Vm_single_opt
	wr_initial[s,:] = wr
        theta_initial[s,:] = theta
        Tgen_initial[s,:] = Tgen
        Power_initial[s,:] = Power
        Mz_initial[s,:] = Mz
        Fz_initial[s,:] = Fz
        Mr_initial[s,:] = Mr
        xfa_initial[s,:] = xfa
        vel_xfa_initial[s,:] = vel_xfa
        MyTB_initial[s,:] = MyTB
        lambda_eff_initial[s,:] = lambda_eff
        Ct_initial[s,:] = Ct
        Cm_initial[s,:] = Cm
        wg_initial[s,:] = wg

	
        delta_wg_integral_initial[s,:] = delta_wg_integral
        s = s + 1
    end
    bin_max_mean_initial[i-Vbin[1]+1] = bm
end

s = 1
for i in Vbin
    bm = bin_max_mean_initial[i-Vbin[1]+1]
    for j in seedSet
    	MyTB = Array(Float64, Nt)
        MyTB[1:end] = MyTB_initial[s,:]'
        s = s + 1
    end
end

save("Result/initial.jld", "Vm", Vm, "wr_initial", wr_initial, "theta_initial", theta_initial, "Tgen_initial", Tgen_initial, "Power_initial", Power_initial, "Fz_initial", Fz_initial, "Mz_initial", Mz_initial, "Mr_initial", Mr_initial, "xfa_initial", xfa_initial, "vel_xfa_initial", vel_xfa_initial, "MyTB_initial", MyTB_initial, "lambda_eff_initial", lambda_eff_initial, "Ct_initial", Ct_initial, "Cm_initial", Cm_initial, "wg_initial", wg_initial, "delta_wg_integral_initial", delta_wg_integral_initial)

