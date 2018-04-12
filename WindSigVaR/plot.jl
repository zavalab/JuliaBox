# solution of windTurbine problem
# subject to load constraint
# compute optimal setpoint for theta and Tgen
# pitch and torque controller is not included
# Yankai Cao, Victor Zavala
# UW-Madison, 2016

using Ipopt
using Distributions
include("simulator.jl")
include("createModel_svar.jl")            # scenario model building function
using MAT
using PyPlot
using JLD


# scenarios
Vbin = 3:25                               # mean wind speed set
seedSet = 0:9                             # seed set for each mean wind speed
nseed = length(seedSet)                   # number of seed per bin
S = length(Vbin)*length(seedSet)          # number of scenarios
SCEN = 1:S                                # scenario set
Vave = 10
deltaV = 1
threshold = 6                             # load threshold [1e7 N m]


# Axes font size
fsa = 16
# Label font size
fsl = 20


d = load("Result/initialCVar.jld")
Vm = d["Vm"]
wr_initial = d["wr_initial"]
theta_initial = d["theta_initial"]
Tgen_initial = d["Tgen_initial"]
Power_initial = d["Power_initial"]
Fz_initial = d["Fz_initial"]
Mz_initial = d["Mz_initial"]
Mr_initial = d["Mr_initial"]
xfa_initial = d["xfa_initial"]
vel_xfa_initial = d["vel_xfa_initial"]
MyTB_initial = d["MyTB_initial"]
local_max_initial = d["local_max_initial"]


Loadnbins = 40:5:200
Powernbins = 0:0.25:9
fig = figure("CVaRLoad",figsize=(10,10))
plt[:hist](local_max_initial*10, bins = Loadnbins)   #1e7 to 1e6
xlabel(L"Load $y_L$ (MNm)",fontsize=fsl)
xlim(40, 200)
ylim(0, 140)
ax=gca()
setp(ax[:get_yticklabels](),fontsize=fsa)
setp(ax[:get_xticklabels](),fontsize=fsa)
grid("on")
legend()
ax[:legend](loc="best", fontsize=fsa)
savefig("WindCVaRLoad.pdf")


Power= zeros(S)
for s in SCEN
    Power[s] = mean(Power_initial[s,:]/1e6)
end
fig = figure("CVaRPower",figsize=(10,10))
plt[:hist](Power, bins = Powernbins)
xlabel("Power $y_P$ (MW)", fontsize=fsl)
xlim(0, 9)
ylim(0, 70)
ax=gca()
setp(ax[:get_yticklabels](),fontsize=fsa)
setp(ax[:get_xticklabels](),fontsize=fsa)
grid("on")
legend()
ax[:legend](loc="best", fontsize=fsa)
savefig("WindCVaRPower.pdf")




d = load("Result/initialSVar.jld")
Power_initial = d["Power_initial"]
MyTB_initial = d["MyTB_initial"]
local_max_initial = d["local_max_initial"]

fig = figure("SVaRLoad",figsize=(10,10))
plt[:hist](local_max_initial*10, bins = Loadnbins)   #1e7 to 1e6
xlabel(L"Load $y_L$ (MNm)",fontsize=fsl)
xlim(40, 200)
ylim(0, 140)
ax=gca()
setp(ax[:get_yticklabels](),fontsize=fsa)
setp(ax[:get_xticklabels](),fontsize=fsa)
grid("on")
legend()
ax[:legend](loc="best", fontsize=fsa)
savefig("WindSVaRLoad.pdf")

Power= zeros(S)
for s in SCEN
    Power[s] = mean(Power_initial[s,:]/1e6)
end
fig = figure("SVaRPower",figsize=(10,10))
plt[:hist](Power, bins = Powernbins)
xlabel("Power $y_P$ (MW)", fontsize=fsl)
xlim(0, 9)
ylim(0, 70)
ax=gca()
setp(ax[:get_yticklabels](),fontsize=fsa)
setp(ax[:get_xticklabels](),fontsize=fsa)
grid("on")
legend()
ax[:legend](loc="best", fontsize=fsa)
savefig("WindSVaRPower.pdf")