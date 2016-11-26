# Make several plots of data
# Yankai Cao, Siyu Chen, Luis Fuentes
# UW-Madison, 2016

push!(LOAD_PATH,pwd())
using PyPlot
include("data.jl")  #load chp data
NS=188              # number of scenarios


GD = GD/DT/1000     # convert unit from kg/5min to ton/1hour
GDm=zeros(T)  
WDm=zeros(T)
TAMBm=zeros(T)
for i=1:T
    for j=1:NS
        GDm[i]=GDm[i]+GD[i,j]
	WDm[i]=WDm[i]+WD[i,j]
        TAMBm[i]=TAMBm[i]+TAMB[i,j]
    end
    GDm[i]=GDm[i]/NS
    WDm[i]=WDm[i]/NS
    TAMBm[i]=TAMBm[i]/NS
end


# Axes font size
fsa = 16

# Label font size
fsl = 20

#GD versus TIME
TIME = (1:T)*DT
fig = figure("GD",figsize=(10,10))
for s in 1:NS
    if s == 1
        plot(TIME, GD[:,s], color="grey",label="Real realization")
    else
	plot(TIME, GD[:,s], color="grey")
    end
end
plot(TIME, GDm, color="blue", label="Average")
xlabel("Time (h)",fontsize=fsl)
ylabel("Hot water demand (Ton/h)",fontsize=fsl)
xlim(0, 24)
xticks([0:4:24])
yticks([0:5:30])
ax=gca()
setp(ax[:get_yticklabels](),fontsize=fsa)
setp(ax[:get_xticklabels](),fontsize=fsa)
legend()
grid("on")
PyPlot.tight_layout()
savefig("../Figure_Data/Hot_water_demand.pdf")
#savefig("../Figure_Data/Hot_water_demand.eps")



#WD versus TIME
fig = figure("WD",figsize=(10,10))
for s in 1:NS
    if s == 1
        plot(TIME, WD[:,s], color="grey",label="Real realization")
    else
        plot(TIME, WD[:,s], color="grey")
    end
end
plot(TIME, WDm, color="blue", label="Average")
xlabel("Time (h)",fontsize=fsl)
ylabel("Electricity demand (Kw)",fontsize=fsl)
xlim(0, 24)
xticks(collect([0:4:24]))
yticks(collect([0:500:4000]))
ax=gca()
setp(ax[:get_yticklabels](),fontsize=fsa)
setp(ax[:get_xticklabels](),fontsize=fsa)
legend()
grid("on")
PyPlot.tight_layout()
savefig("../Figure_Data/Electricity_demand.pdf")
#savefig("../Figure_Data/Electricity_demand.eps")



#TAMB versus TIME
fig = figure("TAMB",figsize=(10,10))
for s in 1:NS
    if s == 1
        plot(TIME, TAMB[:,s], color="grey",label="Real realization")
    else
        plot(TIME, TAMB[:,s], color="grey")
    end
end
plot(TIME, TAMBm, color="blue", label="Average")
xlabel("Time (h)",fontsize=fsl)
ylabel(L"Ambient temperature ($^{\circ}$C)",fontsize=fsl)
xlim(0, 24)
xticks(collect([0:4:24]))
yticks(collect([5:5:30]))
ax=gca()
setp(ax[:get_yticklabels](),fontsize=fsa)
setp(ax[:get_xticklabels](),fontsize=fsa)
legend()
grid("on")
PyPlot.tight_layout()
savefig("../Figure_Data/Ambient_temperature.pdf")
#savefig("../Figure_Data/Ambient_temperature.eps")



#CostE versus TIME
fig = figure("CostE",figsize=(10,10))
plot(TIME, CostE[1:T]/DT, color="blue")
xlabel("Time (h)",fontsize=fsl)
ylabel("Electricity price (\$/Kwh)",fontsize=fsl)
xlim(0, 24)
xticks(collect([0:4:24]))
yticks(collect([0.06:0.02:0.2]))
ax=gca()
setp(ax[:get_yticklabels](),fontsize=fsa)
setp(ax[:get_xticklabels](),fontsize=fsa)
grid("on")
PyPlot.tight_layout()
savefig("../Figure_Data/Electricity_price.pdf")
#savefig("../Figure_Data/Electricity_price.eps")
