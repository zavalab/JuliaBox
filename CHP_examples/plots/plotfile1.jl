using PyPlot

close("all")

figure(1)
t=linspace(0,7,2016)
WD=readdlm("Electric.dat", ',')
plot(t,WD)
xlabel("Time (days)")
ylabel("Electricity demand (kW)")
grid("on")
#title("Demand of electricity in 7 days")
savefig("Electric.pdf")

figure(2)
t=linspace(0,7,2016)
A=readdlm("Thermal.dat",',')
DHWS=12*A
plot(t,DHWS)
xlabel("Time (days)")
ylabel("Hot water demand (kg)")
grid("on")
#title("Demand of hot water in 7 days")
savefig("Thermal.pdf")

figure(3)
t=linspace(0,7,2016)
TAMB=readdlm("Temperature.dat",',')
plot(t,TAMB)
xlabel("Time (days)")
ylabel("Temperature (°C)")
grid("on")
#title("Local temperature in 7 days")
savefig("Temperature.pdf")

figure(4)
t=linspace(0,24,288)
CostE=readdlm("CostEd.dat",',')
CostEd=CostE[1:288]*12
plot(t,CostEd)
xlabel("Time (h)")
ylabel("Electricity cost (USD/kWh)")
grid("on")
xlim([0,24])
#title("Cost of electricity")
PyPlot.tight_layout()
savefig("Cost.pdf")

figure(5)
t=linspace(0,7,2016)
BioEP=readdlm("Bio7cf4.dat", ',')
BioCS=readdlm("Bio7csf0123.dat", ',')
plot(t,BioEP,label="Monetization")
plot(t,BioCS,label="CS")
xlabel("Time (days)")
ylabel("Biogas Demand (kW)")
#title("Consumption of biogas in 7 days")
grid("on")
ylim([-0.5, 12.5])
legend(loc="upperleft")
savefig("Bio.pdf")

figure(6)
t=linspace(0,7,2016)
GHGTEP=readdlm("GHGT7cf4.dat", ',')
GHGTCS=readdlm("GHGT7csf0123.dat", ',')
plot(t,GHGTEP,label="Monetization")
plot(t,GHGTCS,label="Compromise")
xlabel("Time (days)")
ylabel("GHG Emissions (ton)")
grid("on")
#title("Greenhouse gas emission in 7 days")
legend()
savefig("GHGT.pdf")

figure(7)
t=linspace(0,7,2016)
SWEP=readdlm("SW7cf4.dat", ',')
SWCS=readdlm("SW7csf0123.dat", ',')
plot(t,SWEP,label="Monetization")
plot(t,SWCS,label="Compromise")
xlabel("Time (days)")
ylabel("Water Demand (kg)")
#title("Consumption of water in 7 days")
grid("on")
ylim([0, 1200])
legend()

savefig("SW.pdf")

figure(8)
t=linspace(0,7,2016)
VSTEP=readdlm("VST7cf4.dat", ',')
VSTCS=readdlm("VST7csf0123.dat", ',')
plot(t,VSTEP,label="Monetization")
plot(t,VSTCS,label="Compromise")
xlabel("Time (days)")
ylabel("Storage Level (m3)")
#title("Thermal tank storage level in 7 days")
grid("on")
legend(loc="upper left")
savefig("VST.pdf")

figure(9)
t=linspace(0,7,2016)
WCHPEP=readdlm("WCHP7cf4.dat", ',')
WCHPCS=readdlm("WCHP7csf0123.dat", ',')
plot(t,WCHPEP,label="Monetization")
plot(t,WCHPCS,label="Compromise")
xlabel("Time (days)")
ylabel("CHP Power (kW)")
#title("Electricity generation from CHP in 7 days")
grid("on")
ylim([30, 110])
legend()
savefig("WCHP.pdf")

figure(10)
cost1=3803.372721
cost2=5968.828715
cost3=5854.619235
cost4=3803.419664
cost5=4104.061634
cost6=3803.372721
cost7=5968.828715
GHGT1=17.05673791
GHGT2=13.29304582
GHGT3=14.91733242
GHGT4=17.05673771
GHGT5=15.76387189
GHGT6=13.29304582
GHGT7=17.05673791
SW1=474.0137625
SW2=462.9448857
SW3=459.7285287
SW4=465.0363851
SW5=459.7514692
SW6=459.7285287
SW7=474.0137625
scatter3D(cost1,GHGT1,SW1)#color="red",label="Min cost")
scatter3D(cost2,GHGT2,SW2)#color="green",label="Min GHGT")
scatter3D(cost3,GHGT3,SW3)#,color="brown",label="Min SW")
scatter3D(cost4,GHGT4,SW4)#,color="yellow",label="Min Cost Ext")
scatter3D(cost5,GHGT5,SW5)#,color="pink",label="Compromise")
scatter3D(cost6,GHGT6,SW6)#,color="orange",label="UP")
scatter3D(cost7,GHGT7,SW7)#,label="NS")
xlabel("Cost (USD)")
ylabel("GHG Emissions (ton)")
zlabel("Water Demand (ton)")
legend(loc="upper left")
grid("on")
savefig("CS.pdf")

#figure(11)
#cost1=3803.372721
#cost2=5968.828715
#cost3=5854.619235
#cost4=3803.419664
#cost5=4104.061634
#cost6=3803.372721
#cost7=5968.828715
#GHGT1=17.05673791
#GHGT2=13.29304582
#GHGT3=14.91733242
#GHGT4=17.05673771
#GHGT5=15.76387189
#GHGT6=13.29304582
#GHGT7=17.05673791
#SW1=474.0137625
#SW2=462.9448857
#SW3=459.7285287
#SW4=465.0363851
#SW5=459.7514692
#SW6=459.7285287
#SW7=474.0137625
#scatter3D(cost1,GHGT1,SW1,color="red",label="Min cost")
#scatter3D(cost2,GHGT2,SW2,color="green",label="Min GHGT")
#scatter3D(cost3,GHGT3,SW3,color="brown",label="Min SW")
#scatter3D(cost4,GHGT4,SW4,color="yellow",label="Min Cost Ext")
#scatter3D(cost5,GHGT5,SW5,color="pink",label="Compromise")
#scatter3D(cost6,GHGT6,SW6,color="orange",label="UP")
#scatter3D(cost7,GHGT7,SW7,label="NS")
#xlabel("Cost (USD)")
#ylabel("GHG Emissions (ton)")
#zlabel("Water Demand (ton)")
#legend(loc="upper left")
#grid("on")
#savefig("CScolor.pdf")


