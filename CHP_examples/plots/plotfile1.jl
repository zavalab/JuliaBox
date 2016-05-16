using PyPlot


figure(1)
t=linspace(0,7,2016)
WD=readdlm("Electric.dat", ',')
plot(t,WD)
xlabel("Time (days)")
ylabel("Electricity demand (kW)")
title("Demand of electricity in 7 days")
savefig("Bio.pdf")


figure(2)
t=linspace(0,7,2016)
A=readdlm("Thermal.dat",',')
DHWS=12*A
plot(t,DHWS)
xlabel("Time (days)")
ylabel("Hot water demand (kg)")
title("Demand of hot water in 7 days")
savefig("DHWS.pdf")

figure(3)
t=linspace(0,7,2016)
TAMB=readdlm("Temperature.dat",',')
plot(t,TAMB)
xlabel("Time (days)")
ylabel("Temperature (°C)")
title("Local temperature in 7 days")
savefig("Temp.pdf")

figure(4)
t=linspace(0,24,288)
CostE=readdlm("CostEd.dat",',')
CostEd=CostE[1:288]*12
plot(t,CostEd)
xlabel("Time (h)")
ylabel("Electricity cost (USD/kWh)")
title("Cost of electricity")
savefig("Cost.pdf")

figure(5)
t=linspace(0,7,2016)
BioEP=readdlm("Bio7cf4.dat", ',')
BioCS=readdlm("Bio7csf0123.dat", ',')
plot(t,BioEP,label="Min CostExt")
plot(t,BioCS,label="CS")
xlabel("Time (days)")
ylabel("Consumption of biogas (kW)")
title("Consumption of biogas in 7 days")
legend(loc="upperleft")
savefig("Bio.pdf")

figure(6)
t=linspace(0,7,2016)
GHGTEP=readdlm("GHGT7cf4.dat", ',')
GHGTCS=readdlm("GHGT7csf0123.dat", ',')
plot(t,GHGTEP,label="Min CostExt")
plot(t,GHGTCS,label="CS")
xlabel("Time (days)")
ylabel("carbon dioxide emission(ton)")
title("Greenhouse gas emission in 7 days")
legend()
savefig("GHGT.pdf")

figure(7)
t=linspace(0,7,2016)
SWEP=readdlm("SW7cf4.dat", ',')
SWCS=readdlm("SW7csf0123.dat", ',')
plot(t,SWEP,label="Min CostExt")
plot(t,SWCS,label="CS")
xlabel("Time (days)")
ylabel("Consumption of water (Kg)")
title("Consumption of water in 7 days")
legend()
savefig("SW.pdf")

figure(8)
t=linspace(0,7,2016)
VSTEP=readdlm("VST7cf4.dat", ',')
VSTCS=readdlm("VST7csf0123.dat", ',')
plot(t,VSTEP,label="Min CostExt")
plot(t,VSTCS,label="CS")
xlabel("Time (days)")
ylabel("Storage level (m3)")
title("Thermal tank storage level in 7 days")
legend(loc="upper left")
savefig("VST.pdf")

figure(9)
t=linspace(0,7,2016)
WCHPEP=readdlm("WCHP7cf4.dat", ',')
WCHPCS=readdlm("WCHP7csf0123.dat", ',')
plot(t,WCHPEP,label="Min CostExt")
plot(t,WCHPCS,label="CS")
xlabel("Time (days)")
ylabel("Generation of electricity from CHP (kW)")
title("Electricity generation from CHP in 7 days")
legend()
savefig("WCHP.pdf")

figure(10)
cost1=3803.372721
GHGT1=17.05673791
SW1=474013.7625
scatter3D(cost1,GHGT1,SW1,color="red",label="Min cost")

cost2=5968.828715
cost3=5854.619235
cost4=3803.419664
cost5=4104.061634
cost6=3803.372721
cost7=5968.828715
GHGT2=13.29304582
GHGT3=14.91733242
GHGT4=17.05673771
GHGT5=15.76387189
GHGT6=13.29304582
GHGT7=17.05673791
SW2=462944.8857
SW3=459728.5287
SW4=465036.3851
SW5=459751.4692
SW6=459728.5287
SW7=474013.7625
scatter3D(cost2,GHGT2,SW2,color="green",label="Min GHGT")
scatter3D(cost3,GHGT3,SW3,color="brown",label="Min SW")
scatter3D(cost4,GHGT4,SW4,color="yellow",label="Min Cost Ext")
scatter3D(cost5,GHGT5,SW5,color="pink",label="CS")
scatter3D(cost6,GHGT6,SW6,color="orange",label="UP")
scatter3D(cost7,GHGT7,SW7,label="NS")
xlabel("Cost (USD)")
ylabel("Emission (tons CO2)")
zlabel("Water consumption (kg)")
legend(loc="upper left")
savefig("CS.pdf")
