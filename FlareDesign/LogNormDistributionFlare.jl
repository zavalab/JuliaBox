using JuMP
using Ipopt

FlareDesign=Model(solver=IpoptSolver(tol=1e-6))

#PARAMETERS
s=1000		#Number of senarios
mf=10000	#inlet mass flow (lb/h)
M=46.1		#Molecular mass
Tf=760		#Temperature °R
ComF=1		#Compressibility factor
Hc=21500	#Heat of combustion (Btu/lb)
p2=14.7		#Absolute pressure (psi)
uinf=29.3	#Wind velocity (ft/s)
F=0.3		#Fraction of heat radiated 
to=1		#Emisivility
r=150		#Distance to flare stack (ft)
LHV=M*Hc/379.5				#Heating value BTU/scf	
umax=400				#Maximum allowable flare tip velicity for 1000<LHV
#umax=10^((LHV+1212)/850)		#Maximum allowable flare tip velicity for 300<=LHV<=1000
kmax=2000				#Maximum allowable radiation, BTU/(h*ft^2)

#Random number generation with log-normal distribution
srand(1234);
qm=exp.(randn(1000)*0.3+randn(1000))*mf
H=(qm*Hc)				#Heat liberated BTU/h
L=10.^(0.4507*log10.(H)-1.9885)		#Flare lenght
Qtot=(qm/3600)*(379.1/M)*(Tf/520)	#Vapour volume flow ft3/s

#VARIABLES
@variable(FlareDesign,60/12>=dft>=1/12, start =1)		#Diameter, ft

@variable(FlareDesign,600>=h>=30, start=30)			#Height of flare stack, ft

@variable(FlareDesign,umax>=uj[1:s]>=0, start=181)		#Flare tip exit velocity, ft/s
 
@variable(FlareDesign,0.9>=Ma[1:s]>=0, start=0.2)		#Mach number

@variable(FlareDesign,delx[1:s]>=0)				#Flame distortion measure, ft 

@variable(FlareDesign,dely[1:s]>=0)				#Flame distortion measure, ft

@variable(FlareDesign,hp[1:s]>=30, start=42)			#Height of flame, ft

@variable(FlareDesign,rp[1:s]>=0, start=24)			#Reference radius, ft

@variable(FlareDesign,K[1:s]>=0,start=100)			#Radiation, BTU/(h*ft^2)

@variable(FlareDesign,D[1:s]>=30, start=160)			#Distance from the reference point to flame, ft		

@variable(FlareDesign,cost>=0)					#Flare stack cost, USD 
@variable(FlareDesign, MeanK>=0, start=6)			#Mean of K in all scenarios

@variable(FlareDesign, stdK>=0, start=0.3)			#Standar deviation of K

@variable(FlareDesign,Klog[1:s], start=6)			#log(radiation) 

#EQUATIONS
@NLconstraint(FlareDesign, constr1[j=1:s], (Ma[j]^2)*(dft^2)*p2^2==(1.702e-5)^2*qm[j]^2*(ComF*Tf/M))

@NLconstraint(FlareDesign, constr2[j=1:s], pi*dft^2*uj[j]==4*Qtot[j])

@NLconstraint(FlareDesign, constr3[j=1:s], log(delx[j])==log(L[j]*0.9838)+0.0754*(log(uinf)-log(uj[j])))

@NLconstraint(FlareDesign, constr4[j=1:s], log(dely[j])==log(L[j]*0.0985)-0.705*(log(uinf)-log(uj[j])))

@NLconstraint(FlareDesign, constr6[j=1:s], D[j]^2*4*pi*K[j]==to*F*H[j])

@NLconstraint(FlareDesign, constr7[j=1:s], hp[j]==h+0.5*dely[j])
@NLconstraint(FlareDesign, constr8[j=1:s], rp[j]==r-0.5*delx[j])

@NLconstraint(FlareDesign, constr9[j=1:s], D[j]^2==rp[j]^2+hp[j]^2)

@NLconstraint(FlareDesign, cost==(94.3+11.05*(dft*12)+0.906*h)^2)

#Standard normal distribution Z
#zprob=3.576279e-7		#Z value given a probability of 50%
#zprob=0.841621			#Z value given a probability of 80%
#zprob=1.281551			#Z value given a probability of 90%
zprob=1.644852 			#Z value given a probability of 95%
#zprob=2.053749			#Z value given a probability of 98%
#zprob=3.090232			#Z value given a probability of 99.9%

@NLconstraint(FlareDesign, constr11[j=1:s], Klog[j]==log(K[j]))

@NLconstraint(FlareDesign, MeanK==sum(Klog[i]*(1/s) for i=1:s))

@NLconstraint(FlareDesign, stdK^2== sum((Klog[i]-MeanK)^2 for i=1:s)/(s-1))

@NLconstraint(FlareDesign, MeanK+zprob*stdK<=log(2000)) 

#Objective function
@NLobjective(FlareDesign,Min, cost)
status=solve(FlareDesign)
#Results
println("Mean Cost ", getobjectivevalue(FlareDesign))
println("Diameter ft ",getvalue(dft))
println("Height ",getvalue(h))

Rad=getvalue(K)
CDFRad=zeros(1000)
PDFRad=zeros(1000)
PDFflow=zeros(1000)
CDFRad1=zeros(1000)

#CDF and Empirical CDF calculation
using Stats
CDFRad=map(ecdf(Rad),Rad)
xCDF=sortperm(Rad)
yCDF=sortperm(CDFRad)
Rad[xCDF]
CDFRad[yCDF]

for n=1:s
CDFRad1[n]=0.5+0.5*(erf((log(Rad[n])-getvalue(MeanK))/(sqrt(2)*getvalue(stdK))))
end
xCDF1=sortperm(Rad)
yCDF1=sortperm(CDFRad1)
Rad[xCDF1]
CDFRad1[yCDF1]

for n=1:s
PDFRad[n]=1/(Rad[n]*getvalue(stdK)*sqrt(3.1416*2)) * exp(-0.5*( (log(Rad[n])-getvalue(MeanK)) / getvalue(stdK) )^2)    
end
xPDF=sortperm(Rad)
yPDF=sortperm(PDFRad, rev=true)
Rad[xPDF]
PDFRad[yPDF]

for n=1:s
PDFflow[n]=1/(qm[n]*std(log.(qm))*sqrt(3.1416*2))*exp(-0.5*((log.(qm[n])-mean(log.(qm)))/std(log.(qm)))^2)
end
xPDFflow=sortperm(qm)
yPDFflow=sortperm(PDFflow, rev=true)
qm[xPDFflow]
PDFflow[yPDFflow]

maxFlow=maximum(qm)

#Graphs of results
using PyPlot
x = getvalue(K) 		
nbins = 50			
figure("Radiation",figsize=(10,10)) 
h=plt[:hist](x,nbins,normed="True", color="blue") 	#Normalized histogram
plot(Rad[xPDF],PDFRad[yPDF] ,"--",linewidth=1,color="red") 
plt[:ticklabel_format](style="sci",axis="y",scilimits=(0,0))
grid("on")
xlabel(L"$\ Radiation, \ \frac{BTU}{h*ft^2}$",fontsize=16)
ylabel(L"$\ Probability $",fontsize=16)
xlim(0,6200)
ylim(0,0.003)
ax = gca()
setp(ax[:get_yticklabels](),fontsize=16)
setp(ax[:get_xticklabels](),fontsize=16)
savefig("/home/javier/Documentos/FiguresPaperDistribution/LND_His_Rad_MM_95.pdf")

nbins = 50 			
figure("Inlet Flow",figsize=(10,10)) 
k=plt[:hist](qm,nbins,normed="True",color="blue") 		#Normalized histogram
plot(qm[xPDFflow],PDFflow[yPDFflow] ,"--",linewidth=1,color="red") 
plt[:ticklabel_format](style="sci",axis="x",scilimits=(0,0))
plt[:ticklabel_format](style="sci",axis="y",scilimits=(0,0))
axis("tight")
grid("on")
xlabel(L"$\ Inlet \ flow, \frac{lb}{h}$",fontsize=16)
ylabel(L"$\ Probability$",fontsize=16)
xlim(0,maxFlow)
ax = gca()
setp(ax[:get_yticklabels](),fontsize=16)
setp(ax[:get_xticklabels](),fontsize=16)
savefig("/home/javier/Documentos/FiguresPaperDistribution/LND_His_InletFlow_95.pdf")

figure("CDF",figsize=(10,10))
plot(Rad[xCDF],CDFRad[yCDF] ,"-", label="MM",linewidth=1,color="blue")
#plot(Rad[xCDF1],CDFRad1[yCDF1] ,"-",label="CDF",linewidth=2,color="blue")
#plot(Rad[xCDF],CDFRad[yCDF] ,"-.", label="Empirical CDF",linewidth=2,color="black")
legend(bbox_to_anchor=(0.6, 0.35), loc=2,fontsize=14)
plot([2000,2000], [0,0.95] ,"--",linewidth=1, color="red")
plot([0,2000], [0.95,0.95] ,"--",linewidth=1, color="red")
axis("tight")
xlabel(L"$\ Radiation, \frac{BTU}{h*ft^2}$",fontsize=16)
ylabel(L"$\ Cumulative \ Probability$",fontsize=16)
xlim(0,6200)
ylim(0,1)
grid("on")
ax = gca()
setp(ax[:get_yticklabels](),fontsize=16)
setp(ax[:get_xticklabels](),fontsize=16)
savefig("/home/javier/Documentos/FiguresPaperDistribution/LND_CDF_MM_95.pdf")
