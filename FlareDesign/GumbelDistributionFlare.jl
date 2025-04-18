using JuMP
using Ipopt

FlareDesign=Model(solver=IpoptSolver(max_iter=5000,tol=1e-6))

#PARAMETERS
s=1000		#Senarios
mf=100000	#Mass flow (lb/h)
srand(1234);

M=46.1		#Molecular mass
Tf=760		#Temperature °R
ComF=1		#Compressibility factor
Hc=21500	#Heat of combustion (Btu/lb)
p2=14.7		#Absolute pressure (psi)
uinf=29.3	#Wind velocity (ft/s)
F=0.3		#Fraction of heat radiated 
Msteam=46	#Steam requirement (lbs/hr)
to=1		#emisivility
r=150		#Distance to flare stack ft SE tiene que cambiar A VARIABLE

#GUMBEL DISTRIBUTION
	
qm=mf*0.15*(-log.(-log.(rand(1000))))+mf	

#CALCULATED PARAMETERS
H=(qm*Hc)				#Heat liberated BTU/h
L=10.^(0.4507*log10.(H)-1.9885)		#Flare lenght
Qtot=(qm/3600)*(379.1/M)*(Tf/520)	#Vapour volume flow ft3/s
LHV=M*Hc/379.5				#Heating value BTU/scf	
umax=400				#Maximum allowable flare tip velicity for 1000<LHV
#umax=10^((LHV+1212)/850)		#Maximum allowable flare tip velicity for 300<=LHV<=1000

#VARIABLES

@variable(FlareDesign,60/12>=dft>=1/12, start =1)		#Diameter in ft

@variable(FlareDesign,600>=h>=30, start=30)			#Height of flare stack, ft

@variable(FlareDesign,umax>=uj[1:s]>=0, start=181)		#Flare tip exit velocity, ft/s
 
@variable(FlareDesign,0.9>=Ma[1:s]>=0, start=0.2)		#Mach number

@variable(FlareDesign,delx[1:s]>=0)				#Flame distortion measure, ft 

@variable(FlareDesign,dely[1:s]>=0)				#Flame distortion measure, ft

@variable(FlareDesign,hp[1:s]>=30, start=42)			#Height of flame, ft

@variable(FlareDesign,rp[1:s]>=0, start=24)			#Reference radius, ft

@variable(FlareDesign,K[1:s]>=0,start=100)			#Max allowable radiation, BTU/(h*ft^2)

@variable(FlareDesign,D[1:s]>=30, start=160)					

@variable(FlareDesign,cost[1:s]>=0)				#Flare stack cost

@variable(FlareDesign,Mcost>=0)					#expected cost

#EQUATIONS
@NLconstraint(FlareDesign, constr1[j=1:s], (Ma[j]^2)*(dft^2)*p2^2==(1.702e-5)^2*qm[j]^2*(ComF*Tf/M))

@NLconstraint(FlareDesign, constr2[j=1:s], pi*dft^2*uj[j]==4*Qtot[j])

@NLconstraint(FlareDesign, constr3[j=1:s], log(delx[j])==log(L[j]*0.9838)+0.0754*(log(uinf)-log(uj[j])))

@NLconstraint(FlareDesign, constr4[j=1:s], log(dely[j])==log(L[j]*0.0985)-0.705*(log(uinf)-log(uj[j])))

@NLconstraint(FlareDesign, constr6[j=1:s], D[j]^2*4*pi*K[j]==to*F*H[j])

@NLconstraint(FlareDesign, constr7[j=1:s], hp[j]==h+0.5*dely[j])
@NLconstraint(FlareDesign, constr8[j=1:s], rp[j]==r-0.5*delx[j])

@NLconstraint(FlareDesign, constr9[j=1:s], D[j]^2==rp[j]^2+hp[j]^2)

@NLconstraint(FlareDesign, constr10[j=1:s], cost[j]==(94.3+11.05*(dft*12)+0.906*h)^2)

@NLconstraint(FlareDesign, Mcost==sum(cost[i]*(1/s) for i=1:s))

#GUMBEL DISTRIBUTION. Probability that Radiation will be lower or equal than 2000
#prob=0.5				#probability 50%
#prob=0.8				#probability 80%
#prob=0.9				#probability 90%
prob=0.95 				#probability 95%
#prob=0.98				#probability 98%
#prob=0.999				#probability 99.9%
@variable(FlareDesign, MeanK, start=100)	#Mean of K in different scenarios
@variable(FlareDesign, stdK, start=300)		#Standar deviation of K

@variable(FlareDesign, beta)	#Mean of K in different scenarios
@variable(FlareDesign, mu)		#Standar deviation of K


@NLconstraint(FlareDesign, MeanK==sum(K[i]*(1/s) for i=1:s))
@NLconstraint(FlareDesign, stdK^2== sum((K[i]-MeanK)^2 for i=1:s)/(s-1))

@NLconstraint(FlareDesign, beta==stdK*2.449489742783178/3.1416)
@NLconstraint(FlareDesign, mu==MeanK-0.5772*beta)

@NLconstraint(FlareDesign, mu-beta*log(-log(prob))<=2000) 

#Objective FUNCTION
@NLobjective(FlareDesign,Min, Mcost)

status=solve(FlareDesign)

println("Mean Cost ", getobjectivevalue(FlareDesign))
println("Diameter ft ",getvalue(dft))
println("Height ",getvalue(h))

#NORMAL DISTRIBUTION
ExitVel=getvalue(uj)
Delx=getvalue(delx)
Dely=getvalue(dely)
Rad=getvalue(K)
prom=getvalue(mu)
BETA=getvalue(beta)
CDFRad=zeros(1000)

for n=1:s
CDFRad[n]=exp(-exp(-(Rad[n]-prom)/BETA))
end
xCDF=sortperm(Rad)
yCDF=sortperm(CDFRad)
Rad[xCDF]
CDFRad[yCDF]

using PyPlot

x = getvalue(K) 		#Values
nbins = 75 			#Number of bins
figure("Radiation",figsize=(10,10)) 	
h=plt[:hist](x,nbins,normed="True") 	#Histogram
grid("on")
xlabel(L"$\ Radiation, \ K = \frac{BTU}{h*ft^2}$",fontsize=14)
ylabel("Y",fontsize=14)
savefig("/home/javier/Documentos/FiguresPaperDistribution/GDHisRad95.jpg")

nbins = 75 			#Number of bins
figure("Inlet flow",figsize=(10,10)) 						
k=plt[:hist](qm,nbins,normed="True") 	#Histogram
grid("on")
xlabel(L"$\ Inlet \ flow, \ Q = \frac{lb}{h}$",fontsize=14)
ylabel("Y",fontsize=14)
savefig("/home/javier/Documentos/FiguresPaperDistribution/GDHisFlow95.jpg")


figure("CPD",figsize=(10,10))
plot(Rad[xCDF],CDFRad[yCDF] ,"-", label=L"$\ p (K \leq 2000) \geq 0.95 $",linewidth=1)
legend(bbox_to_anchor=(0.1, 1), loc=2,fontsize=14)
xlabel(L"$\ Radiation, \ K = \frac{BTU}{h*ft^2}$",fontsize=14)
ylabel(L"$\ Cumulative \ Probability$",fontsize=14)
xlim(1050,2500)
ylim(0,1)
grid("on")
savefig("/home/javier/Documentos/FiguresPaperDistribution/GDCPD95.jpg")

