using JuMP
using Ipopt

#PARAMETERS
s=2000		#Number of senarios
mf=21000	#inlet mass flow (lb/h)
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

#Random number generation with log-normal distribution
srand(1234);
#qm=exp.(randn(s)*0.3+randn(s))*mf
qm=-mf*log.(1-rand(s))			#Random inlet mass flow (lb/h)
H=(qm*Hc)				#Heat liberated BTU/h
L=10.^(0.4507*log10.(H)-1.9885)		#Flare lenght
Qtot=(qm/3600)*(379.1/M)*(Tf/520)	#Vapour volume flow ft3/s

alpha = 0.05;           		#Probability level cvar
threshold = 2000;			#Maximum allowable radiation, BTU/(h*ft^2)

#Set CVaR formulation
FlareDesign=Model(solver=IpoptSolver(max_iter=200,tol=1e-6))

#VARIABLES
@variable(FlareDesign,60/12>=dft>=1/12, start =1)		#Diameter in ft

@variable(FlareDesign,600>=h>=30, start=30)			#Height of flare stack, ft

@variable(FlareDesign,umax>=uj[1:s]>=0, start=181)		#Flare tip exit velocity, ft/s
 
@variable(FlareDesign,0.9>=Ma[1:s]>=0, start=0.2)		#Mach number

@variable(FlareDesign,delx[1:s]>=0)				#Flame distortion measure, ft 

@variable(FlareDesign,dely[1:s]>=0)				#Flame distortion measure, ft

@variable(FlareDesign,hp[1:s]>=30, start=42)			#Height of flame, ft

@variable(FlareDesign,rp[1:s]>=0, start=24)			#Reference radius, ft

@variable(FlareDesign,6000>=K[1:s]>=0,start=100)		#Radiation, BTU/(h*ft^2)

@variable(FlareDesign,D[1:s]>=30, start=160)			#Distance from the reference point to flame, ft		

@variable(FlareDesign,cost>=0)					#Flare stack cost, USD

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

#Probability Constraint
@variable(FlareDesign, phi[1:s] >= 0, start=100)    	#CVaR auxiliary variable

@variable(FlareDesign, VaR>=0, start=100)            	#CVaR auxiliary variable

@variable(FlareDesign, 0<= CVaR <=threshold, start=100)

@constraint(FlareDesign, cvar[j in 1:s], K[j]-VaR <= phi[j])

@constraint(FlareDesign, CVaR == VaR + (1/s)*sum((1/alpha)*phi[j] for j=1:s))

#Objective function
@NLobjective(FlareDesign,Min, cost)
status=solve(FlareDesign)
#Results
println("Mean Cost ", getobjectivevalue(FlareDesign))
println("Diameter ft ",getvalue(dft))
println("Height ",getvalue(h))
lastobj = getobjectivevalue(FlareDesign)


Kvalue = getvalue(K)
Kvalue = sort(Kvalue)
index = round(Integer,(1-alpha)*s)
println("Var from quntile  ", Kvalue[index])
println(sum((Kvalue-threshold).<=0)/s)


SVar_obj = []
SVar_status = []
SVar_prob =[]
SVar_Var =[]
SVar_b =[]
SVar_c =[]
b= 2.506
a = -1/(getvalue(VaR)-threshold)
#a = -1/(VaRs - threshold)
c= a*(b+1)/2
for iter = 1:100
    dfts = getvalue(dft);
    hs = getvalue(h);
    ujs = getvalue(uj);
    Mas = getvalue(Ma);
    delxs = getvalue(delx);
    delys = getvalue(dely);
    hps   = getvalue(hp);
    rps   = getvalue(rp);
    Ks   = getvalue(K);
    Ds   = getvalue(D);
    costs   = getvalue(cost);

    FlareDesign=Model(solver=IpoptSolver(max_iter=500,tol=1e-6))

    #VARIABLES
    @variable(FlareDesign,60/12>=dft>=1/12, start =dfts)             #Diameter, ft

    @variable(FlareDesign,600>=h>=30, start=hs)                      #Height of flare stack, ft

    @variable(FlareDesign,umax>=uj[j in 1:s]>=0, start=ujs[j])       #Flare tip exit velocity, ft/s

    @variable(FlareDesign,0.9>=Ma[j in 1:s]>=0, start=Mas[j])        #Mach number

    @variable(FlareDesign,delx[j in 1:s]>=0, start=delxs[j])         #Flame distortion measure, ft

    @variable(FlareDesign,dely[j in 1:s]>=0, start=delys[j])         #Flame distortion measure, ft

    @variable(FlareDesign,hp[j in 1:s]>=30, start=hps[j])            #Height of flame, ft

    @variable(FlareDesign,rp[j in 1:s]>=0, start=rps[j])             #Reference radius, ft

    @variable(FlareDesign,6000>=K[j in 1:s]>=0,start=Ks[j])          #Radiation, BTU/(h*ft^2)

    @variable(FlareDesign,D[j in 1:s]>=30, start=Ds[j])	             #Distance from the reference point to flame, ft

    @variable(FlareDesign,cost>=0, start=costs)         #Flare stack cost

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

    #Objective function
    @NLobjective(FlareDesign,Min, cost)

    #Probability Constraint
    @variable(FlareDesign, phi[1:s] >= 0, start=100)    # svar auxiliary variable
    @variable(FlareDesign, z[1:s])           # svar auxiliary variable
    @constraint(FlareDesign, zdef[j=1:s], K[j] - threshold == z[j])
    @NLconstraint(FlareDesign, phidef[j=1:s],   2*(1+b)/(b+exp(-c*z[j])) - 1 <= phi[j])
   
    @constraint(FlareDesign, sum{phi[j], j in 1:s} <= alpha*s)
    for j in 1:s
        setvalue(z[j], Ks[j]-threshold);
        setvalue(phi[j], max(2*(1+b)/(b+exp(-c*(Ks[s]-threshold))),0) );
    end


    status=solve(FlareDesign)
    if status != :Optimal
        break
    end

	println("Mean Cost ", getobjectivevalue(FlareDesign))
	println("Diameter ft ",getvalue(dft))
	println("Height ",getvalue(h))

    Kvalue = getvalue(K)
    Kvalue = sort(Kvalue)
    index = round(Integer,(1-alpha)*s)
    println("Var from quntile  ", Kvalue[index])

    push!(SVar_status, status)
    push!(SVar_obj, getobjectivevalue(FlareDesign))
    push!(SVar_Var, Kvalue[index])
    push!(SVar_prob, sum((Kvalue-threshold).<=0)/s)
    push!(SVar_b, b)
    push!(SVar_c, c)
    b = 2*b
    c= a*(b+1)/2

    newobj = getobjectivevalue(FlareDesign)
    if abs(newobj-lastobj)<=0.001
        break
    end
    lastobj = newobj
end

#Results
println(SVar_status)
println(SVar_obj)
println(SVar_Var)
println(SVar_prob)
println(SVar_b)
println(SVar_c)

