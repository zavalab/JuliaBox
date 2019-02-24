# solution of windTurbine problem
# subject to load constraint
# compute optimal setpoint for theta and Tgen
# pitch and torque controller is not included
# Yankai Cao, Victor Zavala
# UW-Madison, 2016

using Ipopt
using Plasmo
using JuMP
using Distributions
include("simulator.jl")
include("createModel_svar.jl")            # scenario model building function
using MAT
using JLD


include("setup.jl")



# get initial value created by initial.jl
d = load("Result/initialCVar.jld")
Vm = d["Vm"]
wr_initial = d["wr_initial"]
theta_initial = d["theta_initial"]
Tgen_initial = d["Tgen_initial"]
Power_initial = d["Power_initial"]
Fz_initial = d["Fz_initial"]
Mz_initial = d["Mz_initial"]
xfa_initial = d["xfa_initial"]
vel_xfa_initial = d["vel_xfa_initial"]
MyTB_initial = d["MyTB_initial"]
lambda_eff_initial = d["lambda_eff_initial"]
Ct_initial = d["Ct_initial"]
Cm_initial = d["Cm_initial"]
vel_xfa_dot_initial = d["vel_xfa_dot_initial"]
local_max_initial = d["local_max_initial"]
VaRs = d["Var_value"]
costvalue = local_max_initial - threshold



a_ = 1.0/(threshold-VaRs)
b= 2.506
c= (b+1)/2*a_
s = 1
if(rank == 0)
    temp = 0
    for i in Vbin
        for j in seedSet
            MyTB = Array{Float64}(Nt)
            MyTB[1:end] = MyTB_initial[s,:]'/1e7
            temp += max(0, 2*(1+b)/(b+exp(-c*(maximum(MyTB) - threshold))) - 1)
            s = s + 1
        end
    end
end
local_max_list = copy(local_max_initial)
sortedZ = sort(local_max_list)
index = round(Integer,(1-alpha)*S)



SVar_obj = []
SVar_prob =[]
SVar_Var =[]
SVar_x =[]
SVar_b =[]
SVar_c =[]

for iter = 1:100
println(b, "    ", c,"    ", a_)

graph = Plasmo.PlasmoGraph()
graph.solver = Ipopt.IpoptSolver()

m = Model()
master_node = add_node(graph,m)

@variable(m,  Tgen[t in firstTIMEG],  start = mean(Tgen_initial[:,t])/1e4)                       #Rated power [1e4 W]
@variable(m,  theta[t in firstTIMEG], start = mean(theta_initial[:,t]))                          #Reference generator speed in pitch controller [rad/s]
bl=Array{JuMP.Model}(S)
scen_nodes = Array{NodeOrEdge}(S)
owned = []
s = 1
for i in Vbin
    for j in seedSet
        push!(owned, s)
        MyTB = Array{Float64}(Nt)
        MyTB[1:end] = MyTB_initial[s,:]'/1e7

        # get scenario model
        bl[s] = createModel(s:s)
	node = add_node(graph,bl[s])
        scen_nodes[s] = node

        # add constraints couple first-stage and second-stage variables

	@linkconstraint(graph,[ t in firstTIMEG],   bl[s][:Tgen][s,t]==Tgen[t])
        @linkconstraint(graph,[ t in firstTIMEG],   bl[s][:theta][s,t]==theta[t])

        # add second-stage variables and constraints
        @variable(bl[s], local_max<=20, start =maximum(MyTB))
        @constraint(bl[s], local_max_def[t in TIMEG], local_max >=  bl[s][:MyTB][s,t])
       	@variable(bl[s], -threshold<=z<=20-threshold, start = maximum(MyTB) - threshold)
        @variable(bl[s], phi>=0, start =  max(0, 2*(1+b)/(b+exp(-c*(maximum(MyTB) - threshold))) - 1)     )    # cvar auxiliary variable
        @NLconstraint(bl[s], z  == local_max - threshold)
	@NLconstraint(bl[s], phi   >= 2*(1+b)/(b+exp(-c*z)) - 1)

        # add second stage objective
        @objective(bl[s], Min, -100*sum(deltaP[round(Int, ceil(s/length(seedSet)))]*bl[s][:Power][s,t] for t = TIMEG)/Nt/length(seedSet))
    	s = s + 1
    end
end	
@linkconstraint(graph, sum(bl[s][:phi] for s in owned)<=alpha*S)


status = Plasmo.solve(graph)

if status != :Optimal
    break
end

#print and store data if coreid is 0
    s=1
    Power = 0
    for i in 1:length(Vbin)
        Pmean = 0
        for j in seedSet
            P = bl[s].colVal[(getvariable(bl[1], :Power)[1,1].col):(getvariable(bl[1], :Power)[1,Nt].col)]
            Pmean += mean(P)
            s= s + 1
        end
        Pmean = Pmean/length(seedSet)
        Power += Pmean*deltaP[i]
    end
    println("Optimization average Power:  ", Power/1e2, " MW")
    if iter >= 2  && Power/1e2 <= SVar_obj[end]
       break
    end   
    push!(SVar_obj, Power/1e2)


    MyTB_list = Array{Float64}(length(SCEN), Nt)
    wr_list = Array{Float64}(length(SCEN), Nt)
    Power_list = Array{Float64}(length(SCEN), Nt)
    theta_list = Array{Float64}(length(SCEN), Nt)
    Tgen_list = Array{Float64}(length(SCEN), Nt)
    local_max_list = Array{Float64}(length(SCEN))

    Fz_list = Array{Float64}(S, Nt)
    Mz_list = Array{Float64}(S, Nt)
    xfa_list = Array{Float64}(S, Nt, ncp)
    vel_xfa_list = Array{Float64}(S, Nt, ncp)
    lambda_eff_list = Array{Float64}(S, Nt)
    Ct_list = Array{Float64}(S, Nt)
    Cm_list = Array{Float64}(S, Nt)
    vel_xfa_dot_list = Array{Float64}(S, Nt, ncp)

    s=1
    for i in 1:length(Vbin)
        for j in seedSet
            MyTB = zeros(Nt)
            MyTB[1:end] = bl[s].colVal[(getvariable(bl[1], :MyTB)[1,1].col):(getvariable(bl[1], :MyTB)[1,Nt].col)]
            MyTB_list[s, :] = 1e7*MyTB
            local_max_value =  bl[s].colVal[(getvariable(bl[1], :local_max).col)]
            local_max_list[s] = maximum(MyTB)
  

            wr = zeros(Nt)
            wr[1:end] = bl[s].colVal[(getvariable(bl[1], :wr)[1,1].col):(getvariable(bl[1], :wr)[1,Nt].col)]
            wr_list[s, :] = wr
            Power = zeros(Nt)
            Power[1:end] = 1e4*bl[s].colVal[(getvariable(bl[1], :Power)[1,1].col):(getvariable(bl[1], :Power)[1,Nt].col)]
            Power_list[s, :] = Power
            theta = zeros(Nt)
            theta[1:end] = bl[s].colVal[(getvariable(bl[1], :theta)[1,1].col):(getvariable(bl[1], :theta)[1,Nt].col)]
            theta_list[s, :] = theta
            Tgen = zeros(Nt)
            Tgen[1:end] = 1e4*bl[s].colVal[(getvariable(bl[1], :Tgen)[1,1].col):(getvariable(bl[1], :Tgen)[1,Nt].col)]
            Tgen_list[s, :] = Tgen

            Fz_list[s,2:end] = 1e5*copy(bl[s].colVal[(getvariable(bl[1], :Fz)[1,2].col):(getvariable(bl[1], :Fz)[1,Nt].col)])'
            Mz_list[s,:] = 1e5*copy(bl[s].colVal[(getvariable(bl[1], :Mz)[1,1].col):(getvariable(bl[1], :Mz)[1,Nt].col)])'
            lambda_eff_list[s,:] = copy(bl[s].colVal[(getvariable(bl[1], :lambda_eff)[1,1].col):(getvariable(bl[1], :lambda_eff)[1,Nt].col)])'
            Ct_list[s,:] = copy(bl[s].colVal[(getvariable(bl[1], :Ct)[1,1].col):(getvariable(bl[1], :Ct)[1,Nt].col)])'
            Cm_list[s,:] = copy(bl[s].colVal[(getvariable(bl[1], :Cm)[1,1].col):(getvariable(bl[1], :Cm)[1,Nt].col)])'
            for t = 1:Nt
                for c in 1:ncp
                    xfa_list[s,t,c] = bl[s].colVal[(bl[1][:xfa][1,t,c].col)]
                    vel_xfa_list[s,t,c] = bl[s].colVal[(bl[1][:vel_xfa][1,t,c].col)]
                end
            end
            for t = 2:Nt
                for c in 1:ncp
                    vel_xfa_dot_list[s,t,c] = bl[s].colVal[(bl[1][:vel_xfa_dot][1,t,c].col)]
                end
            end
            s = s + 1
        end
    end
    sortedZ = sort(local_max_list)
    index = round(Integer,(1-alpha)*S)
    println("Var from quntile  ", sortedZ[index], "     ",index)
    push!(SVar_Var, sortedZ[index])
    costvalue2 = sortedZ-threshold
    push!(SVar_prob, sum(costvalue2.<=0)/S)
    push!(SVar_b, b)
    push!(SVar_c, c)
    save("Result/initialSVar.jld", "Vm", Vm, "wr_initial", wr_list, "theta_initial", theta_list, "Tgen_initial", Tgen_list, "Power_initial", Power_list, "Fz_initial", Fz_list, "Mz_initial", Mz_list, "xfa_initial", xfa_list, "vel_xfa_initial", vel_xfa_list, "MyTB_initial", MyTB_list, "lambda_eff_initial", lambda_eff_list, "Ct_initial", Ct_list, "Cm_initial", Cm_list, "vel_xfa_dot_initial",vel_xfa_dot_list, "local_max_initial",local_max_list)



d = load("Result/initialSVar.jld")
wr_initial = d["wr_initial"]
theta_initial = d["theta_initial"]
Tgen_initial = d["Tgen_initial"]
Power_initial = d["Power_initial"]
Fz_initial = d["Fz_initial"]
Mz_initial = d["Mz_initial"]
xfa_initial = d["xfa_initial"]
vel_xfa_initial = d["vel_xfa_initial"]
MyTB_initial = d["MyTB_initial"]
lambda_eff_initial = d["lambda_eff_initial"]
Ct_initial = d["Ct_initial"]
Cm_initial = d["Cm_initial"]
vel_xfa_dot_initial = d["vel_xfa_dot_initial"]
local_max_initial = d["local_max_initial"]

b = b*2
c= a_*(b+1)/2
end


println(SVar_obj)
println(SVar_Var)
println(SVar_prob)
println(SVar_b)
println(SVar_c)
println(SVar_x)
