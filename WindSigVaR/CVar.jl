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
include("createModel.jl")         # scenario model building function
using MAT
using JLD


include("setup.jl")


# get initial value created by initial.jl
d = load("Result/initial.jld")
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



graph = Plasmo.PlasmoGraph()
graph.solver = Ipopt.IpoptSolver()

m = Model()
master_node = add_node(graph,m)

@variable(m,  Tgen[t in firstTIMEG],  start = mean(Tgen_initial[:,t])/1e4)                       #Rated power [1e4 W]
@variable(m,  theta[t in firstTIMEG], start = mean(theta_initial[:,t]))                          #Reference generator speed in pitch controller [rad/s]
@defVar(m, VaR)            # cvar auxiliary variable
@defVar(m, CVaR<=threshold)
bl=Array(JuMP.Model, S)
s=1
scen_nodes = Array{NodeOrEdge}(S)
owned = []
for i in Vbin
    for j in seedSet
        push!(owned, s)
        MyTB = Array(Float64, Nt)
        MyTB[1:end] = MyTB_initial[s,:]'

        # get scenario model
        bl[s] = createModel(s:s)
	node = add_node(graph,bl[s])
        scen_nodes[s] = node

        # add constraints couple first-stage and second-stage variables
	@linkconstraint(graph,[ t in firstTIMEG],   bl[s][:Tgen][s,t]==Tgen[t])
        @linkconstraint(graph,[ t in firstTIMEG],   bl[s][:theta][s,t]==theta[t])

	@variable(bl[s], local_VaR)
	@linkconstraint(graph,  local_VaR == VaR)


        # add second-stage variables and constraints
        @variable(bl[s], local_max<=20, start =maximum(MyTB))
        @constraint(bl[s], local_max_def[t in TIMEG], local_max >=  bl[s][:MyTB][s,t])
        #@constraint(bl[s], local_max_def[t in TIMEG], local_max >=  -bl[s][:MyTB][s,t])

	@variable(bl[s], phi >= 0, start = 10)    # cvar auxiliary variable
	@addConstraint(bl[s], cvar, local_max -local_VaR <= phi)	

        # add second stage objective
        @objective(bl[s], Min, -100*sum(deltaP[round(Int, ceil(s/length(seedSet)))]*bl[s][:Power][s,t] for t = TIMEG)/Nt/length(seedSet))
	#bl[s].solver=IpoptSolver(linear_solver="ma57", constr_mult_init_max=0.0, nlp_scaling_method="none")	    
    	s = s + 1
    end
end

@linkconstraint(graph, m[:CVaR] == m[:VaR] + (1/alpha)*(1/length(seedSet))*sum(deltaP[round(Int, ceil(s/length(seedSet)))]*bl[s][:phi] for s in owned))


Plasmo.solve(graph)

    s=1
    Power = 0
    println("VaR:  ", getvalue(getvariable(m, :VaR)))	   
    println("CVaR:  ", getvalue(getvariable(m, :CVaR)))
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

    MyTB_list = Array(Float64, length(SCEN), Nt)
    wr_list = Array(Float64, length(SCEN), Nt)
    Power_list = Array(Float64, length(SCEN), Nt)
    theta_list = Array(Float64, length(SCEN), Nt)
    Tgen_list = Array(Float64, length(SCEN), Nt)
    local_max_list = Array(Float64, length(SCEN))


    Fz_list = Array(Float64, S, Nt)
    Mz_list = Array(Float64, S, Nt)
    xfa_list = Array(Float64, S, Nt, ncp)
    vel_xfa_list = Array(Float64, S, Nt, ncp)
    lambda_eff_list = Array(Float64, S, Nt)
    Ct_list = Array(Float64, S, Nt)
    Cm_list = Array(Float64, S, Nt)
    vel_xfa_dot_list = Array(Float64, S, Nt, ncp)


    s=1
    for i in 1:length(Vbin)
        for j in seedSet
            MyTB = zeros(Nt)
            MyTB[1:end] = bl[s].colVal[(getvariable(bl[1], :MyTB)[1,1].col):(getvariable(bl[1], :MyTB)[1,Nt].col)]
            MyTB_list[s, :] = MyTB
	    local_max_value =  1e7*bl[s].colVal[(getvariable(bl[1], :local_max).col)]
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
	    #xfa_list[s,:] = copy(bl[s].colVal[(getvariable(bl[1], :xfa)[1,1].col):(getvariable(bl[1], :xfa)[1,Nt].col)])'
	    #vel_xfa_list[s,:] = copy(bl[s].colVal[(getvariable(bl[1], :vel_xfa)[1,1].col):(getvariable(bl[1], :vel_xfa)[1,Nt].col)])'
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
    VaRs = sortedZ[index]

    save("Result/initialCVar.jld", "Vm", Vm, "wr_initial", wr_list, "theta_initial", theta_list, "Tgen_initial", Tgen_list, "Power_initial", Power_list, "Fz_initial", Fz_list, "Mz_initial", Mz_list, "xfa_initial", xfa_list, "vel_xfa_initial", vel_xfa_list, "MyTB_initial", MyTB_list, "lambda_eff_initial", lambda_eff_list, "Ct_initial", Ct_list, "Cm_initial", Cm_list, "vel_xfa_dot_initial",vel_xfa_dot_list, "local_max_initial",local_max_list, "Var_value", VaRs)

