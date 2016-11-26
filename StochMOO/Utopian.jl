# solution of utopian points 
# Yankai Cao, Siyu Chen, Luis Fuentes
# UW-Madison, 2016

push!(LOAD_PATH,pwd())
using Ipopt
using Plasmo
using JuMP

include("data.jl")  #load chp data
NS=188              # number of scenarios
NB =188             # number of block
SPB=NS/NB           # number of scenarios per block


include("chp_model_up.jl")  #scenario model building function
Cost = Array(Float64, NB*3)  
GHGE = Array(Float64, NB*3)
SW = Array(Float64, NB*3)
s= 1
for j in 1:NB
     for n=1:3
     	 # get scenario model
         m = get_scenario_model((j-1)*SPB+1:j*SPB,1,NS)
         if n==1
             @objective(m, Min,getvariable(m,:TCost))
         end
         if n==2
             @objective(m,Min,getvariable(m,:GHGT))
         end
         if n==3
             @objective(m, Min, getvariable(m,:SWT))
         end
         solve(m)
         Cost[s]=getvalue(getvariable(m,:TCost))
         GHGE[s]=getvalue(getvariable(m,:GHGT))
         SW[s]=getvalue(getvariable(m,:SWT))
	 s = s + 1
     end
end

#get utopia point
CostUP = minimum(Cost)
GHGEUP = minimum(GHGE)
SWUP = minimum(SW)

println("CostUP: ", CostUP)
println("GHGEUP: ", GHGEUP)
println("SWUP ",SWUP)

#write results to files
filename = string("MatlabFigures/FormulationC/Cost_UP.txt")
writedlm(filename, Cost)
filename = string("MatlabFigures/FormulationC/GHGE_UP.txt")
writedlm(filename, GHGE)
filename = string("MatlabFigures/FormulationC/SW_UP.txt")
writedlm(filename, SW)
filename = string("MatlabFigures/FormulationD/Cost_UP.txt")
writedlm(filename, Cost)
filename = string("MatlabFigures/FormulationD/GHGE_UP.txt")
writedlm(filename, GHGE)
filename = string("MatlabFigures/FormulationD/SW_UP.txt")
writedlm(filename, SW)