#workspace()
println("-------------------------------------------");
println("Multi-Product Framework");
println("Victor M. Zavala 2016 - UW Madison");
println("Apoorva Sampat 2016 - UW Madison");
println("-------------------------------------------\n");

#println("Reading data, ");
#@time(include("MultProd_Data_v2.jl"));



println("Reading model, ");
@time(include("MultProd_v2.jl"));



# SOLVE
println("-------------------------------------------");
println("SOLVER BEGINS");
println("-------------------------------------------");

#print("Including Phosphorus balances, ");
#@time(include("scaling.jl"));

open("dissatisfaction_results_"*"$(file_str2)"*".csv", "a") do dp
open("summary_results_"*"$(file_str2)"*".csv", "a") do rp
#open("tradeoff.csv", "a") do tp

#for alpha in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]	# This loop does not work!
#Test mode#
solve(m)

println("-------------------------------------------");
println("SOLVER ENDS");
println("-------------------------------------------");

include("MultProd_Plot.jl") #Plots dissatisfaction histograms

include("Technologies_Sited.jl")
include("MultProd_Flows.jl")

#Writing results


print(dp, alpha)
print(dp,",")
for sh in SHOLDERS
print(dp, getvalue(dissatisfaction[sh]))
print(dp,",")
end
print(dp, "\n")


println(rp, "Investment Cost (\$)", ",", getvalue(invcost))
println(rp, "Transportation Cost (\$/day)", ",",  getvalue(transcost))
println(rp, "Operating Cost (\$/yr)", ",", getvalue(opcost))
println(rp, "Total Daily Cost (\$/day)", ",",  getvalue(daily_cost))
println(rp, "Social Welfare Function (\$/day)", ",",  getvalue(swf))
println(rp, "Total Technologies Installed", ",",  getvalue(sum(y)))
println(rp, "Total Struvite Recovered (kg/day)", ",",  getvalue(struvite_total))
println(rp, "Total Biogas Recovered (kg/day)", ",",  getvalue(biogas_total))
println(rp, "Total Unprocessed Manure (kg/day)", ",", getvalue(umanure))
println(rp, "Percentage of Manure Unprocessed (%)", ",", getvalue(umanure)/10181035.74*100)
println(rp, "epsilon", ",",  epsilon)
println(rp, "Total Manure Hauled (kg/day)", ",", getvalue(tot_flow["p1"]))
println(rp, "Total Struvite Hauled (kg/day)", ",", getvalue(tot_flow["p2"]))
println(rp, "Total Digestate Hauled (kg/day)", ",", getvalue(tot_flow["p3"]))
println(rp, "Average Haul Distance for Manure (km/day)",",", getvalue(dist_flow["p1"])/getvalue(tot_flow["p1"]))
println(rp, "Average Haul Distance for Struvite (km/day)",",", getvalue(dist_flow["p2"])/getvalue(tot_flow["p2"]))
println(rp, "Average Haul Distance for Digestate (km/day)",",", getvalue(dist_flow["p3"])/getvalue(tot_flow["p3"]))
println(rp, "CVaR", ",",  getvalue(cvar))
println(rp, "alpha", ",",  alpha)
println(rp, "Objective",",",file_str)
print(rp, "\n")

#=
print(tp, "\n")
print(tp, epsilon)
print(tp,",")
print(tp, getvalue(invcost))
print(tp,",")
print(tp, getvalue(transcost))
print(tp,",")
print(tp, getvalue(opcost))
print(tp,",")
#print(tp, getvalue(cost_total))
print(tp,",")
print(tp, getvalue(umanure))
print(tp,",")
print(tp, getvalue(umanure)/10181035.74*100)
print(tp,",")
=#

println("Investment Cost (\$) 				: ", getvalue(invcost))
println("Transportation Cost (\$/day) 			: ", getvalue(transcost))
println("Operating Cost (\$/yr))  			: ", getvalue(opcost))
println("Total Daily Cost (\$/day) 			: ", getvalue(daily_cost))
println("Social Welfare Function (\$/day)  		: ", getvalue(swf))
println("Total Technologies Installed 			: ", getvalue(sum(y)))
println("Total Struvite Recovered (kg/day)		: ", getvalue(struvite_total))
println("Total Biogas Recovered (kg/day)		: ", getvalue(biogas_total))
println("Total Unprocessed Manure (kg/day)		: ", getvalue(umanure))
println("Percentage of Manure Unprocessed (%)		: ", getvalue(umanure)/10181035.74*100)
println("epsilon  					: ", epsilon)
println("obj1_scaled  					: ", getvalue(obj1_scaled))
println("obj2_scaled  					: ", getvalue(obj2_scaled))
println("Total Manure Hauled (kg/day)			: ", getvalue(tot_flow["p1"]))
println("Total Struvite Hauled (kg/day)			: ", getvalue(tot_flow["p2"]))
println("Total Digestate Hauled (kg/day)		: ", getvalue(tot_flow["p3"]))
println("Average Haul Distance for Manure (km/day)	: ", getvalue(dist_flow["p1"])/getvalue(tot_flow["p1"]))
println("Average Haul Distance for Struvite (km/day)	: ", getvalue(dist_flow["p2"])/getvalue(tot_flow["p2"]))
println("Average Haul Distance for Digestate (km/day)	: ", getvalue(dist_flow["p3"])/getvalue(tot_flow["p3"]))
println("CVaR				: ", getvalue(cvar))
println("alpha				: ", alpha)
println("Objective					: ", file_str)

#end
end
end
