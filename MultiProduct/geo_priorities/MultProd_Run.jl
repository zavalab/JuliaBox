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

#open("demand_results_w_invcap.csv", "a") do dp
open("summary_results_w_invcap.csv", "a") do rp
open("tradeoff.csv", "a") do tp

#for alpha in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]	# This loop does not work!
#Test mode#
solve(m)

println("-------------------------------------------");
println("SOLVER ENDS");
println("-------------------------------------------");

#include("MultProd_Print.jl")

include("Technologies_Sited.jl")
include("MultProd_Flows.jl")

#Writing results

#=
for nn in NODES
print(dp, getvalue(demtot[nn, "p1"]))
print(dp,",")
end
print(dp, "\n")
=#

println(rp, "Investment Cost (\$)", ",", getvalue(invcost))
println(rp, "Transportation Cost (\$/day)", ",",  getvalue(transcost))
println(rp, "Operating Cost (\$/yr)", ",", getvalue(opcost))
println(rp, "Total Daily Cost (\$/day)", ",",  getvalue(daily_cost))
println(rp, "Social Welfare Function (\$/day)", ",",  getvalue(swf))
println(rp, "Total Technologies Installed", ",",  getvalue(sum(y)))
println(rp, "Total Struvite Recovered (kg/day)", ",",  getvalue(struvite_total))
println(rp, "Total Unprocessed Manure (kg/day)", ",", getvalue(umanure))
println(rp, "Percentage of Manure Unprocessed (%)", ",", getvalue(umanure)/10181035.74*100)
println(rp, "Total Manure Hauled (kg/day)", ",", getvalue(tot_flow["p1"]))
println(rp, "Total Struvite Hauled (kg/day)", ",", getvalue(tot_flow["p2"]))
println(rp, "Total Digestate Hauled (kg/day)", ",", getvalue(tot_flow["p3"]))
println(rp, "Average Haul Distance for Manure (km/day)",",", getvalue(dist_flow["p1"])/getvalue(tot_flow["p1"]))
println(rp, "Average Haul Distance for Struvite (km/day)",",", getvalue(dist_flow["p2"])/getvalue(tot_flow["p2"]))
println(rp, "Average Haul Distance for Digestate (km/day)",",", getvalue(dist_flow["p3"])/getvalue(tot_flow["p3"]))
println(rp, "epsilon", ",",  epsilon)
#println(rp, "CVaR", ",",  getvalue(cvar))
#println(rp, "alpha", ",",  alpha)
print(rp, "\n")


print(tp, "\n")
print(tp, epsilon)
print(tp,",")
print(tp, getvalue(invcost))
print(tp,",")
print(tp, getvalue(transcost))
print(tp,",")
print(tp, getvalue(opcost))
print(tp,",")
print(tp, getvalue(daily_cost))
print(tp,",")
print(tp, getvalue(umanure))
print(tp,",")
print(tp, getvalue(umanure)/10181035.74*100)
print(tp,",")
print(tp, getvalue(sum(y)))
print(tp, ",")
print(tp, getvalue(struvite_total))
print(tp, ",")
print(tp, getvalue(dist_flow["p1"])/getvalue(tot_flow["p1"]))
print(tp, ",")
print(tp, getvalue(dist_flow["p2"])/getvalue(tot_flow["p2"]))
print(tp, ",")

println("Investment Cost (\$) 				: ", getvalue(invcost))
println("Transportation Cost (\$/day) 			: ", getvalue(transcost))
println("Operating Cost (\$/yr))  			: ", getvalue(opcost))
println("Total Daily Cost (\$/day) 			: ", getvalue(daily_cost))
println("Social Welfare Function (\$/day)  		: ", getvalue(swf))
println("Total Technologies Installed 			: ", getvalue(sum(y)))
println("Total Struvite Recovered (kg/day)		: ", getvalue(struvite_total))
println("Total Unprocessed Manure (kg/day)		: ", getvalue(umanure))
println("Percentage of Manure Unprocessed (%)		: ", getvalue(umanure)/10181035.74*100)
println("Total Manure Hauled (kg/day)			: ", getvalue(tot_flow["p1"]))
println("Total Struvite Hauled (kg/day)			: ", getvalue(tot_flow["p2"]))
println("Total Digestate Hauled (kg/day)		: ", getvalue(tot_flow["p3"]))
println("Average Haul Distance for Manure (km/day)	: ", getvalue(dist_flow["p1"])/getvalue(tot_flow["p1"]))
println("Average Haul Distance for Struvite (km/day)	: ", getvalue(dist_flow["p2"])/getvalue(tot_flow["p2"]))
println("Average Haul Distance for Digestate (km/day)	: ", getvalue(dist_flow["p3"])/getvalue(tot_flow["p3"]))
println("epsilon  					: ", epsilon)
#println("CVaR				: ", getvalue(cvar))
#println("alpha				: ", alpha)

#end
end
end
