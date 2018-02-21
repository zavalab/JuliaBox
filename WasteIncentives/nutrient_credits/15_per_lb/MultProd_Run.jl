#workspace()
println("-------------------------------------------");
println("Multi-Product Framework");
println("Victor M. Zavala 2016 - UW Madison");
println("Apoorva Sampat 2016 - UW Madison");
println("-------------------------------------------\n");

#println("Reading data, ");
#@time(include("MultProd_Data_v2.jl"));

include("demand_writer.jl")

println("Reading model, ");
@time(include("MultProd_v2.jl"));


# SOLVE
println("-------------------------------------------");
println("SOLVER BEGINS");
println("-------------------------------------------");

#print("Including Phosphorus balances, ");
#@time(include("scaling.jl"));

#open("demand_results_w_invcap.csv", "a") do dp
open("./Summary_Output/summary_results_w_invcap.csv", "a") do rp
#open("./Summary_Output/tradeoff.csv", "a") do tp

solve(m)

println("-------------------------------------------");
println("SOLVER ENDS");
println("-------------------------------------------");

#include("MultProd_Print.jl")

include("Technologies_Sited.jl")
include("MultProd_Flows.jl")

roi =  (getvalue(prod_revenue) - (getvalue(invcost)/365/20) - (getvalue(opcost)/365) - getvalue(transcost) + getvalue(p_credit))*365/(getvalue(invcost))*100

#Writing results

println(rp, "Investment Cost (\$)", ",", thoudand_to_dollar*getvalue(invcost))
println(rp, "Transportation Cost (\$/day)", ",",  thoudand_to_dollar*getvalue(transcost))
println(rp, "Operating Cost (\$/yr)", ",", thoudand_to_dollar*getvalue(opcost))
println(rp, "Revenue from Products (\$/day)", ",",  thoudand_to_dollar*getvalue(prod_revenue))
println(rp, "Revenue from P Credit (\$/day)", ",",  thoudand_to_dollar*getvalue(p_credit))
println(rp, "Return on Investment (\%)", ",", roi)
println(rp, "Payback Period (yrs)", ",", 100/roi)
#println(rp, "Revenue from P credit (\$/day)", ",",  getvalue(p_credit_value))
println(rp, "Total Technologies Installed", ",",  getvalue(sum(y)))
println(rp, "Percentage of Manure Unprocessed (%)", ",", tons_to_kg*getvalue(umanure)/sup_tot*100)
# Printing total product recovered
for p in PRODS
       println(rp, "Total $(prod_alias[p]) recovered ($(prod_units[p])/day)", ",", tons_to_kg*getvalue(prod_recovered[p]))
end
# Printing Mass Distance for all products
for p in PRODS
	println(rp, "Mass-Distance for $(prod_alias[p]) ($(prod_units[p]) km/day)", ",", tons_to_kg*getvalue(dist_flow[p]))
end
println(rp, "Percentage of Manure Moved ", ",", tons_to_kg*getvalue(tot_flow[prod_num["Manure"]])/sup_tot*100 )
println(rp, "CO2 emissions	(kg/day)", ",", getvalue(co2_total))
println(rp, "Objective Value (\$/day)", ",",  thoudand_to_dollar*getvalue(obj))
println(rp, "epsilon", ",",  epsilon)
print(rp, "\n")

# Printing on terminal
println("Investment Cost (\$) 				: ", thoudand_to_dollar*getvalue(invcost))
println("Transportation Cost (\$/day) 			: ", thoudand_to_dollar*getvalue(transcost))
println("Operating Cost (\$/yr))  			: ", thoudand_to_dollar*getvalue(opcost))
println("Revenue from Products (\$/day)     		: ", thoudand_to_dollar*getvalue(prod_revenue))
println("Revenue from P Credits (\$/day)			: ", thoudand_to_dollar*getvalue(p_credit))
println("Return on Investment (\%)         		: ", roi)
println("Payback Period (yrs)               		: ", 100/roi)
#println("Revenue from P credit (\$/day)                : ", getvalue(p_credit_value))
println()
println("Total Technologies Installed 			: ", getvalue(sum(y)))
# Printing total product recovered
for p in PRODS
       println("Total $(prod_alias[p]) recovered ($(prod_units[p])/day)       :",  tons_to_kg*getvalue(prod_recovered[p]))
end
println("Percentage of Manure Unprocessed (%)		: ", tons_to_kg*getvalue(umanure)/sup_tot*100)
println()

# Printing Mass Distance for all products
for p in PRODS
	println("Mass-Distance for $(prod_alias[p]) ($(prod_units[p]) km/day)		:", getvalue(dist_flow[p]))
end
println()
println("Percentage of Manure Moved			: ", tons_to_kg*getvalue(tot_flow[prod_num["Manure"]])/sup_tot*100)
println("CO2 emissions	(g/day)				: ", getvalue(co2_total))
println("Objective Value (\$/day)			: ", thoudand_to_dollar*getvalue(obj))
println("epsilon  					: ", epsilon)

#end
#end
end
