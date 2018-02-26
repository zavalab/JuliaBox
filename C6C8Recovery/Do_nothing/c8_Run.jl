#workspace()
#println("-------------------------------------------");
#println("Multi-Product Framework");
#rintln("Victor M. Zavala 2016 - UW Madison");
#println("Apoorva Sampat 2016 - UW Madison");
println("-------------------------------------------\n");

# Read the model
println("Reading model, ");
@time(include("c8_mix.jl"));

# SOLVE
println("-------------------------------------------");
println("SOLVER BEGINS");
println("-------------------------------------------");

#open("demand_results_w_invcap.csv", "a") do dp
open("summary_results_w_invcap.csv", "a") do rp # writes a summary of results, similar to the ones printed on console
#open("tradeoff.csv", "a") do tp # writes a file used in plotting the tradeoff curve


solve(m)

println("-------------------------------------------");
println("SOLVER ENDS");
println("-------------------------------------------");

#include("MultProd_Print.jl")

#include("Technologies_Sited.jl")  # writes a file listing the technologies sited
#include("MultProd_Flows.jl")  # Used to write files useful in plotiing flows across the map of Wisconsin



println(rp, "Investment Cost (Million \$)", ",", getvalue(invcost))
println(rp, "Transportation Cost (Million \$/yr)", ",",  getvalue(transcost))
println(rp, "Operating Cost (Million \$/yr)", ",", getvalue(opcost))
println(rp, "Total Annual Cost (Million \$/yr)", ",",  getvalue(transcost+opcost+invcost/20))
println(rp, "Total Revenue (Million \$/yr)", ",",  -getvalue(swf))
println(rp, "Total Annual Profit (Million \$/yr)", ",",  getvalue(profit))        
println(rp, "Total Technologies Installed", ",",  getvalue(sum(y)))
println(rp, "Total Biogas Technologies Installed", ",", getvalue(Bgtech))
println(rp, "Total C6/C8 Technologies Installed", ",", getvalue(C8tech))        
println(rp, "Total C8 Generated (Million t/year)", ",",  getvalue(totalprod["p4"]))
println(rp, "Total C6 Generated (Million t/year)", ",",  getvalue(totalprod["p5"]))        
println(rp, "Total Biogas Generated (Million t/year)", ",", getvalue(totalprod["p11"]))
println(rp, "Percentage of Unprocessed Manure %", ",", 100*getvalue(unamp))
println(rp, "Percentage of Unprocessed Sludge %", ",", 100*getvalue(unsp))
println(rp, "Percentage of Unprocessed Food Waste %", ",", 100*getvalue(unfwp))
println(rp, "CH4 Emitted by Unprocessed Waste (Million t/yr)",",",getvalue(totalch4)) 
println(rp, "CO2 Emitted by Transportation and Processing(Million t/yr)",",",getvalue(totalco2))
println(rp, "CO2 Saved by produing Biogas (Million t/yr)", ",", getvalue(SCO2))
println(rp, "Total Social Cost of Carbon (Million \$/yr)", ",", getvalue(TSCC))
println(rp, "Saved Social Cost of Carbon (Million \$/yr)", ",", getvalue(SSCC))                    
println(rp, "Net Social Cost of Carbon (Million \$/yr)",",",getvalue(NSCC))
println(rp, "Total Social Welfare (Million \$/yr)",",",getvalue(profit-NSCC))        
        

print(rp, "\n")


#print(tp, "\n")
#print(tp, epsilon)
#print(tp,",")
#print(tp, getvalue(invcost))
#print(tp,",")
#print(tp, getvalue(transcost))
#print(tp,",")
#print(tp, getvalue(opcost))
#print(tp,",")
#print(tp, getvalue(transcost+opcost+invcost/20))
#print(tp,",")
#print(tp, getvalue(totalprod["p1"]))
#print(tp,",")
#print(tp, getvalue(totalprod["p2"]))
#print(tp,",")
#print(tp, getvalue(totalprod["p3"]))
#print(tp,",")
#print(tp, getvalue(unamp))
#print(tp,",")
#print(tp, getvalue(unsp))
#print(tp,",")
#print(tp, getvalue(unfwp))
#print(tp,",")
#print(tp, getvalue(unamp+unsp+unfwp))
#print(tp,",")        
#print(tp, getvalue(profit))
#print(tp,",")
#print(tp, getvalue(umanure)/10181035.74*100)
#print(tp,",")


println("Investment Cost (Million \$) 					: ", getvalue(invcost))
println("Transportation Cost (Million \$/yr) 				: ", getvalue(transcost))
println("Operating Cost (Million \$/yr))  				: ", getvalue(opcost))
println("Total Annual Cost (Million \$/yr) 				: ", getvalue(transcost+opcost+invcost/20))
println("Total Revenue (Million \$/yr)  					: ", -getvalue(swf))
println("Total Annual Income (Million \$/yr)				: ", getvalue(profit))
println("Total Technologies Installed 				: ", getvalue(sum(y)))
println("Total Biogas Technologies Installed			: ", getvalue(Bgtech))
println("Total C6/C8 Technologies Installed			: ", getvalue(C8tech))
println("Total C8 Generated (Million t/year)				: ", getvalue(totalprod["p4"]))
println("Total C6 Generated (Million t/year)				: ", getvalue(totalprod["p5"]))        
println("Total Biogas Generated (Million t/year)				: ", getvalue(totalprod["p11"]))
println("Percentage of Unprocessed Manure %			: ", 100*getvalue(unamp))
println("Percentage of Unprocessed Sludge %			: ", 100*getvalue(unsp))
println("Percentage of Unprocessed Food Waste %			: ", 100*getvalue(unfwp))
println("CH4 Emitted by Unprocessed Waste (Million t/yr)			: ", getvalue(totalch4)) 
println("CO2 Emitted by Transportation and Processing(Million t/yr)	: ", getvalue(totalco2))
println("CO2 Saved by produing Biogas (Million t/yr)			: ", getvalue(SCO2))
println("Total Social Cost of Carbon (Million \$/yr)			: ", getvalue(TSCC))
println("Saved Social Cost of Carbon (Million \$/yr)			: ", getvalue(SSCC))                    
println("Net Social Cost of Carbon (Million \$/yr)			: ", getvalue(NSCC))
println("Total Social Welfare (Million \$/yr)				: " ,getvalue(profit-NSCC))        
            
        


#end
#end
end
