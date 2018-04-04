println("-------------------------------------------\n");

# Read the model
println("Reading model, ");
@time(include("c8_mix.jl"));

# SOLVE
println("-------------------------------------------");
println("SOLVER BEGINS");
println("-------------------------------------------");


solve(m)

println("-------------------------------------------");
println("SOLVER ENDS");
println("-------------------------------------------");


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
