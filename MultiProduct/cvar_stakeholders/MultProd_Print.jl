println("PRINTING RESULT FILES");
println("-------------------------------------------");

println("Printing Technology Siting");
open("technology_results.csv", "w") do fp
println(fp,"Node,t1,t2,t3,t4")

		for j in NODES
		       println(fp,j,",",getvalue(y["t1",j]),",",getvalue(y["t2", j]),",", getvalue(y["t3",j]),",", getvalue(y["t4",j]))
	        end
end

println("Printing Flow Matrices")
open("flow_p1_results.csv", "w") do f1
open("flow_p2_results.csv", "w") do f2
open("flow_p3_results.csv", "w") do f3
open("flow_p4_results.csv", "w") do f4

	for f in [f1,f2,f3,f4]

		if f == f1
			p = "p1"
		elseif f == f2
			p = "p2"
		elseif f == f3
			p = "p3"
		elseif f == f4
			p = "p4"
		end

		print(f,",")

		for j in NODES
			print(f,j,",")
		end
		println(f,"")
	
		for j in NODES
			print(f,j,",")
		for k in NODES
		       print(f,getvalue(flow[j,k,p]),",")
	        end

	        println(f)
	        end
	end

end
end
end
end

println("Printing Demand Values")
open("demand_results.csv", "w") do dp
println(dp,"Demand Number,Node,Product,Demand Value")
		for d in DEMS
		       println(dp,d,",",dem_node[d],",",dem_prod[d],",", getvalue(dem[d]),",")
	        end
end

println("Printing Summary File")

open("summary_results.csv", "w") do sp
println(sp,"Investment Cost,Transportation Cost,Social Welfare Function,Total Technologies Installed")
println(sp,getvalue(invcost),",",getValue(transcost),",",getValue(swf),",", getValue(sum(y)),",")
end




println("-------------------------------------------");
println("PRINTING COMPLETE");
