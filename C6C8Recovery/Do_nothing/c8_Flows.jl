println("Printing Flow Matrices")
open("flow_p1_results_"*"$(epsilon)"*".csv", "w") do f1
open("flow_p2_results_"*"$(epsilon)"*".csv", "w") do f2
open("flow_p3_results_"*"$(epsilon)"*".csv", "w") do f3
open("flow_p4_results_"*"$(epsilon)"*".csv", "w") do f4
open("flow_p5_results_"*"$(epsilon)"*".csv", "w") do f5
open("flow_p6_results_"*"$(epsilon)"*".csv", "w") do f6                    
open("flow_p7_results_"*"$(epsilon)"*".csv", "w") do f7
open("flow_p8_results_"*"$(epsilon)"*".csv", "w") do f8
open("flow_p9_results_"*"$(epsilon)"*".csv", "w") do f9
open("flow_p10_results_"*"$(epsilon)"*".csv", "w") do f10
open("flow_p11_results_"*"$(epsilon)"*".csv", "w") do f11
open("flow_p12_results_"*"$(epsilon)"*".csv", "w") do f12
open("flow_p13_results_"*"$(epsilon)"*".csv", "w") do f13
open("flow_p14_results_"*"$(epsilon)"*".csv", "w") do f14
open("flow_p15_results_"*"$(epsilon)"*".csv", "w") do f15
open("flow_p16_results_"*"$(epsilon)"*".csv", "w") do f16  
                                                                
	for ff in [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16]

		if ff == f1
			p = "p1"
		elseif ff == f2
			p = "p2"
		elseif ff == f3
			p = "p3"
		elseif ff == f4
			p = "p4"
		elseif ff == f5
			p = "p5"
        elseif ff == f6
			p = "p6" 
        elseif ff == f7
			p = "p7"
		elseif ff == f8
			p = "p8"
		elseif ff == f9
			p = "p9"
        elseif ff == f10
			p = "p10"
        elseif ff == f11
			p = "p11"                                                                
        elseif ff == f12
			p = "p12"
		elseif ff == f13
			p = "p13"
		elseif ff == f14
			p = "p14"
		elseif ff == f15
			p = "p15"
        elseif ff == f16  
			p = "p16"
		end
		#show(p)
		print(ff,",")

		for j in NODES			# Prints the header with node index
			print(ff,j,",")
		end
		println(ff,"")			# Used to enter next line
	
		for j in NODES			# Prints the 1st row entry i.e. the sender node
			print(ff,j,",")
		for k in NODES			# Prints the flow value from node j to node k with product p
		       print(ff,getvalue(f[j,k,p]),",")
	        end

	        println(ff)
	        end
	end

end
end
end
end
end
end
end
end
end
end
end
end
end
end
end
end                        