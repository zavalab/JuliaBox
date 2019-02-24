

open("technologies_sited_"*"$(epsilon)"*".csv", "w") do lp
	println(lp, "#Node Number", ",", "Latitude",",","Longitude",",", "Technology", ",", "Capacity")
	for n in NODES
		for t in TECHS
			if round(Int, getvalue(y[t, n])) == 1
			println(lp, n, ",", node_lat[n], ",", node_long[n],",",t, ",", tech_cap[t])
			end
		end
	end
end

