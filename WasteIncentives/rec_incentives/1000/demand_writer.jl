# Apoorva Sampat
# July, 2017

# Read the input price values
prod_price_matrix = readdlm("./InputData/product_prices.csv", ',')
product_list = prod_price_matrix[:, 1]
price_list = prod_price_matrix[:, 2]

# Demand nodes
node_matrix = 	readdlm("./InputData/node_matrix.csv"	,',')
NODES 	    = 	node_matrix[:,1]		# set of nodes

# Writing the file
open("./InputData/demand_matrix.csv", "w") do fp
	println(fp,"# DEMS", ",", "dem_node", ",", "dem_prod", ",", "dem_cap", ",", "dem_price")
	
	d_count = 1 # Counter for demand number
	p_count  = 1 # Counter for product 

	for p in product_list
	 price =  price_list[p_count]
	
		for j in NODES
			 println(fp,"d",d_count,",",j, ",",p, ",", 60000, ",", price)
			 d_count = d_count + 1
		end

	p_count = p_count + 1
	end


end
