##Facility location problem
##Modeled in JuMP, Solver = Gurobi

println("-------------------------------------------");
println("READING MODEL");

using JuMP
using Gurobi

## Import Data
customer_matrix = readcsv("Customer.csv");           # Customers information
source_matrix = readcsv("Source.csv");               # Potential sources information
facility_matrix = readcsv("Facilities.csv");         # Data about facilities
demand_matrix = readcsv("Demand.csv");               # Demands in 3 scenarios
probability_matrix = readcsv("Probability.csv");     # Probability of 3 scenarios

## Define Sets
customer = customer_matrix[:,1];                     # Set of custormer nodes
source = source_matrix[:,1];                         # Set of source nodes
facility = facility_matrix[:,1];                     # Set of facilities we can intall
scenario = 1:length(probability_matrix[:,1]);        # Set of scenarios

## Define Normal Dictionaries 
Xc = Dict(zip(customer, customer_matrix[:,2]));      # Match each customer with its location(x,y)[=](km,km)
Yc = Dict(zip(customer, customer_matrix[:,3]));
Xs = Dict(zip(source, source_matrix[:,2]));          # Match each source node with its location(x,y)[=](km,km)
Ys = Dict(zip(source, source_matrix[:,3]));
Cap= Dict(zip(facility, facility_matrix[:,2]));      # Match each facility with its maximum production capacity(tonne/year)
Inv= Dict(zip(facility, facility_matrix[:,3]));      # Match each facility with its investment cost($)
Op = Dict(zip(facility, facility_matrix[:,4]));      # Match each facility with its operational cost($/year)
Pr = Dict(zip(scenario, probability_matrix[:,1]));   # Match each scenario with its probability

## Define Two-key Dictionaries
distance = Dict(("A", "I") => 60.21);                # Distances beween source and customer, the first line is an initiator
for i in source
    for j in customer
       distance[(i, j)] = sqrt((Xs[i]-Xc[j])^2 + (Ys[i]-Yc[j])^2);
    end
end

demand = Dict(("I",1) => 100.5);                     # Demands of customers in each scenario
for j in 1:length(customer)
    for s in scenario
        demand[(customer[j],s)] = demand_matrix[j,s];
    end
end

## Establish Model

m = Model(solver = GurobiSolver());

@variable(m, y[source,facility], Bin);               # Binary variables: technology t is installed in source i <=> y[i,t] = 1
@variable(m, f[source,customer,scenario]>=0);        # Flow of product from source i to customer j in scenario s
@variable(m, z[customer,scenario]>=0);               # Unsatisfied demand of customer j in scenario s

# At most one facility can be placed at one source
@constraint(m, onetech[i in source], sum(y[i,t] for t in facility) <=1);

# Balance constraint in each customer node
@constraint(m, balance[j in customer, s in scenario], sum(f[i,j,s] for i in source) + z[j,s] >= demand[j,s]);

# Capacity constraint in each source node
@constraint(m, capacity[i in source, s in scenario], sum(f[i,j,s] for j in customer) <= sum(y[i,t]*Cap[t] for t in facility));

# Define some intermediate quantities
@variable(m, investment>=0);                         # Total annual investment cost
@variable(m, operation>=0);                          # Total operational cost
@variable(m, transcost>=0);                          # Expected total transportation cost (summed and weighted by the probability)
@variable(m, penalty>=0);                            # Expected total penalty caused by unsatisfied demand

# Calculate total investement cost (annualized by life = 20 year)
@constraint(m, investment == sum(y[i,t]*Inv[t]/20 for i in source for t in facility));

# Calculate total operational cost
@constraint(m, operation == sum(y[i,t]*Op[t] for i in source for t in facility));
                     
# Calculate expected total transportation cost with transportation cost = $0.3/tonne/km                              
@constraint(m, transcost == sum(f[i,j,s]*distance[i,j]*0.3*Pr[s] for i in source for j in customer for s in scenario));     

# Calculate expected total penalty with $40 per unit of unsatisfied demand                                                  
@constraint(m, penalty == sum(Pr[s]*40*z[j,s] for j in customer for s in scenario));   
    
# Objective function is to minimize total annual expense                                                                    
@objective(m, Min, investment + operation + transcost + penalty); 
                                                                    
## Solve the model and output results
println("-------------------------------------------");
println("SOLVER BEGINS");
println("-------------------------------------------");
solve(m);
println("-------------------------------------------");
println("SOLVER ENDS");
println("-------------------------------------------\n"); 
println("Facility Installation Result:\n");
for i in source
for t in facility
    if getvalue(y[i, t]) >= 0.5
      println("	",t," is installed at ", i, "(",Xs[i],",",Ys[i],")");
    end
end
end
println("\n");
println("Expected Total Annual Cost(\$):\n	", getobjectivevalue(m));