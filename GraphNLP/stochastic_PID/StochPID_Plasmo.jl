using Plasmo, JuMP, MadNLP, MadNLPGraph
     
# sets
NS= 5      # number of scenarios
N=100;     # number of timesteps
Tf=10;     # final time
h=Tf/N;    # time step
T=1:N;     # set of times
Tm=1:N-1;  # set of times minus one
S=1:NS     # set of scenarios

# set time vector
time=zeros(N);
for t=1:N
  time[t] = h*(t-1);
end

# Define parameters used within each scenario
K    = 1.0;            # gain
x0   = 0.0;            # starting point
Kd   = 0.5;            # disturbance gain
tau  = 1.0             # inverse time constant
d    = fill(-1.0,5);   # disturbance
xsp  = [-2.0, -1.5, -0.5, 0.5, 1.0] # set point

# define a function to create a graph for each scenario
function get_scenario_model(s)
   graph = OptiGraph() # create OptiGraph
 
   @optinode(graph,n[T]) # add N nodes, n[1:100], to the optigraph
 
   # add decision variables to each node
   for node in n
 
     # define variables on node n
     @variable(node,-2.5<= x <=2.5)
     @variable(node,-2.0<= u <=2.0)
     @variable(node, int )
     @variable(node, Kc)
     @variable(node, tauI)
     @variable(node, tauD)
 
     # define dummy variables to avoid nonlinearity in @linkconstraint
     @variable(node, Kcx) 
     @variable(node, tauIint)
     @variable(node, tauDx)
     
     # define dummy variable constraints
     @constraint(node, Kcx == Kc * x)
     @constraint(node, tauIint == tauI * int)
     @constraint(node, tauDx  == tauD * x) 
 
     # objective function
     @objective(node, Min, (1/(NS))*(100*(xsp[s]-x)^2 + 0.01*u^2))
   end
 
   # set initial conditions
   @constraint(n[1], eqinix,  n[1][:int] == 0)
   @constraint(n[1], eqinit,  n[1][:x] == x0)
 
   # constraints
   @linkconstraint(graph, [t=Tm], (1/tau)*(n[t+1][:x]-n[t][:x])/h + n[t+1][:x]
                                                                              == K*n[t+1][:u] + Kd*d[s]);
   @linkconstraint(graph, [t=Tm], n[t+1][:u] == n[t][:Kc]*xsp[s] - n[t][:Kcx] 
                                              + n[t+1][:tauIint] + n[t+1][:tauDx]/h - n[t][:tauDx]/h);
   @linkconstraint(graph, [t=Tm], (n[t+1][:int]-n[t][:int])/h == (xsp[s]-n[t+1][:x]));
 
   #link stage one decision variables
   @linkconstraint(graph, [t=Tm], n[t][:Kc] == n[t+1][:Kc])
   @linkconstraint(graph, [t=Tm], n[t][:tauI] == n[t+1][:tauI])
   @linkconstraint(graph, [t=Tm], n[t][:tauD] == n[t+1][:tauD])
 
   return graph, n
 end
 PID=OptiGraph() # create initial high-level OptiGraph
 
 @optinode(PID,master) # create node called "master" for first stage variables
 
 # define tuning parameters (first stage variables)
 @variable(master, -10<= Kc <=10)
 @variable(master,-100<= tauI <=100)
 @variable(master,-100<= tauD <=100)
 
 for s in 1:NS           
   graph,nodes = get_scenario_model(s) # get OptiGraph and OptiNode for each Scenario
   add_subgraph!(PID,graph) # add scenario subgraph to high-level OptiGraph |\label{line:add_subgraphs}|
 
   # link constraints from subgraph to master node
   @linkconstraint(PID, nodes[1][:Kc]==master[:Kc])
   @linkconstraint(PID, nodes[1][:tauI]==master[:tauI])
   @linkconstraint(PID, nodes[1][:tauD]==master[:tauD])    
 end


# create reference map for partition
hypergraph, refmap = gethypergraph(PID)

# define a node membership vector that assigns each node index to a partition
node_membership_vector = Array{Int64,1}(undef,length(all_nodes(PID))) 

# fill node membership vector; value of vector entry corresponds to partition number
node_membership_vector[[1:26; 102:126; 202:226; 302:326; 402:426]]   .= 1 
node_membership_vector[[27:51; 127:151; 227:251; 327:351; 427:451]]  .= 2
node_membership_vector[[52:76; 152:176; 252:276; 352:376; 452:476]]  .= 3
node_membership_vector[[77:101; 177:201; 277:301; 377:401; 477:501]] .= 4 

# create partition
PID_partition = Partition(node_membership_vector,refmap) 

# repartition subgraphs according to the partition
make_subgraphs!(PID, PID_partition) 

MadNLP.optimize!(PID; linear_solver=MadNLPUmfpack)