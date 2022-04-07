using Plasmo
using JLD2
using MadNLP
using MadNLPGraph
using MadNLPHSL
using Statistics
using DelimitedFiles
using JuMP

nt = 24 # number of time points
nx = 10 # number of space points per pipeline
horizon = 24*3600 # time horizon in seconds
dt = horizon / (nt - 1) # time delta

# Load in the necessary files
include("modelfunctions_stoch.jl")
include("load_data_stoch.jl")

# Define a function to time the solvers; Function 'timed_gasnetwork' receives the number of scenarios to run with a solver_num that indicates the solver to use
function timed_gasnetwork(scenario_num,solver_num)
    #Model Data
    global NS
    NS = scenario_num       #Number of scenarios       
    S = 1:NS



    # Make overall optigraph
    gas_network_stoch = OptiGraph() 

    # Define master node and master variables
    @optinode(gas_network_stoch,master)
    @variable(master,power[1:11,1:24])

    # Create array for subgraphs; each subgraph will contain its own scenario
    gas_networks = Array{Any,1}(undef,NS)

    # Build subgraph for each scenario and add the subgraph to the overall OptiGraph
    for l in 1:NS
        gas_networks[l] = make_stoch_subgraph(l)
        add_subgraph!(gas_network_stoch,gas_networks[l])
    end

    stoch_subs = getsubgraphs(gas_network_stoch)
    subs = Array{Any,1}(undef,NS)

    for j in S
        subs[j] = getsubgraphs(stoch_subs[j])
    end

    # Link the master variables to each scenario's compressor power variables
    for i in 1:11
        for j in 1:24
            for k in S
                @linkconstraint(gas_network_stoch, subs[k][i+13][:time_nodes][j][:power] == master[:power][i,j])
            end
        end
    end

    # Aggregate the subgraphs into nodes
    global gas_node
    gas_node,dict = aggregate(gas_network_stoch,0)

    
    global xip
    xip = zeros(3)

    # Solve the problem 4 times; the first time will be omitted to avoid compilation time
    for i in 1:4
        println("Running iteration $i")
        # Output files contained the solver (p or s indicates parallel or serial for Ma57) along with the number of scenarios, 
        # number of spacial discretization points, and which iteration the solution corresponds to (value from 1-4)

        # solver_num == 1: Run MadNLP in parallel with Ma57; this code only changes the output file name rather than the number of threads
        # to run in parallel, the user must set the number of threads that Julia is accessing; run Threads.nthreads() to identify number of threads
        if solver_num == 1
            x = @elapsed MadNLP.optimize!(gas_node;linear_solver=MadNLPMa57,max_iter=300,output_file="output_files/MadNLP_Ma57pfinal_$NS $nx $i")

        # solver_num == 2: Run MadNLP in serial with Ma57; see comment above 
        elseif solver_num == 2
            x = @elapsed MadNLP.optimize!(gas_node;linear_solver=MadNLPMa57,max_iter=300,output_file="output_files/MadNLP_Ma57sfinal_$NS $nx $i")
            
        # solver_num == 3: Run MadNLP with Schur solver
        elseif solver_num == 3
            print("Running with ", Threads.nthreads(), " threads")
            x = @elapsed MadNLP.optimize!(gas_node; linear_solver=MadNLPSchur, max_iter=300, schur_custom_partition=true, schur_subproblem_solver=MadNLPMa57, output_file="output_files/MadNLP_Schurfinal_$NS $nx $i")

        # solver_num == 4: Run MadNLP with PardisoMKL
        elseif solver_num == 4
            println("Running with ", Threads.nthreads(), " threads")
            x = @elapsed MadNLP.optimize!(gas_node, linear_solver=MadNLPPardisoMKL, max_iter=300, pardisomkl_num_threads = Threads.nthreads(), output_file="output_files/MadNLP_PardisoMKL_$NS $nx $i")
        end

        # If i != 1, save the elapsed time
        if i != 1
            xip[i-1] = x
	    end
    end

    # Save the number of iterations used
    xiter   = gas_node.optimizer.cnt.k
    # Save the objective value
    obj     = gas_node.optimizer.obj_val    
    # Save the average solution time
    xip_avg = mean(xip)
    # Return all solution times, average solution times, the number of iterations, and the objective value
    return xip, xip_avg, xiter, obj
end

# Define a function to run a given solver for a set of scenarios
function run_timed_gas(solver_num)
    # Define the number of scenarios to use
    scenario_array = [1,5,10,25,75,100,150]
    println("Running with solver number ", solver_num, "!")
    len = length(scenario_array)
    # Define an array to take each runs solution times
    # xips will contain the solver time for runs 2-4 for each set of scenarios
    xips = zeros((3,len))
    # Define an array to take additional information from the timed_gasnetwork function
    # xs will contain the average solver time, the number of iterations, the objective value,
    # the average time per iteration, the number of scenarios, and the solver number indicator
    xs = Array{Any,2}(undef,(6, len))

    # Iterate over all scenarios for a given solver. 
    for (i,j) in enumerate(scenario_array)
        println("Running iteration $i out of ", length(scenario_array), " for   ", j, "  scenarios")
        xips[:,i], xs[1,i],xs[2,i],xs[3,i] = timed_gasnetwork(j,solver_num)
        xs[4,i] = xs[1,i]/xs[2,i]
        xs[5,i] = j
        xs[6,i] = solver_num
    end

    return xips, xs
end

# Run the gas network for a given solver. This will run over the scenarios defined within run_timed_gas
# The values for total solver time and linear solver time were ultimately the values used within the output files
all_solver_times, solver_results = run_timed_gas(1)

