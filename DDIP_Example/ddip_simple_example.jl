# A Simple Example of Implementation of Dual Dynamic Integer Programming (DDIP)
# A Problem with 1 State, 2 Controls, 2 Stages.
# Written in Julia 1.0.5, JuMP 0.21
# Developed by: Ranjeet Kumar and Victor M. Zavala, UW-Madison
# Journal Article: Dual Dynamic Programming for Multi-Scale Mixed-Integer MPC, Computers and Chemical Engineering, 2020.

using JuMP; using Gurobi;

maxiter = 5; # Maximum number of DDIP iterations
tolerance = 0.01; # Percentage tolerance stopping criteria

stages = 2; # Total number of stages
blocksize = 1; # Length of each block partition
nblocks = Int(stages/blocksize); # Number of block partitions

x0 = 1; # Initial state
A = 1; B = -1; P = [0.5,0.2]; Q = [0.6,0.3]; D = 0.5; L = [0.7,0.4] # System matrices

function MPC_model(xi,L,P,Q,MILP_LP)
    N = length(L) # Horizon length of this problem
    m = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @variable(m, -1 <= u[k=1:N] <= 1)
    if (MILP_LP == "MILP")
        @variable(m, v[k=1:N], Bin)
    elseif (MILP_LP == "LP")
        @variable(m, 0 <= v[k=1:N] <= 1)
    end
    @variable(m, 0 <= x[k=1:N+1] <= 1)
    @variable(m, z[k=1:N])
    @variable(m, Cost >= 0)
    @variable(m, theta >=0)
    @constraint(m, InitialState, x[1] == xi)
    @constraint(m, Auxiliary, z[1] == xi)
    @constraint(m, state_dynamics1[k = 1], x[k+1] == A*z[k] + B*u[k])
    @constraint(m, state_dynamics2[k = 2:N], x[k+1] == A*x[k] + B*u[k])
    @constraint(m, Load[k = 1:N], u[k] + D*v[k] == L[k])
    @constraint(m, COST, Cost == sum(P[k]*u[k] + Q[k]*v[k] for k = 1:N))
    @objective(m, Min, Cost + theta)
    return m
end

##########  Defining Arrays and Dictionaries for the DDIP algorithm  ###########

LB_Vector = Vector();
UB_Vector = Vector();
best_UB_Vector = Vector();
gap_Vector = Vector();
best_gap_Vector = Vector();
mu_Vector =  Vector();
m_f = Array{Model}(undef, nblocks);
m_b = Array{Model}(undef, nblocks);
Cost_f = Dict();
x0_f = Dict();
x_f_blk = Dict();
mu = Dict();
obj = Dict();
################################################################################
#############  Dual Dynamic Integer Programming Algorithm Starts  ##############
LB = -1e10; UB = 1e10; best_UB = 1e10;
start_time = time()
# First iteration of Forward and backward sweeps - Need to initialize models in forward and backward sweeps in first iteration before the 'while' loop
j = 1; # Iteration number
x0_f[j,1] = x0;
# Forward sweep for first iteration starting
for p in 1:nblocks # Forward pass
    L_ = L[(p-1)*blocksize+1:p*blocksize];
    P_ = P[(p-1)*blocksize+1:p*blocksize];
    Q_ = Q[(p-1)*blocksize+1:p*blocksize];
    m_f[p] = MPC_model(x0_f[j,p],L_,P_,Q_,"MILP");
    optimize!(m_f[p]);
    x0_f[j,p+1] = value.(m_f[p][:x])[end];
    x_f_blk[j,p] = value.(m_f[p][:x]);
    Cost_f[j,p] = value.(m_f[p][:Cost]);
end # End forward pass
UB = sum(Cost_f[j,p] for p = 1:nblocks);
push!(UB_Vector,UB);
best_UB = minimum(UB_Vector);
push!(best_UB_Vector,best_UB);
# Forward pass for first iteration ended
# Backward pass starting for first iteration
obj[j,nblocks+1] = 0; mu[j,nblocks+1] = 0;
for p in nblocks:-1:1 # Backward pass
    L_ = L[(p-1)*blocksize+1:p*blocksize];
    P_ = P[(p-1)*blocksize+1:p*blocksize];
    Q_ = Q[(p-1)*blocksize+1:p*blocksize];
    m_b[p] = MPC_model(x0_f[j,p],L,P,Q,"LP");
    optimize!(m_b[p]);
    mu[j,p] = dual(m_b[p][:Auxiliary]);
    obj[j,p] = objective_value(m_b[p]);
end # Backward pass end
LB = value.(m_b[1][:Cost]);
push!(LB_Vector,LB)
gap = abs(UB-LB)/UB*100; # Percentage gap
push!(gap_Vector,gap);
best_gap = abs(best_UB-LB)/best_UB*100; # Percentage gap
push!(best_gap_Vector,best_gap);
println("Iteration : $(j), Current UB = $(UB), Best UB = $(best_UB), Current LB = $(LB), cur_gap = $(round(gap,digits=3))%, best_gap = $(round(best_gap,digits=3))%")
# Backward sweep for first iteration ended

#################  Iterations 2,3,... start  ###################
while best_gap >= tolerance && j < maxiter
    global j += 1; # Iteration number
    for p in 1:nblocks
        @constraint(m_f[p], m_f[p][:theta] >= obj[j-1,p+1] + mu[j-1,p+1]*(m_f[p][:x][end] - x0_f[j-1,p+1])) # Cut added to forward problem
        @constraint(m_b[p], m_b[p][:theta] >= obj[j-1,p+1] + mu[j-1,p+1]*(m_b[p][:x][end] - x0_f[j-1,p+1])) # Cut added to backward problem
    end
    x0_f[j,1] = x0;
    for p in 1:nblocks # Forward pass
        set_normalized_rhs(m_f[p][:Auxiliary], x0_f[j,p])
        optimize!(m_f[p]);
        x0_f[j,p+1] = value.(m_f[p][:x])[end];
        x_f_blk[j,p] = value.(m_f[p][:x]);
        Cost_f[j,p] = value.(m_f[p][:Cost]);
    end # End forward pass
    UB = sum(Cost_f[j,p] for p = 1:nblocks);
    push!(UB_Vector,UB);
    best_UB = minimum(UB_Vector);
    push!(best_UB_Vector,best_UB);
    # Forward pass for iteration j ended
    # Backward pass starting for iteration j
    obj[j,nblocks+1] = 0; mu[j,nblocks+1] = 0;
    for p in nblocks:-1:1 # Backward pass
        set_normalized_rhs(m_b[p][:Auxiliary], x0_f[j,p])
        optimize!(m_b[p]);
        mu[j,p] = dual(m_b[p][:Auxiliary]);
        obj[j,p] = objective_value(m_b[p]);
    end # Backward pass end
    LB = value.(m_b[1][:Cost]);
    push!(LB_Vector,LB)
    gap = abs(UB-LB)/UB*100; # Percentage gap
    push!(gap_Vector,gap);
    global best_gap = abs(best_UB-LB)/best_UB*100;
    push!(best_gap_Vector,best_gap);
    println("Iteration : $(j), Current UB = $(UB), Best UB = $(best_UB), Current LB = $(LB), cur_gap = $(round(gap,digits=3))%, best_gap = $(round(best_gap,digits=3))%")
    # Backward pass for j iteration ended
end # Iteration j ends
time_taken_dual_dynamic = start_time - time();
################################################################################
(best_ub,best_index) = findmin(UB_Vector);
# DDIP solution - State evolution
x_ddip = Vector(undef,stages+1)
x_ddip[1] = x0;
for p in 1:nblocks
    x_ddip[2+(p-1)*blocksize:1+p*blocksize] = x_f_blk[best_index,p][2:end];
end
################################################################################
#########  FOR COMPARISON WITH OPTIMAL SOLUTION WITH FULL PROBLEM  #############
m_opt = MPC_model(x0,L,P,Q,"MILP")
optimize!(m_opt)
time_taken_opt = solve_time(m_opt)
Cost_opt = value.(m_opt[:Cost])
x_opt = value.(m_opt[:x])
################################################################################
gap_opt = abs((UB_Vector[end]-Cost_opt)/Cost_opt*100); # Percentage gap
best_gap_opt = abs((best_UB_Vector[end]-Cost_opt)/Cost_opt*100); # Percentage gap
################################################################################
