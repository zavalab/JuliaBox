# Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (zavalatejeda@wisc.edu)
# Functions for the parameter estimation of Lotka-Volterra model

# Makes JuMP/Plasmo model for the parameter estimation
function param(data::Dict,args::Dict)

    # Count total number of data points
    n_datapoint = sum(data[exp][i][ii][:n_time]*data[exp][i][ii][:n_species] for exp in keys(data) for i in 1:length(data[exp]) for ii in 1:length(data[exp][i]))
    
    # Create a graph-based parameter estimation model
    graph=Plasmo.PlasmoGraph()

    # Create a parent model
    m_parent=Model()

    # Define paramets in parent node
    @variable(m_parent, -args[:maxr]/args[:r_scale]<=r[i=1:args[:n_total_species],j=0:args[:n_total_species]]<=args[:maxr]/args[:r_scale], start=args[:r_start][i,j]/args[:r_scale])

    # Additional parameter if the model is Saturable form
    if args[:model] == :Sat
        @variable(m_parent, 0<=K[i=1:args[:n_total_species],j=0:args[:n_total_species]]<=args[:maxr], start=args[:r_start][i,j])
        @expression(m_parent, reg_K, args[:lambda]*sum((K[i,j])^2 for i=1:args[:n_total_species], j=1:args[:n_total_species]))
    else
        reg_K = .0
    end

    # Additional slack variable if CVaR formulation
    if args[:min_norm]==:CVaR
        @variable(m_parent,t,start=0.0)
    else
        t = 0.0
    end

    # Priors
    if args[:reg_norm]==:L1
        @variable(m_parent, r_minnorm[i=1:args[:n_total_species],j=0:args[:n_total_species]], start=abs(args[:r_start][i,j])/args[:r_scale])
        @constraint(m_parent, [i=1:args[:n_total_species],j=0:args[:n_total_species]], r_minnorm[i,j]*args[:r_scale]>=r[i,j]*args[:r_scale]-args[:r_prior][i,j])
        @constraint(m_parent, [i=1:args[:n_total_species],j=0:args[:n_total_species]], r_minnorm[i,j]*args[:r_scale]>=-r[i,j]*args[:r_scale]+args[:r_prior][i,j])
        @expression(m_parent, reg_r, args[:lambda]*sum(r_minnorm[i,j] for i=1:args[:n_total_species], j=0:args[:n_total_species]))
    elseif args[:reg_norm]==:L2
        @expression(m_parent, reg_r, args[:lambda]*sum(((r[i,j]*args[:r_scale]-args[:r_prior][i,j]))^2 for i=1:args[:n_total_species], j=0:args[:n_total_species]))
    else
        reg_r = 0
    end

    # Parent model objective function (include priors and t in CVaR formulation)
    @objective(m_parent,Min, reg_r + reg_K + t * n_datapoint)

    # Connect the parent model to the graph model
    node_parent=add_node!(graph,m_parent)

    # Define models and nodes for children models and nodes
    m_children=Dict(exp=>[] for exp in keys(data))
    node_children=Array{NodeOrEdge,1}()

    # For each experiment,
    for exp in keys(data)
        for i in 1:length(data[exp])
            # Create a children model
            m = Model()
            d = data[exp][i]

            # Define parameteargs[:r_scale] in children node
            @variable(m, -args[:maxr]/args[:r_scale]<=r[j in d[1][:species] ,k in [0;d[1][:species]]]<=args[:maxr]/args[:r_scale],
                      start=args[:r_start][j,k]/args[:r_scale])
            if args[:model] == :Sat
                @variable(m, args[:eps]<=K[j in d[1][:species] ,k in [0;d[1][:species]]]<=args[:maxr], start=1)
            end
            # Define state variables in children node
            @variable(m, y[ii in 1:length(data[exp][i]), j in 1:d[ii][:n_species], k in 1:d[ii][:n_disc_time]],
                      start=args[:y_start][exp,i,ii,j,k]/args[:y_scale][exp][i][ii,j,k])

            # Discretized nonlinear ODE
            if args[:model]==:Sat
                @NLconstraint(m, [ii in 1:length(data[exp][i]), j in 1:d[ii][:n_species], k in 2:d[ii][:n_disc_time]],
                              y[ii,j,k]*args[:y_scale][exp][i][ii,j,k]==y[ii,j,k-1]*args[:y_scale][exp][i][ii,j,k-1]
                              +(r[d[ii][:species][j],0]*args[:r_scale] + sum(r[d[ii][:species][j],d[ii][:species][p]]*args[:r_scale]*y[ii,p,k]*args[:y_scale][exp][i][ii,p,k]
                                                                            /(K[d[ii][:species][j],d[ii][:species][p]]+y[ii,p,k]*args[:y_scale][exp][i][ii,p,k])
                                                                            for p in 1:d[ii][:n_species]))
                              *y[ii,j,k]*args[:y_scale][exp][i][ii,j,k]*(d[ii][:disc_time][k]-d[ii][:disc_time][k-1]))
            else
                @NLconstraint(m,
                              [ii in 1:length(data[exp][i]), j in 1:d[ii][:n_species], k in 2:d[ii][:n_disc_time]],
                              y[ii,j,k]*args[:y_scale][exp][i][ii,j,k]==y[ii,j,k-1]*args[:y_scale][exp][i][ii,j,k-1]
                              +(r[d[ii][:species][j],0]*args[:r_scale]+sum(r[d[ii][:species][j],d[ii][:species][p]]*y[ii,p,k]*args[:r_scale]*
                                                                           args[:y_scale][exp][i][ii,p,k] for p in 1:d[ii][:n_species]))
                              *y[ii,j,k]*args[:y_scale][exp][i][ii,j,k]*(d[ii][:disc_time][k]-d[ii][:disc_time][k-1]))
            end
            # Initial conditions
            @constraint(m, [ii=1, j in 1:d[ii][:n_species]], y[1,j,1]*args[:y_scale][exp][i][1,j,1]==d[ii][:y0][j])

            # Initial conditions for special cases (diluting the solution)
            @constraint(m, [ii in 2:length(data[exp][i]),
                            j in 1:d[ii][:n_species]], y[ii,j,1]*args[:y_scale][exp][i][ii,j,1]==
                        y[ii-1,j,d[ii][:n_disc_time]]*args[:y_scale][exp][i][ii-1,j,d[ii][:n_disc_time]]*args[:dilution_rate])
            
            # Objective function for children nodes (mean square error)
            if args[:min_norm]==:L1
                @objective(m, Min,
                           sum(((y[ii,j,(kk-1)*args[:n_disc][exp]+1]*args[:y_scale][exp][i][ii,j,(kk-1)*args[:n_disc][exp]+1]-d[ii][:y][j][kk])
                                /args[:y_sig_abs][exp,i,ii,j,kk]/args[:y_sig_rel])^2
                               for ii in 1:length(data[exp][i]) for j in 1:d[1][:n_species] for kk in 1:d[ii][:n_time]))
            elseif args[:min_norm]==:CVaR
                @variable(m,t,start=0)
                @variable(m,ts[ii in 1:length(data[exp][i]),j in 1:d[1][:n_species],kk in 1:d[ii][:n_time]]>=0,start=0)
                @constraint(m,[ii in 1:length(data[exp][i]), j in 1:d[1][:n_species],kk in 1:d[ii][:n_time]],
                            ts[ii,j,kk]>=
                            ((y[ii,j,(kk-1)*args[:n_disc][exp]+1]*args[:y_scale][exp][i][ii,j,(kk-1)*args[:n_disc][exp]+1]-d[ii][:y][j][kk])
                             /args[:y_sig_abs][exp,i,ii,j,kk]/args[:y_sig_rel])^2-t)
                @objective(m, Min, sum(1/(1-args[:beta])*ts[ii,j,kk] for ii in 1:length(data[exp][i]) for j in 1:d[1][:n_species] for kk in 1:d[ii][:n_time]))
            end
            # Save the children model
            push!(m_children[exp],m)

            # Connect the children model to the graph model
            push!(node_children,add_node!(graph,m))
            
            # Linking constraints (keep parameter values the same)
            @linkconstraint(graph,[j in d[1][:species], k in [0;d[1][:species]]],
                            m_parent[:r][j,k]==m[:r][j,k])

            if args[:model] == :Sat
                @linkconstraint(graph,[j in d[1][:species], k in d[1][:species]],
                                m_parent[:K][j,k]==m[:K][j,k])
            end
            if args[:min_norm] == :CVaR
                @linkconstraint(graph,  m_parent[:t]==m[:t])
            end

        end
    end

    return graph,m_parent,m_children,node_parent,node_children
end

# Performs inference analysis and saves the samples from posterior
function inference_analysis(data::Dict,args::Dict,opts::Dict,max_sample::Int64)
    
    # Solve the original problem (to provide warm start)
    mkpath(opts[:outputpath])
    graph,m_parent,m_children = param(data,args)
    graph.solver=opts[:solver]
    Plasmo.solve(graph)

    # Warm start using the solution of original problem
    args[:r_start]=Dict((i,j)=>getvalue(m_parent[:r][i,j]) for i=1:args[:n_total_species] for j=0:args[:n_total_species])
    args[:y_start]=Dict((exp,i,ii,j,k)=>getvalue(m_children[exp][i][:y][ii,j,k])
                 for exp in keys(data) for i in 1:length(data[exp]) for ii in 1:length(data[exp][i])
                 for j in 1:data[exp][i][ii][:n_species] for k in 1:data[exp][i][ii][:n_disc_time])

    # Repeat sampling
    for sampleindex in 1:max_sample
        
        # Seeding random variable
        srand()

        # Duplicate data dictionary
        data_randomized = deepcopy(data)
        
        # Randomize the data
        for exp in keys(data)
            for i in 1:length(data[exp])
                for ii in 1:length(data[exp][i])
                    for j=1:data[exp][i][ii][:n_species]
                        for k=1:data[exp][i][ii][:n_time]
                            data_randomized[exp][i][ii][:y][j][k] +=
                                max(data_randomized[exp][i][ii][:y][j][k],args[:y_sig_min]) * args[:y_sig_rel] * randn()
                        end
                    end
                end
            end
        end

        # Randomize the prior
        for i=1:args[:n_total_species]
            for j=0:args[:n_total_species]
                args[:r_prior][i,j] = randn()/sqrt(args[:lambda])
            end
        end

        # Solve the problem
        graph,m_parent,m_children = param(data,args)
        graph.solver=IpoptSolver(linear_solver="ma57",max_iter=30)
        Plasmo.solve(graph)
        
        # Save output
        rValue=getvalue(m_parent[:r])[:,:]'[:]*args[:r_scale]        
        mkpath(opts[:outputpath]*"posterior")
        ind = length(readdir(opts[:outputpath]*"posterior"))+1 # This allows running multiple scripts at the same time
        writecsv(opts[:outputpath]*"posterior/param$ind.csv",rValue[:,:])
    end
end


# Saves output as *.csv files
function output(data::Dict,args::Dict,opts::Dict,outs::Dict)
    # Get the solution
    rValue=getvalue(outs[:m_parent][:r])[:,:]'[:]*args[:r_scale]
    yValue=Dict(exp=>[[[[getvalue(outs[:m_children][exp][i][:y][ii,j,k] * args[:y_scale][exp][i][ii,j,k])
                         for k in 1:data[exp][i][ii][:n_disc_time]]
                        for j in 1:data[exp][i][ii][:n_species]]
                       for ii=1:length(data[exp][i])] for i=1:length(data[exp])]
                for exp in keys(data))

    # Save the parameters
    writecsv(opts[:outputpath]*"param.csv",rValue[:,:])
    if args[:model] == :Sat
        kValue=getvalue(outs[:m_parent][:K])[:,:]'[:]
        writecsv(opts[:outputpath]*"param_K.csv",kValue[:,:])
    end

    # Save the states
    if opts[:save_y]
        yout=Dict(exp=>Array{Any}(length(data[exp])) for exp in keys(data))
        yexp=Dict(exp=>Array{Any}(length(data[exp])) for exp in keys(data))
        tout=Dict(exp=>Array{Any}(length(data[exp])) for exp in keys(data))
        texp=Dict(exp=>Array{Any}(length(data[exp])) for exp in keys(data))

        for exp in keys(data)
            for i in 1:length(data[exp])
                yexp[exp][i]=[data[exp][i][1][:y][j] for j in 1:data[exp][i][1][:n_species]]
                yout[exp][i]=[yValue[exp][i][1][j] for j in 1:data[exp][i][1][:n_species]]
                texp[exp][i]=data[exp][i][1][:time]
                tout[exp][i]=data[exp][i][1][:disc_time]
                for ii in 2:length(data[exp][i])
                    append!(texp[exp][i], data[exp][i][ii][:time]+texp[exp][i][end]+1e-10)
                    append!(tout[exp][i], data[exp][i][ii][:disc_time]+tout[exp][i][end]+1e-10)
                    for j=1:length(data[exp][i][ii][:species])
                        append!(yexp[exp][i][j],data[exp][i][ii][:y][j])
                        append!(yout[exp][i][j],yValue[exp][i][ii][j])
                    end
                end

                youtPath=opts[:outputpath]*"yout/"*string(exp)*"/y$i/"
                mkpath(youtPath)
                writecsv(youtPath*"species.csv",data[exp][i][1][:species][:,:])
                writecsv(youtPath*"timePoints.csv",texp[exp][i][:,:])
                writecsv(youtPath*"timePointsExtended.csv",tout[exp][i][:,:])
                for j=1:length(data[exp][i][1][:species])
                    writecsv(youtPath*"yexp$j.csv",yexp[exp][i][j][:,:])
                    writecsv(youtPath*"yout$j.csv",yout[exp][i][j][:,:])
                end
            end
        end
    end

    # Save the prediction errors
    if opts[:save_e]
        mses=[]
        mses_exp=[]
        mse = Dict()
        mse_exp=Dict()
        for exp in keys(data)
            for i in 1:length(data[exp])
                # mse=0
                for j in 1:data[exp][i][1][:n_species]
                    mse_exp_1 = 0
                    for ii in 1:length(data[exp][i])
                        for kk in 1:data[exp][i][ii][:n_time]
                            y_val=yValue[exp][i][ii][j][(kk-1)*args[:n_disc][exp]+1]
                            mse[exp,i,ii,j,kk]=((y_val-data[exp][i][ii][:y][j][kk])/args[:y_sig_abs][exp,i,ii,j,kk]/args[:y_sig_rel])^2
                            mse_exp_1+=mse[exp,i,ii,j,kk]
                            push!(mses,mse[exp,i,ii,j,kk])
                        end
                    end
                    mse_exp[exp,i,j] = mse_exp_1
                    push!(mses_exp,mse_exp_1)
                end
            end
        end
        mse_cr = sort(mses)[end-10]
        mse_exp_cr = sort(mses_exp)[end-10]
        mse_ind = []
        mse_exp_ind=[]

        for exp in keys(data)
            for i in 1:length(data[exp])
                for j in 1:data[exp][i][1][:n_species]
                    for ii in 1:length(data[exp][i])
                        for kk in 1:data[exp][i][ii][:n_time]
                            if mse[exp,i,ii,j,kk] > mse_cr
                                push!(mse_ind,
                                      "exp:"*string(exp)*"i:"*string(i)*",species:"*string(data[exp][i][ii][:species])*
                                      ",phase:"*string(ii)*",spind:"*string(j)*",timeind:"*string(kk))
                            end
                        end
                    end
                    if mse_exp[exp,i,j] > mse_exp_cr
                        push!(mse_exp_ind,
                              "exp:"*string(exp)*"i:"*string(i)*",species:"*string(data[exp][i][1][:species])*
                              ",spind:"*string(j))
                    end
                end
            end
        end
        writecsv(opts[:outputpath]*"mses.csv",mses[:,:])
    end

    # Save additional infos
    touch(opts[:outputpath]*"out.out")
    out=open(opts[:outputpath]*"out.out","w")
    println(out, "* Elapsed Time:")
    println(out, outs[:time_elapsed]/1.0e9) # Computation time
    println(out, "* Large error indices:") # Large prediction error indices
    for i=1:length(mse_ind)
        println(out,mse_ind[i])
    end
    println(out, "* Large total error indices:") # Large totla prediction error indices
    for i=1:length(mse_ind)
        println(out,mse_exp_ind[i])
    end
    close(out)
end

# Adds information about discretization mesh on data dictionary.
function add_discretize_info!(data::Dict,args::Dict)
    for exp in keys(data)
        for i=1:length(data[exp])
            for ii=1:length(data[exp][i])
                disc_time=zeros((data[exp][i][ii][:n_time]-1)*args[:n_disc][exp]+1)
                for t in 1:data[exp][i][ii][:n_time]-1
                    disc_time[(t-1)*args[:n_disc][exp]+1:t*args[:n_disc][exp]+1]=linspace(data[exp][i][ii][:time][t],data[exp][i][ii][:time][t+1],args[:n_disc][exp]+1)
                end
                data[exp][i][ii][:disc_time] = disc_time
                data[exp][i][ii][:n_disc_time] = (data[exp][i][ii][:n_time]-1)*args[:n_disc][exp]+1
            end
        end
    end
end
