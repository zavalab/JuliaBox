# This was a function to make the 13 pipe plot for both the optimization and simulation results.
# The simulation stuff is commented out.

using PyPlot

flow_profile = readdlm("./results/linkflowresults.tab")
press_profile = readdlm("./results/linkpressresults.tab")

# sim_press_profile10 = readdlm("./results/pressure_dmnetwork10.dat")
# sim_press_profile50 = readdlm("./results/pressure_dmnetwork50.dat")

function plot_13_pipe(flow_profile,Nx,Ntime;label = nothing,color = nothing,p = :line)
    #flowdata is a nlinkxntime matrix
    dims = size(flow_profile)
    n_links = dims[1]
    @assert dims[2] == Nx*Ntime
    flow_matrix = zeros(Ntime,Nx*n_links)
    j = 0
    for t = 1:Ntime
        istart = 1 + j*Nx
        iend = istart+Nx - 1
        full_link = vec(flow_profile[:,istart:iend]')
        flow_matrix[t,:] = full_link
        j += 1
    end
    total_distance = n_links*100
    n_points = size(flow_matrix)[2]
    x_plot = collect(linspace(0,total_distance,n_points))
    for t = 1:1
        if p == :line
            plot(x_plot,flow_matrix[t,:]', linewidth = 2,label = label, color =color)
        else
            scatter(x_plot,flow_matrix[t,:]', linewidth = 2,label = label, color =color)
        end
    end
    grid("on")
    xlabel("Distance [km]", fontsize = 18)
    ylabel("Pressure [bar]", fontsize = 18)
    ylim(ymin=32)
    # if sim_profile != nothing && sim_nx != nothing && sim_Ntime != nothing
    #     dims = size(sim_profile)
    #     n_links = dims[1]
    #     @assert dims[2] == sim_nx*sim_Ntime
    #     flow_matrix = zeros(sim_Ntime,sim_nx*n_links)
    #     j = 0
    #     for t = 1:Ntime
    #         istart = 1 + j*sim_nx
    #         iend = istart+sim_nx - 1
    #         full_link = vec(sim_profile[:,istart:iend]')
    #         flow_matrix[t,:] = full_link
    #         j += 1
    #     end
    #     total_distance = n_links*100
    #     n_points = size(flow_matrix)[2]
    #     x_plot = collect(linspace(0,total_distance,n_points))
    #     for t = 1:1
    #         if t == 1
    #             plot(x_plot,flow_matrix[t,:]', linewidth = 2 ,label = "Petsc Verification")
    #         else
    #             plot(x_plot,flow_matrix[t,:]', linewidth = 2 )
    #         end
    #     end
    #     legend()
    # end
    legend()
end

#plot_13_pipe(press_profile,length(x_grid),length(time_grid), label = "Optimization - 3 points", color = "blue",p = :scatter)
# plot_13_pipe(sim_press_profile10,10,48;label = "Simulation - 10 points", color = "green")
# plot_13_pipe(sim_press_profile50,50,48;label = "Simulation - 50 points", color = "red")
#plot_13_pipe(press_profile,length(x_grid),length(time_grid),sim_press_profile,20,48)
