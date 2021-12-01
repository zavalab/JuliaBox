using InfiniteOpt, Gurobi, GaussianRandomFields, Interpolations, Random, Statistics, Plots, LaTeXStrings, JSON

function main()
     Random.seed!(42)

     # Prepare the random field
     num_grids = 31
     num_samples = 7
     pts = range(-1, stop = 1, length = num_grids)
     cov = CovarianceFunction(2, Matern(1/4, 1.5))
     grf = GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=625)
     samples = [max.(sample(grf) .+ 0.5, 0.0) for k in 1:num_samples]
     Ds = map(g -> LinearInterpolation((pts, pts), g), samples)

     # Plot the fields 
     hs = [heatmap(pts, pts, samples[k], legend = false, xlabel = L"$x_1$", 
               ylabel = L"$x_2$", clim = (0, 4)) for k in 1:4]
     plot(hs..., heatmap((0:0.01:1).*ones(101,1), legend=:none, xticks=:none, 
          yticks=(1:10:101, string.(0:0.4:4))), layout = @layout[grid(2,2) a{0.05w}])
     # savefig("diffusion_fields.png")

     # Solve the relaxation problem for each pareto pair
     data = Dict()
     for (i, a) in enumerate([0.0947])
          model = InfiniteModel(Gurobi.Optimizer)
          set_time_limit_sec(model, 1800.0)

          @infinite_parameter(model, t in [0, 1], num_supports = 10)
          @infinite_parameter(model, x[1:2] in [-1, 1], supports = collect(pts), independent = true)

          @variable(model, yc[1:num_samples] >= 0, Infinite(t, x))
          @variable(model, 0 <= yg <= 0.1, Infinite(t, x))
          @variable(model, 0 <= q[1:num_samples] <= 1)
          @variable(model, yc_max[1:num_samples])

          @parameter_function(model, D[k = 1:num_samples] == x -> Ds[k](x[1], x[2])) 

          @objective(model, Min, 1/num_samples * sum(q))
          @constraint(model, [k = 1:num_samples], yc_max[k] >= yc[k])
          @constraint(model, [k = 1:num_samples], yc_max[k] - 0.25 <= q[k])

          @constraint(model, pde[k = 1:num_samples], ∂(yc[k], t) == D[k] * (@∂(yc[k], x[1]^2) + @∂(yc[k], x[2]^2)) + yg)
          @constraint(model, [k = 1:num_samples], yc[k](t, [-1, x[2]]) == 0)
          @constraint(model, [k = 1:num_samples], yc[k](t, [1, x[2]]) == 0)
          @constraint(model, [k = 1:num_samples], yc[k](t, [x[1], -1]) == 0)
          @constraint(model, [k = 1:num_samples], yc[k](t, [x[1], 1]) == 0)
          @constraint(model, [k = 1:num_samples], yc[k](0, x) == 0)
          @constraint(model, cost, 1/num_samples * sum(∫(∫(∫((yc[k] - 0.2)^2, x[1]), x[2]), t) for k in 1:num_samples) <= a)

          optimize!(model)

          if has_values(model)
               sim_yc = value.(yc, ndarray = true)
               sim_yg = value(yg, ndarray = true)
               data[i] = Dict(:yc => sim_yc, :q => value.(q), :cost => value(cost), :yg => sim_yg)
          end
     end
     # open("diffusion_data.json", "w") do f
     #      JSON.print(f, data)
     # end
     return data
end

data = main()

# Helper function to plot the results as wanted
function data_plot(arr)
     p = heatmap(permutedims(arr, (2, 3, 1)), ticks = false, xlabel = L"$x_1$", 
             ylabel = L"$x_2$", zlabel = L"$t$")
     return p
end
