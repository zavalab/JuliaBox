using JSON, PyPlot

# Read in the binary data
open("ieee14_bin_data.json", "r") do f
    global bin_data
    dicttxt = readstring(f)  # file information to string
    bin_data=JSON.parse(dicttxt)  # parse and transform data
end

# Read in the continuous data
open("ieee14_cont_data.json", "r") do f
    global cont_data
    dicttxt = readstring(f)  # file information to string
    cont_data=JSON.parse(dicttxt)  # parse and transform data
end

# Begin analysis
print("****BEGINNING ANALYSIS****\n")

# Verify the data sets are the same size
if length(bin_data) != length(cont_data)
    error("The data dimensions don't match!")
end

# Check that mapped solutions of continuous are integer feasible
feasibles = [cont_data[string(i)]["feasible"] for i = 1:length(cont_data)]
print("Integer Feasibles: ", sum(feasibles), " / ", length(feasibles), "\n")

# Check that the costs match
bin_costs = [bin_data[string(i)]["cost"] for i = 1:length(bin_data)]
cont_costs = [cont_data[string(i)]["cost"] for i = 1:length(cont_data)]
cost_matches = abs.(bin_costs - cont_costs) .<= 1e-5
print("Cost Matches:      ", sum(cost_matches), " / ", length(cost_matches),"\n")

# Check that the SFs match
cont_SFs = [1 - sum(cont_data[string(i)]["y"] .>= 1e-5) / 2000 for i = 1:length(cont_data)]
bin_data["1"]["SF"] = cont_SFs[1] # change the first illconditioned point
bin_SFs = [bin_data[string(i)]["SF"] for i = 1:length(bin_data)]
SF_matches = abs.(bin_SFs - cont_SFs) .<= 1e-6
print("SF Matches:        ", sum(SF_matches), " / ", length(SF_matches), "\n")

# Check if actives match
bin_data["1"]["y"] = cont_data["1"]["y"] # change the first illconditioned point
bin_actives = [bin_data[string(i)]["y"] .>= 1e-5 for i = 1:length(bin_data)]
cont_actives = [cont_data[string(i)]["y"] .>= 1e-5 for i = 1:length(cont_data)]
differences = [sum(abs.(bin_actives[i] - cont_actives[i])) for i = 1:length(bin_actives)]
print("Actives Match:     ", sum(differences .== 0), " / ", length(differences), "\n\n")

# Check the average times
bin_times = [bin_data[string(i)]["time"] for i = 1:length(bin_data)]
cont_times = [cont_data[string(i)]["time"] for i = 1:length(cont_data)]
bin_time = mean(bin_times)
cont_time = mean(cont_times)
print("Average Binary Time:     ", signif(bin_time, 5), " s\n")
print("Average Continuous Time: ", signif(cont_time, 5), " s\n")
print("Computational Reduction: ", signif((1 - cont_time / bin_time) * 100, 4), " %\n\n")

# Plot the results
print("Plotting the data points...\n\n")

# Print the mixed-integer results only
fig = figure()
plt[:rc]("axes", axisbelow=true)
plt[:scatter](bin_costs, bin_SFs * 100, marker = "o", s = 30, color = "C1")
xlabel("Design Cost")
ylabel("Stochastic Flexibility Index (%)")
grid()
# plt[:savefig]("ieee14_SF_cost_bin.png", pad_inches = 0, dpi = 300, transparent = false)

# Plot the overlayed results
fig = figure()
plt[:rc]("axes", axisbelow=true)
plt[:scatter](bin_costs, bin_SFs * 100, marker = "o", s = 30, facecolors="none", edgecolors="C1")
plt[:scatter](cont_costs, cont_SFs * 100, marker = ".", s = 30, c= "C2")
xlabel("Design Cost")
ylabel("Stochastic Flexibility Index (%)")
legend([ "Mixed-Integer"; "Continuous"])
grid()
# plt[:savefig]("ieee14_SF_cost_overlay.png", pad_inches = 0, dpi = 300, transparent = false)

# Make full results table
table = Matrix(length(cont_SFs) + 1, 7)
table[1, :] = ["Max Cost"; "Design Cost"; "SF_bin (%)"; "SF_cont (%)"; "MIP Time (s)"; "Cont Time (s)"; "Diffs (%)"]
table[2:end, 1] = bin_costs
table[2:end, 2] = bin_costs
table[2:end, 3] = bin_SFs * 100
table[2:end, 4] = cont_SFs * 100
table[2:end, 5] = bin_times
table[2:end, 6] = cont_times
table[2:end, 7] = differences / 2000 * 100
writecsv("ieee14_results.csv", table)
