using Plots, DelimitedFiles
pgfplotsx()

record = reshape(readdlm("output/adchem_record.csv",',',Float64),3,3,4,7)
numvar = reshape(readdlm("output/adchem_numvar.csv",',',Int),4,7)


dataset = [
    [
        numvar[1,:],
        [record[1,:,1,:]' record[1,:,2,:]'],
        "Solution Wall Time (sec)",
        "gas-total"
    ],
    [
        numvar[1,:],
        [record[2,:,1,:]' record[2,:,2,:]'],
        "Linear Solver Wall Time (sec)",
        "gas-linear"
    ],
    [
        numvar[1,:],
        [record[3,:,1,:]' record[3,:,2,:]'],
        "Function Eval. Wall Time (sec)",
        "gas-function"
    ],
    [
        numvar[3,:],
        [record[1,:,3,:]' record[1,:,4,:]'],
        "Solution Wall Time (sec)",
        "power-total"
    ],
    [
        numvar[3,:],
        [record[2,:,3,:]' record[2,:,4,:]'],
        "Linear Solver Wall Time (sec)",
        "power-linear"
    ],
    [
        numvar[3,:],
        [record[3,:,3,:]' record[3,:,4,:]'],
        "Function Eval. Wall Time (sec)",
        "power-function"
    ],
]

label = ["JuMP-Ma57" "JuMP-PardisoMKL" "JuMP-Schwarz/Ma57" "Plasmo-Ma57" "Plasmo-PardisoMKL" "Plasmo-Schwarz/Ma57"]
marker= [:circle :square :diamond :circle :square :diamond]
linestyles=[:solid :solid :solid :dash :dash :dash]
markersize= [5 5 5 3 3 3]

for (xdata,ydata,ylabel,fname) in dataset
    plt = plot(
        xdata,ydata,label=label,xlabel="Number of Variables", size = (500,250), ylabel=ylabel, 
        legend=:topleft, framestyle=:box, xscale=:log10, linestyles=linestyles,
        marker=marker, markersize=markersize,
        xtick=[1e4,1e5,1e6],
        ylims = (0,maximum(ydata)), xlims = (minimum(xdata),maximum(xdata)));
    savefig(plt,"output/adchem-$(fname).pdf")
end

for (xdata,ydata,ylabel,fname) in dataset
    plt = plot(
        ydata,xdata,label=label,ylabel="Number of Variables", size = (500,250), xlabel=ylabel, 
        legend=:bottomright, framestyle=:box, yscale=:log10, linestyles=linestyles,
        marker=marker, markersize=markersize,
        ytick=[1e4,1e5,1e6],
        xlims = (0,maximum(ydata)), ylims = (minimum(xdata),maximum(xdata)));
    savefig(plt,"output/adchem-reverse-$(fname).pdf")
end

# To generate tables
function print_latex_tabular(tbl)
    for i=1:size(tbl,1)
        for j=1:size(tbl,2)
            print(tbl[i,j],j!=size(tbl,2) ? "& " : "\\\\\n\\hline\n")
        end
    end
end

tbl1 = Any[["1 days";"3 days";"7 days";"14 days";"30 days";"60 days";"180 days"] numvar[1,:] round.(reshape(record[:,:,1,:][:],9,:)',digits=2) round.(reshape(record[:,:,2,:][:],9,:)',digits=2)]
tbl2 = Any[["1 days";"3 days";"7 days";"14 days";"30 days";"60 days";"180 days"] numvar[3,:] round.(reshape(record[:,:,3,:][:],9,:)',digits=2) round.(reshape(record[:,:,4,:][:],9,:)',digits=2)]
println("Table 1: ")
print_latex_tabular(tbl1)
println("Table 2: ")
print_latex_tabular(tbl2)

