using PyPlot

# Storing the result in a dummy variable
x = ones(length(SHOLDERS))
i = 1
nbins = 20

for sh in SHOLDERS
x[i] = getvalue(dissatisfaction[sh])
i = i+1
end

plt[:hist](x, nbins)
xlabel("Dissatisfaction")
ylabel("Frequency")
savefig("figure_histogram_alpha_" * "$(file_str)" * ".pdf")
#savefig("Histogram_alpha_" * "$(alpha)" * ".pdf")
close()
