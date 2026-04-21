include(joinpath(@__DIR__, "../src/formulation.jl"))

DAM_10yr = DataFrame(XLSX.readtable(joinpath(@__DIR__, "..", "data/ERCOT_DAM_AVG_2014-2024.xlsx"), "Sheet1"))

spp_arr_DAM = Float64.(DAM_10yr[:, "SettlementPointPrice"])
# spp_arr_DAM = vcat(spp_arr_DAM, spp_arr_DAM)

DAM_10yr_pan = DataFrame(XLSX.readtable(joinpath(@__DIR__, "..", "data/ERCOT_DAM_PAN_2014-2024.xlsx"), "Sheet1"))
spp_arr_DAM_pan = Float64.(DAM_10yr_pan[:, "SettlementPointPrice"])
# spp_arr_DAM_pan = vcat(spp_arr_DAM_pan, spp_arr_DAM_pan)

yearly_avg = [6.91, 6.91, 6.76, 6.88, 6.92, 6.81, 6.67, 7.18, 8.32, 8.04, 8.13].*10 # $/MWh
distrib_10 = repeat(yearly_avg, inner=8760)
distrib_10 = distrib_10[1:length(spp_arr_DAM)]
"""
Here we are just going to be plotting the electricity price time series
"""

using Plots
using Statistics

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false,
        titlefontsize = 20, guidefontsize = 18, tickfontsize = 12, legendfontsize=14)
"""
    plot_spp(spp_arr; kwargs...)

Plot a raw SPP time series. Keyword args are forwarded to Plots.plot.

# Examples
    plot_spp(spp_arr_DAM)
    plot_spp(extendo; title="20-Year DAM Prices", ylims=(-50, 300))
"""
function plot_spp(spp_arr; title="Electricity Price (SPP)", xlabel="Hour", ylabel="Price (\$/MWh)", kwargs...)
    T = 1:length(spp_arr)
    plot(T, spp_arr;
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        legend=false,
        linewidth=0.5,
        kwargs...
    )
end

"""
    plot_spp_slice(spp_arr, t_start, t_end; kwargs...)

Plot a time slice of an SPP array between hours `t_start` and `t_end`.

# Examples
    plot_spp_slice(spp_arr_DAM, 1, 24*7)          # first week
    plot_spp_slice(extendo, 8761, 8761+24*30)      # first month of year 2
"""
function plot_spp_slice(spp_arr, t_start::Int, t_end::Int; title="SPP Slice [$(t_start):$(t_end)]", ylabel="Price (\$/MWh)", kwargs...)
    t_range = t_start:t_end
    plot(collect(t_range), spp_arr[t_range];
        title=title,
        xlabel="Hour",
        ylabel=ylabel,
        legend=false,
        linewidth=0.8,
        kwargs...
    )
end

"""
    plot_spp_yearly_avg(spp_arr; kwargs...)

Bar chart of mean hourly SPP per calendar year. Assumes hourly data starting
from the first hour of year 1.

# Examples
    plot_spp_yearly_avg(extendo)
"""
function plot_spp_yearly_avg(spp_arr; title="Annual Average SPP", ylabel="Mean Price (\$/MWh)", kwargs...)
    n_years = floor(Int, length(spp_arr) / 8760)
    means = [mean(spp_arr[(y-1)*8760+1 : y*8760]) for y in 1:n_years]
    bar(1:n_years, means;
        title=title,
        xlabel="Year",
        ylabel=ylabel,
        legend=false,
        kwargs...
    )
end

"""
    plot_spp_duration_curve(spp_arr; kwargs...)

Price duration curve (sorted descending). Useful for visualizing price
distribution and identifying high/low price hours.

# Examples
    plot_spp_duration_curve(spp_arr_DAM)
    plot_spp_duration_curve(spp_arr_DAM; ylims=(-50, 200))
"""
function plot_spp_duration_curve(spp_arr; title="Price Duration Curve", ylabel="Price (\$/MWh)", kwargs...)
    sorted = sort(spp_arr, rev=true)
    plot(sorted;
        title=title,
        xlabel="Hours (sorted)",
        ylabel=ylabel,
        legend=false,
        linewidth=1.0,
        kwargs...
    )
end

"""
    plot_spp_histogram(spp_arr; kwargs...)

Histogram of SPP values. Clips extreme outliers beyond the 99.5th percentile
by default so the distribution is visible.

# Examples
    plot_spp_histogram(spp_arr_DAM)
    plot_spp_histogram(extendo; bins=100)
"""
function plot_spp_histogram(spp_arr; title="SPP Distribution", clip_quantile=0.995, bins=80, kwargs...)
    upper = quantile(spp_arr, clip_quantile)
    clipped = filter(x -> x <= upper, spp_arr)
    histogram(clipped;
        title=title,
        xlabel="Price (\$/MWh)",
        ylabel="Count",
        bins=bins,
        legend=false,
        kwargs...
    )
end

"""
    plot_spp_compare(arrays, labels; kwargs...)

Overlay multiple SPP arrays on one plot. Arrays need not be the same length;
each is plotted against its own 1:length index.

# Examples
    plot_spp_compare(
        [spp_arr_DAM, spp_arr_DAM_pan],
        ["DAM Average", "Panhandle"]
    )
"""
function plot_spp_compare(arrays::Vector, labels::Vector{String}; title="SPP Comparison", ylabel="Price (\$/MWh)", kwargs...)
    p = plot(; xlabel="Time (Hour)", ylabel=ylabel, kwargs...)
    for (arr, lbl) in zip(arrays, labels)
        plot!(p, 1:length(arr), arr; label=lbl, linewidth=2)
    end
    p
end


function plot_operation_compare(arrays::Vector, labels::Vector{String}; ylabel="Price (\$/MWh)", kwargs...)
    p = plot(; xlabel="Time (Hour)", ylabel=ylabel, kwargs...)
    for (arr, lbl) in zip(arrays, labels)
        plot!(p, 1:length(arr), arr; label=lbl, linewidth=2, 
        fillrange    = 0,
        fillalpha    = 0.15,
        fillcolor    = :red)
    end
    p
end

function load_result(path::String)
    d = JSON.parsefile(path)
    fvr = d["final_value_results"]
    A         = Float64.(fvr["A"])
    z_replace = Float64.(fvr["z_replace"])
    return A, z_replace
end