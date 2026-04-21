"""
Plot ERCOT DAM electricity price data for a specific date range (2022–2024).
Reuses plot functions from electricity_markets.jl.

Run from repo root:
    julia --project=. plots/plot_year_slice.jl
"""

include("electricity_markets.jl")   # loads spp_arr_DAM, spp_arr_DAM_pan

using Plots
using XLSX, DataFrames, Dates

fig_dir = joinpath(@__DIR__, "figures/price")
mkpath(fig_dir)

# ── Load date column (same row order as spp_arr_DAM / spp_arr_DAM_pan) ────────
_df_dates = DataFrame(XLSX.readtable(
    joinpath(@__DIR__, "../data/ERCOT_DAM_2014-2024.xlsx"), "Sheet1"))
dates_all = Date.(_df_dates[:, "DeliveryDate"])

# ── Helper: first and last row index for a set of calendar years ───────────────
function year_range_indices(dates::Vector{Date}, years)
    mask    = [Dates.year(d) in years for d in dates]
    idxs    = findall(mask)
    return first(idxs), last(idxs)
end

# ── Helper: tick positions & labels at Jan 1 of each year in the slice ─────────
function year_ticks(dates_slice::Vector{Date})
    seen  = Set{Int}()
    pos   = Int[]
    labs  = String[]
    for (i, d) in enumerate(dates_slice)
        y = Dates.year(d)
        if y ∉ seen
            push!(seen, y); push!(pos, i); push!(labs, string(y))
        end
    end
    return pos, labs
end

# ── Target date range ──────────────────────────────────────────────────────────
target_years = 2022:2024
i_start, i_end = year_range_indices(dates_all, target_years)

dam_slice   = spp_arr_DAM[i_start:i_end]
pan_slice   = spp_arr_DAM_pan[i_start:i_end]
dates_slice = dates_all[i_start:i_end]
n_slice     = length(dam_slice)
week        = 24 * 7

shared_ylims_full = (minimum(min.(dam_slice, pan_slice)),
                     maximum(max.(dam_slice, pan_slice)))

# ── Plot 1: Full 2022–2024 overlay with year-boundary ticks ───────────────────
tick_pos, tick_labs = year_ticks(dates_slice)

p1 = plot_spp_compare([dam_slice, pan_slice], ["DAM Average", "DAM Panhandle"];
    title   = "ERCOT Prices 2022–2024",
    ylims   = shared_ylims_full,
    legend  = :topleft)
plot!(p1; xticks = (tick_pos, tick_labs), xlabel = "Year")
savefig(p1, joinpath(fig_dir, "2022_2024_full.png"))

# ── Plot 2: First week of 2022 ────────────────────────────────────────────────
fw_dam = dam_slice[1:week]
fw_pan = pan_slice[1:week]
shared_ylims_fw = (minimum(min.(fw_dam, fw_pan)), maximum(max.(fw_dam, fw_pan)))

p2 = plot_spp_compare([fw_dam, fw_pan], ["DAM Average", "DAM Panhandle"];
    title  = "First Week of 2022",
    ylims  = shared_ylims_fw,
    legend = :topleft)
savefig(p2, joinpath(fig_dir, "2022_first_week.png"))

# ── Plot 3: Last week of 2024 ─────────────────────────────────────────────────
lw_dam = dam_slice[end-week+1:end]
lw_pan = pan_slice[end-week+1:end]
shared_ylims_lw = (minimum(min.(lw_dam, lw_pan)), maximum(max.(lw_dam, lw_pan)))

p3 = plot_spp_compare([lw_dam, lw_pan], ["DAM Average", "DAM Panhandle"];
    title  = "Last Week of 2024",
    ylims  = shared_ylims_lw,
    legend = false)
savefig(p3, joinpath(fig_dir, "2024_last_week.png"))

# ── Plot 4: First & last week shared y-limits ─────────────────────────────────
shared_ylims_weeks = (minimum(min.(fw_dam, fw_pan, lw_dam, lw_pan)),
                      maximum(max.(fw_dam, fw_pan, lw_dam, lw_pan)))

p4_L = plot_spp_compare([fw_dam, fw_pan], ["DAM Average", "DAM Panhandle"];
    title  = "First Week of 2022",
    ylims  = shared_ylims_weeks,
    legend = :topleft)
p4_R = plot_spp_compare([lw_dam, lw_pan], ["DAM Average", "DAM Panhandle"];
    title  = "Last Week of 2024",
    ylims  = shared_ylims_weeks,
    legend = false)
savefig(p4_L, joinpath(fig_dir, "2022_first_week_shared.png"))
savefig(p4_R, joinpath(fig_dir, "2024_last_week_shared.png"))

println("Year-slice plots saved to plots/figures/price/")
