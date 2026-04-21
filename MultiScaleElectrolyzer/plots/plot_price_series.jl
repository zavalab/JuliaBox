"""
Plot time series for spp_arr_DAM, spp_arr_DAM_pan, and distrib_10.
Uses functions defined in electricity_markets.jl.

Run from repo root:
    julia --project=. plots/plot_price_series.jl
"""

include("electricity_markets.jl")  

using Plots

fig_dir = joinpath(@__DIR__, "figures/price")

# ── 1. Full time series: DAM average ─────────────────────────────────────────
p1 = plot_spp(spp_arr_DAM; title="ERCOT DAM Average Price (10-Year)")
savefig(p1, joinpath(fig_dir, "dam_avg_full.png"))

# ── 2. Full time series: Panhandle ────────────────────────────────────────────
p2 = plot_spp(spp_arr_DAM_pan; title="ERCOT DAM Panhandle Price")
savefig(p2, joinpath(fig_dir, "dam_pan_full.png"))

# ── 3. distrib_10: annual-average-based price proxy (11 years × 8760 hrs) ────
p3 = plot_spp(distrib_10; title="Annual Average Price Distribution (distrib\\_10)",
              ylabel="Price (\$/MWh)", linewidth=1.2)
savefig(p3, joinpath(fig_dir, "distrib10_full.png"))

# ── 4. First-week overlay: DAM avg, Panhandle, distrib_10 ────────────────────
week = 24 * 7
p4 = plot_spp_compare(
    [spp_arr_DAM[1:week], spp_arr_DAM_pan[1:week], distrib_10[1:week]],
    ["DAM Average", "DAM Panhandle", "Distribution Pricing"];
    title="Electricity Prices – First Week",
)
savefig(p4, joinpath(fig_dir, "first_week_comparison.png"))

# ── 4b. First week / last week (individual plots) ────────────────────────────
n_min  = min(length(spp_arr_DAM), length(spp_arr_DAM_pan), length(distrib_10))
fw_arrs = [spp_arr_DAM[1:week], spp_arr_DAM_pan[1:week], distrib_10[1:week]]
lw_arrs = [spp_arr_DAM[n_min-week+1:n_min], spp_arr_DAM_pan[n_min-week+1:n_min], distrib_10[n_min-week+1:n_min]]
shared_ylims = (minimum(minimum.(vcat(fw_arrs, lw_arrs))), maximum(maximum.(vcat(fw_arrs, lw_arrs))))
p4b_L = plot_spp_compare(fw_arrs, ["DAM Average", "DAM Panhandle", "Distribution Pricing"];
    title="First Week", ylims=shared_ylims, legend=:topleft)
p4b_R = plot_spp_compare(lw_arrs, ["DAM Average", "DAM Panhandle", "Distribution Pricing"];
    title="Last Week", ylims=shared_ylims, legend=:topleft)
savefig(p4b_L, joinpath(fig_dir, "first_week.pdf"))
savefig(p4b_R, joinpath(fig_dir, "last_week.pdf"))

# data_dir = joinpath(@__DIR__, "..", "data")
# A_avg, zr_avg = load_result(joinpath(data_dir, "test", "basecase_avg.JSON"))

# p4c_L = plot_operation_compare

# # ── 5. Overlay: DAM avg vs Panhandle (first 8760 hrs) ────────────────────────
# n = min(length(spp_arr_DAM), length(spp_arr_DAM_pan), 8760)
# p5 = plot_spp_compare(
#     [spp_arr_DAM[1:n], spp_arr_DAM_pan[1:n]],
#     ["DAM Average", "Panhandle"];
#     title="DAM Average vs Panhandle (Year 1)"
# )
# savefig(p5, joinpath(fig_dir, "dam_avg_vs_pan_yr1.png"))

# # ── 6. Duration curves ────────────────────────────────────────────────────────
# p6a = plot_spp_duration_curve(spp_arr_DAM;     title="Duration Curve – DAM Average")
# p6b = plot_spp_duration_curve(spp_arr_DAM_pan; title="Duration Curve – Panhandle")
# p6c = plot_spp_duration_curve(distrib_10;      title="Duration Curve – distrib\\_10")
# p6  = plot(p6a, p6b, p6c; layout=(1,3), size=(1400, 400))
# savefig(p6, joinpath(fig_dir, "duration_curves.png"))

# # ── 7. Histograms ─────────────────────────────────────────────────────────────
# p7a = plot_spp_histogram(spp_arr_DAM;     title="Histogram – DAM Average")
# p7b = plot_spp_histogram(spp_arr_DAM_pan; title="Histogram – Panhandle")
# p7c = plot_spp_histogram(distrib_10;      title="Histogram – distrib\\_10")
# p7  = plot(p7a, p7b, p7c; layout=(1,3), size=(1400, 400))
# savefig(p7, joinpath(fig_dir, "histograms.png"))

# # ── 8. Annual average bar charts ──────────────────────────────────────────────
# p8a = plot_spp_yearly_avg(spp_arr_DAM;  title="Annual Avg – DAM Average")
# p8b = plot_spp_yearly_avg(distrib_10;   title="Annual Avg – distrib\\_10")
# p8  = plot(p8a, p8b; layout=(1,2), size=(1000, 400))
# savefig(p8, joinpath(fig_dir, "annual_averages.png"))

# println("All plots saved to plots/figures/")
