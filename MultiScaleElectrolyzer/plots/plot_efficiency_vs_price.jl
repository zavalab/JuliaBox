"""
Plot cumulative efficiency loss ΔA(t) = A[1] - A[t] (left y-axis) against
electricity market prices (right y-axis), reading A from JSON result files
and prices from electricity_markets.jl.

Run from repo root:
    julia --project=. plots/plot_efficiency_vs_price.jl
"""

include("electricity_markets.jl")   # loads spp_arr_DAM, spp_arr_DAM_pan, distrib_10

using Plots
using JSON
using XLSX, DataFrames, Dates

fig_dir = joinpath(@__DIR__, "figures/eff_price")
data_dir = joinpath(@__DIR__, "..", "data")

# ── Date index helper (shared with plot_year_slice.jl) ───────────────────────
# Load DeliveryDate column once; row order matches spp_arr_DAM / spp_arr_DAM_pan
const _dates_all = Date.(DataFrame(XLSX.readtable(
    joinpath(@__DIR__, "../data/ERCOT_DAM_AVG_2014-2024.xlsx"), "Sheet1"))[:, "DeliveryDate"])

"""
    year_range_indices(years) -> (i_start, i_end)

Return the first and last row indices (into spp_arr_DAM / A) that belong to
the given calendar `years` (e.g. `2022:2024`).
"""
function year_range_indices(years)
    mask = [Dates.year(d) in years for d in _dates_all]
    idxs = findall(mask)
    return first(idxs), last(idxs)
end

# ── Helper: load A and z_replace from a JSON result file ─────────────────────
function load_result(path::String)
    d = JSON.parsefile(path)
    fvr = d["final_value_results"]
    A         = Float64.(fvr["A"])
    z_replace = Float64.(fvr["z_replace"])
    return A, z_replace
end

function load_power(path::String)
    d = JSON.parsefile(path)
    fvr = d["final_value_results"]
    e         = Float64.(fvr["e_tot"])
    # z_replace = Float64.(fvr["z_replace"])
    return e
end

# ── Helper: compute cumulative efficiency loss ΔA(t) = A[1] - A[t] ──────────
cumulative_loss(A::Vector{Float64}) = (A[1] .- A) .* 3.336

# ── Helper: dual-axis plot of ΔA(t) vs price ─────────────────────────────────
"""
    plot_A_vs_price(A, price; kwargs...)

Plot cumulative efficiency loss ΔA(t) = A[1] - A[t] on the left y-axis and
electricity `price` on the right y-axis over a common time index. Replacement
events can optionally be marked by passing a `z_replace` vector.

Keyword arguments:
- `title`       – plot title
- `t_start`     – first hour index (default 1)
- `t_end`       – last hour index (default min of both lengths)
- `dates`       – Vector{Date} slice _dates_all[t_start:t_end]; when provided
                  the x-axis shows calendar years instead of hour indices
- `z_replace`   – binary replacement vector; events shown as vertical lines
- `price_label` – legend label for the price series
- `price_color` – colour for the price series (default :steelblue)
- `A_color`     – colour for the ΔA series (default :crimson)
"""
function plot_A_vs_price(A::Vector{Float64}, price::Vector{Float64};
        title        = "Cumulative Efficiency Loss vs. Electricity Price",
        t_start::Int = 1,
        t_end::Int   = min(length(A), length(price)),
        dates        = nothing,
        z_replace    = nothing,
        price_label  = "Electricity Price (\$/MWh)",
        price_color  = :steelblue,
        A_color      = :crimson)

    t  = t_start:t_end
    ΔA = cumulative_loss(A)

    # Build x-axis ticks from the date slice
    if !isnothing(dates)
        n_years = length(unique(Dates.year.(dates)))
        if n_years == 1
            # Single year: one tick per month (Jan–Dec)
            seen_months = Set{Int}()
            tick_pos    = Int[]
            tick_labs   = String[]
            month_names = ["Jan","Feb","Mar","Apr","May","Jun",
                           "Jul","Aug","Sep","Oct","Nov","Dec"]
            for (i, d) in enumerate(dates)
                m = Dates.month(d)
                if m ∉ seen_months
                    push!(seen_months, m)
                    push!(tick_pos,  t_start + i - 1)
                    push!(tick_labs, month_names[m])
                end
            end
            xlabel_str = string(Dates.year(first(dates)))
        else
            # Multi-year: one tick at Jan 1 of each calendar year
            seen_years = Set{Int}()
            tick_pos   = Int[]
            tick_labs  = String[]
            for (i, d) in enumerate(dates)
                y = Dates.year(d)
                if y ∉ seen_years
                    push!(seen_years, y)
                    push!(tick_pos,  t_start + i - 1)
                    push!(tick_labs, string(y))
                end
            end
            xlabel_str = "Year"
        end
        xticks_val = (tick_pos, tick_labs)
    else
        xticks_val = :auto
        xlabel_str = "Time (Hour)"
    end

    # Left axis: cumulative efficiency loss with fill below
    p = plot(t, ΔA[t];
        ylabel             = "Cumulative Efficiency Loss (% LHV)",
        xlabel             = xlabel_str,
        xticks             = xticks_val,
        title              = title,
        label              = false,
        color              = A_color,
        linewidth          = 1.2,
        legend             = false,
        fillrange          = 0,
        fillalpha          = 0.2,
        fillcolor          = A_color,
        yguidefontcolor    = A_color,
        right_margin       = 12Plots.mm,
        left_margin        = 8Plots.mm)

    # Mark replacement events on the left-axis series.
    # z_replace is now yearly (length n_years); year y's replacement is applied
    # at the boundary hour T_star[y] = y*8760 + 1.
    if !isnothing(z_replace)
        replace_hours = [y * 8760 + 1 for y in 1:length(z_replace) if z_replace[y] > 0.5]
        filter!(h -> t_start <= h <= t_end, replace_hours)
        if !isempty(replace_hours)
            vline!(p, replace_hours;
                label     = false,
                color     = :black,
                linestyle = :dash,
                linewidth = 0.8,
                alpha     = 0.7)
        end
    end

    # Right axis: price
    p2 = twinx(p)
    plot!(p2, t, price[t];
        ylabel          = "Electricity Price (\$/MWh)",
        label           = false,
        color           = price_color,
        linewidth       = 0.6,
        alpha           = 0.5,
        legend          = false,
        yguidefontcolor = price_color)

    return p
end

# ── 1. Full horizon: DAM average price vs A (basecase_avg) ───────────────────
A_avg, zr_avg = load_result(joinpath(data_dir, "ercot", "base_avg.JSON"))
# T_common_avg  = min(length(A_avg), length(spp_arr_DAM))

# p1 = plot_A_vs_price(A_avg, spp_arr_DAM;
#     title       = "Cumulative Efficiency Loss vs. DAM Average Price",
#     t_end       = T_common_avg,
#     dates       = _dates_all[1:T_common_avg],
#     z_replace   = zr_avg,
#     price_label = "DAM Average (\$/MWh)")
# savefig(p1, joinpath(fig_dir, "efficiency_vs_dam_avg.png"))

# ── 2. Full horizon: DAM panhandle price vs A (basecase_pan) ─────────────────
# A_pan, zr_pan = load_result(joinpath(data_dir, "test", "basecase_pan.JSON"))
# T_common_pan  = min(length(A_pan), length(spp_arr_DAM_pan))

# p2 = plot_A_vs_price(A_pan, spp_arr_DAM_pan;
#     title       = "Cumulative Efficiency Loss vs. DAM Panhandle Price",
#     t_end       = T_common_pan,
#     dates       = _dates_all[1:T_common_pan],
#     z_replace   = zr_pan,
#     price_label = "DAM Panhandle (\$/MWh)",
#     price_color = :darkorange)
# savefig(p2, joinpath(fig_dir, "efficiency_vs_dam_pan.png"))

# ── 3. Year 1 zoom: DAM average ───────────────────────────────────────────────
# p3 = plot_A_vs_price(A_avg, spp_arr_DAM;
#     title       = "Cumulative Efficiency Loss vs. DAM Average Price – Year 1",
#     t_end       = 8760,
#     dates       = _dates_all[1:8760],
#     z_replace   = zr_avg,
#     price_label = "DAM Average (\$/MWh)")
# savefig(p3, joinpath(fig_dir, "efficiency_vs_dam_avg_yr1.png"))

# ── 4. Year 1 zoom: Panhandle ─────────────────────────────────────────────────
# p4 = plot_A_vs_price(A_pan, spp_arr_DAM_pan;
#     title       = "Cumulative Efficiency Loss vs. DAM Panhandle Price – Year 1",
#     t_end       = 8760,
#     dates       = _dates_all[1:8760],
#     z_replace   = zr_pan,
#     price_label = "DAM Panhandle (\$/MWh)",
#     price_color = :darkorange)
# savefig(p4, joinpath(fig_dir, "efficiency_vs_dam_pan_yr1.png"))

# ── 5. Side-by-side: avg vs panhandle full horizon ────────────────────────────
# p5 = plot(p1, p2; layout=(1, 2), size=(1600, 500))
# savefig(p5, joinpath(fig_dir, "efficiency_vs_price_both.png"))

# ── 2022 only: DAM hub average ───────────────────────────────────────────────
y22_start, y22_end = year_range_indices(2022:2022)

# p_2022 = plot_A_vs_price(A_avg, spp_arr_DAM;
#     title       = "Cumulative Efficiency Loss vs. DAM Average Price (2022)",
#     t_start     = y22_start,
#     t_end       = y22_end,
#     dates       = _dates_all[y22_start:y22_end],
#     z_replace   = zr_avg,
#     price_label = "DAM Average (\$/MWh)")
# savefig(p_2022, joinpath(fig_dir, "efficiency_vs_dam_avg_2022.png"))

# ── 2022–2024 slice ───────────────────────────────────────────────────────────
# ys_start, ys_end = year_range_indices(2022:2024)

# ── 6. 2022–2024: DAM average price vs A ─────────────────────────────────────
# p6 = plot_A_vs_price(A_avg, spp_arr_DAM;
#     title       = "Cumulative Efficiency Loss vs. DAM Average Price (2022–2024)",
#     t_start     = ys_start,
#     t_end       = ys_end,
#     dates       = _dates_all[ys_start:ys_end],
#     z_replace   = zr_avg,
#     price_label = "DAM Average (\$/MWh)")
# savefig(p6, joinpath(fig_dir, "efficiency_vs_dam_avg_2022_2024.png"))

# ── 7. 2022–2024: DAM panhandle price vs A ───────────────────────────────────
# p7 = plot_A_vs_price(A_pan, spp_arr_DAM_pan;
#     title       = "Cumulative Efficiency Loss vs. DAM Panhandle Price (2022–2024)",
#     t_start     = ys_start,
#     t_end       = ys_end,
#     dates       = _dates_all[ys_start:ys_end],
#     z_replace   = zr_pan,
#     price_label = "DAM Panhandle (\$/MWh)",
#     price_color = :darkorange)
# savefig(p7, joinpath(fig_dir, "efficiency_vs_dam_pan_2022_2024.png"))

# ── 8. Side-by-side: 2022–2024 avg vs panhandle ───────────────────────────────
# p8 = plot(p6, p7; layout=(1, 2), size=(1600, 500))
# savefig(p8, joinpath(fig_dir, "efficiency_vs_price_both_2022_2024.png"))

# ── Helper: overlay cumulative efficiency loss for two operation cases ─────────
"""
    plot_A_compare(A_flex, A_const, price; kwargs...)

Overlay cumulative efficiency loss for a flexible and a constant (non-flexible)
operation case on the left y-axis, with an electricity price series on the right
y-axis. Both ΔA curves are zeroed at t_start so they share a common baseline.

Keyword arguments:
- `title`        – plot title
- `t_start`      – first hour index (default 1)
- `t_end`        – last hour index (default min of all lengths)
- `dates`        – Vector{Date} slice for calendar x-axis ticks
- `flex_color`   – colour for the flexible series     (default :crimson)
- `const_color`  – colour for the constant series     (default :steelblue)
- `price_color`  – colour for the price series        (default :darkorange)
- `flex_label`   – legend label for flexible case     (default "Flexible")
- `const_label`  – legend label for constant case     (default "Non-flexible")
"""
function plot_A_compare(A_flex::Vector{Float64}, A_const::Vector{Float64},
                        price::Vector{Float64};
        title        = "",
        t_start::Int = 1,
        t_end::Int   = min(length(A_flex), length(A_const), length(price)),
        dates        = nothing,
        flex_color   = :crimson,
        const_color  = :steelblue,
        price_color  = :darkorange,
        flex_label   = "Flexible",
        const_label  = "Non-flexible")

    t = t_start:t_end

    # Zero both curves at t_start so they share a common baseline
    ΔA_full_flex  = cumulative_loss(A_flex)
    ΔA_full_const = cumulative_loss(A_const)
    ΔA_flex_slice  = ΔA_full_flex[t]  .- ΔA_full_flex[t_start]
    ΔA_const_slice = ΔA_full_const[t] .- ΔA_full_const[t_start]

    # Build x-axis ticks from the date slice
    if !isnothing(dates)
        n_years = length(unique(Dates.year.(dates)))
        if n_years == 1
            seen_months = Set{Int}()
            tick_pos    = Int[]
            tick_labs   = String[]
            month_names = ["Jan","Feb","Mar","Apr","May","Jun",
                           "Jul","Aug","Sep","Oct","Nov","Dec"]
            for (i, d) in enumerate(dates)
                m = Dates.month(d)
                if m ∉ seen_months
                    push!(seen_months, m)
                    push!(tick_pos,  t_start + i - 1)
                    push!(tick_labs, month_names[m])
                end
            end
            xlabel_str = string(Dates.year(first(dates)))
        else
            seen_years = Set{Int}()
            tick_pos   = Int[]
            tick_labs  = String[]
            for (i, d) in enumerate(dates)
                y = Dates.year(d)
                if y ∉ seen_years
                    push!(seen_years, y)
                    push!(tick_pos,  t_start + i - 1)
                    push!(tick_labs, string(y))
                end
            end
            xlabel_str = "Year"
        end
        xticks_val = (tick_pos, tick_labs)
    else
        xticks_val = :auto
        xlabel_str = "Time (Hour)"
    end

    # Left axis: both efficiency loss curves
    p = plot(t, ΔA_flex_slice;
        ylabel       = "Cumulative Efficiency Loss (% LHV)",
        xlabel       = xlabel_str,
        xticks       = xticks_val,
        title        = title,
        label        = flex_label,
        color        = flex_color,
        linewidth    = 1.2,
        fillrange    = 0,
        fillalpha    = 0.15,
        fillcolor    = flex_color,
        left_margin  = 10Plots.mm,
        right_margin = 12Plots.mm, 
        size= (800,600))

    plot!(p, t, ΔA_const_slice;
        label     = const_label,
        color     = const_color,
        linewidth = 1.2,
        fillrange = 0,
        fillalpha = 0.15,
        fillcolor = const_color)

    # Right axis: electricity price
    p2 = twinx(p)
    plot!(p2, t, price[t];
        ylabel          = "Electricity Price (\$/MWh)",
        label           = false,
        color           = price_color,
        linewidth       = 0.6,
        alpha           = 0.5,
        legend          = false,
        yguidefontcolor = price_color)

    return p
end

# ── 9. 2022: flexible vs. non-flexible degradation (DAM hub average) ──────────
A_avg_const, _ = load_result(joinpath(data_dir, "ercot", "base_avg_constant.JSON"))

p9 = plot_A_compare(A_avg, A_avg_const, spp_arr_DAM;
    # title       = "Cumulative Efficiency Loss: Flexible vs. Non-Flexible (2022)",
    t_start     = y22_start,
    t_end       = y22_end,
    dates       = _dates_all[y22_start:y22_end])
savefig(p9, joinpath(fig_dir, "efficiency_flex_vs_const_2022.svg"))

println("All efficiency-vs-price plots saved to plots/figures/")


"""
    plot_e_compare(e_flex, e_const, price; kwargs...)

Overlay power consumption for a flexible and a non-flexible operation case on
the left y-axis, with an electricity price series on the right y-axis.

For weekly windows (≤ 168 hours), x-axis ticks are placed at each day boundary
and labelled with the day name (Mon, Tue, …). For longer windows, month or year
ticks are used depending on the span.

Keyword arguments:
- `title`        – plot title
- `t_start`      – first hour index (default 1)
- `t_end`        – last hour index (default min of all lengths)
- `dates`        – Vector{Date} slice aligned to the hour range
- `flex_color`   – colour for the flexible series  (default :crimson)
- `const_color`  – colour for the constant series  (default :steelblue)
- `price_color`  – colour for the price series     (default :darkorange)
- `flex_label`   – legend label for flexible case  (default "Flexible")
- `const_label`  – legend label for constant case  (default "Non-flexible")
"""
function plot_e_compare(e_flex::Vector{Float64}, e_const::Vector{Float64},
                        price::Vector{Float64};
        title        = "",
        t_start::Int = 1,
        t_end::Int   = min(length(e_flex), length(e_const), length(price)),
        dates        = nothing,
        flex_color   = :crimson,
        const_color  = :steelblue,
        price_color  = :darkorange,
        flex_label   = "Flexible",
        const_label  = "Non-flexible")

    t = t_start:t_end
    n_hours = length(t)

    # Build x-axis ticks
    if !isnothing(dates)
        if n_hours <= 168
            # Weekly view: one tick per day, regardless of calendar alignment.
            # Branch on window length rather than unique-date count so that
            # windows not starting at midnight (≤8 unique dates) still work.
            day_names = ["Mon","Tue","Wed","Thu","Fri","Sat","Sun"]
            seen_days = Set{Date}()
            tick_pos  = Int[]
            tick_labs = String[]
            for (i, d) in enumerate(dates)
                if d ∉ seen_days
                    push!(seen_days, d)
                    push!(tick_pos,  t_start + i - 1)
                    push!(tick_labs, day_names[Dates.dayofweek(d)])
                end
            end
            xlabel_str = string(first(dates)) * " – " * string(last(dates))
        elseif length(unique(Dates.year.(dates))) == 1
            # Single year: one tick per month
            seen_months = Set{Int}()
            tick_pos    = Int[]
            tick_labs   = String[]
            month_names = ["Jan","Feb","Mar","Apr","May","Jun",
                           "Jul","Aug","Sep","Oct","Nov","Dec"]
            for (i, d) in enumerate(dates)
                m = Dates.month(d)
                if m ∉ seen_months
                    push!(seen_months, m)
                    push!(tick_pos,  t_start + i - 1)
                    push!(tick_labs, month_names[m])
                end
            end
            xlabel_str = string(Dates.year(first(dates)))
        else
            # Multi-year: one tick per year
            seen_years = Set{Int}()
            tick_pos   = Int[]
            tick_labs  = String[]
            for (i, d) in enumerate(dates)
                y = Dates.year(d)
                if y ∉ seen_years
                    push!(seen_years, y)
                    push!(tick_pos,  t_start + i - 1)
                    push!(tick_labs, string(y))
                end
            end
            xlabel_str = "Year"
        end
        xticks_val = (tick_pos, tick_labs)
    else
        xticks_val = :auto
        xlabel_str = "Time (Hour)"
    end

    # Left axis: power consumption (MW)
    p = plot(t, e_flex[t];
        ylabel       = "Power Consumption (MW)",
        xlabel       = xlabel_str,
        xticks       = xticks_val,
        title        = title,
        legend = false,
        label        = flex_label,
        yguidefontcolor = flex_color,
        color        = flex_color,
        linewidth    = 1.2,
        fillrange    = 0,
        fillalpha    = 0.15,
        fillcolor    = flex_color,
        left_margin  = 10Plots.mm,
        right_margin = 12Plots.mm,
        bottom_margin = 10Plots.mm,
        size         = (900, 600))

    # plot!(p, t, e_const[t];
    #     label     = const_label,
    #     color     = const_color,
    #     linewidth = 1.2,
    #     fillrange = 0,
    #     fillalpha = 0.15,
    #     fillcolor = const_color)
    # plot!(p, t, e_const[t];
    #     label     = const_label,
    #     color     = const_color,
    #     linewidth = 1.2,)

    # Right axis: electricity price
    p2 = twinx(p)
    plot!(p2, t, price[t];
        ylabel          = "Electricity Price (\$/MWh)",
        label           = false,
        color           = price_color,
        linewidth       = 1.0,
        alpha           = 0.6,
        legend          = false,
        yguidefontcolor = price_color)

    return p
end

# ── Load power consumption series ────────────────────────────────────────────
e_avg_flex  = load_power(joinpath(data_dir, "ercot", "base_avg.JSON"))
e_avg_const = load_power(joinpath(data_dir, "ercot", "base_avg_constant.JSON"))

# ── 10. Representative week: flexible vs. non-flexible power consumption ──────
# Use the first full week of 2022 (168 hours starting at y22_start)
(y24_start, y24_end) = year_range_indices(2023:2023)
# wk_start = y24_end -167 -1000
# wk_end   = y24_end -1000
i_start = 82588
wk_start = i_start
wk_end = i_start + 167

p10 = plot_e_compare(e_avg_flex, e_avg_const, spp_arr_DAM;
    title   = "",
    t_start = wk_start,
    t_end   = wk_end,
    dates   = _dates_all[wk_start:wk_end],
    price_color = :black
    )
savefig(p10, joinpath(fig_dir, "power_flex_vs_const_week.pdf"))

