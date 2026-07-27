using DataFrames
using XLSX
using Dates
using Statistics

if !@isdefined(DATA_DIR)
    const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
end

# (filename, sheet name) for each 2024 monthly AER report, in calendar order
const MONTH_FILES_2024 = [
    ("aer-2024-01-01.xlsx", "AER_SUMMARY"),
    ("aer-2024-02-01.xlsx", "SUM_ALL"),
    ("aer-2024-03-01_v2.xlsx", "AER_SUMMARY"),
    ("average-emissions-rate-report-april-2024.xlsx", "AER_SUMMARY"),
    ("average-emissions-rate-report-may-2024.xlsx", "AER_SUMMARY"),
    ("average-emissions-rate-report-june-2024.xlsx", "AER_SUMMARY"),
    ("average-emissions-rate-report-july-2024.xlsx", "AER_SUMMARY"),
    ("average-emissions-rate-report-august-2024.xlsx", "AER_SUMMARY"),
    ("average-emissions-rate-report-september-2024.xlsx", "AER_SUMMARY"),
    ("average-emissions-rate-report-october-2024.xlsx", "AER_SUMMARY"),
    ("average-emissions-rate-report-november-2024.xlsx", "AER_SUMMARY"),
    ("average-emissions-rate-report-december-2024.xlsx", "AER_SUMMARY"),
]

const MONTH_NAMES = [
    "january", "february", "march", "april", "may", "june",
    "july", "august", "september", "october", "november", "december",
]

# 2025 monthly reports all follow the same naming/sheet convention
month_files_2025() = [
    ("average-emissions-rate-report-$(m)-2025.xlsx", "AER_SUMMARY") for m in MONTH_NAMES
]

"""
Load and concatenate a year's worth of monthly AER reports from data/<year>/.
Returns `nothing` (with a warning) if any monthly file is missing.
"""
function load_monthly_year(year::Int, files)
    year_dir = joinpath(DATA_DIR, string(year))
    dfs = DataFrame[]
    for (fname, sheet) in files
        path = joinpath(year_dir, fname)
        if !isfile(path)
            @warn "Skipping $year emissions data: missing $path"
            return nothing
        end
        push!(dfs, DataFrame(XLSX.readtable(path, sheet)))
    end
    return vcat(dfs...)
end

function parse_trade_date(x)
    if x isa Date
        return x
    end

    s = String(x)

    for fmt in (dateformat"yyyy-mm-dd", dateformat"m/d/y", dateformat"m/d/yyyy")
        try
            return Date(s, fmt)
        catch
        end
    end

    error("Could not parse date: $s")
end

function clean_emissions(df)
    vals = Float64[]

    df2 = df[.!ismissing.(df.TRADE_DT), :]

    df2.TRADE_DT = [parse_trade_date(x) for x in df2.TRADE_DT]

    for d in sort(unique(df2.TRADE_DT))


        day_tf = df2[df2.TRADE_DT .== d, :]

        for h in 1:24
            rows = day_tf[day_tf.TRADE_HR .== h, :]

            if nrow(rows) == 1 && !ismissing(rows.AVG_EM_RATE[1])
                push!(vals, Float64(rows.AVG_EM_RATE[1]))

            elseif nrow(rows) > 1
                push!(vals, mean(skipmissing(rows.AVG_EM_RATE)))

            else
                before = day_tf[day_tf.TRADE_HR .== h - 1, :AVG_EM_RATE]
                after  = day_tf[day_tf.TRADE_HR .== h + 1, :AVG_EM_RATE]

                neigh = collect(skipmissing(vcat(before, after)))

                push!(vals,
                    isempty(neigh) ? vals[end] : mean(neigh)
                )
            end
        end
    end

    @assert length(vals) >= 8760 "Only got $(length(vals)) hourly values"

    return vals[1:8760]
end

# Reads a single wide AER_SUMMARY sheet (2022/2023 report layout: one file per year,
# fixed A1:K8761 range, header row may contain blank/missing columns).
function read_wide_aer_sheet(path)
    xf = XLSX.readxlsx(path)
    sheet = xf["AER_SUMMARY"]

    data = sheet["A1:K8761"]

    raw_headers = vec(data[1, :])
    headers = Symbol.([
        ismissing(h) ? "blank_$i" : String(h)
        for (i, h) in enumerate(raw_headers)])

    vals = data[2:end, :]
    return DataFrame(vals, headers)
end

function load_2022()
    path = joinpath(DATA_DIR, "2022", "2022_averageemissionsrate.xlsx")
    if !isfile(path)
        @warn "Skipping 2022 emissions data: missing $path"
        return nothing
    end
    return clean_emissions(read_wide_aer_sheet(path))
end

function load_2023()
    path = joinpath(DATA_DIR, "2023", "aer_summary_2023_full_year_w_ca_attr_final.xlsx")
    if !isfile(path)
        @warn "Skipping 2023 emissions data: missing $path"
        return nothing
    end

    y2023 = read_wide_aer_sheet(path)

    df2 = y2023[.!ismissing.(y2023.TRADE_DT), :]

    df2.TRADE_DT = [
        x isa Date ? x : Date(String(x))
        for x in df2.TRADE_DT
    ]

    # March 11, 2023 is missing from the source report; interpolate it as the
    # average of the surrounding days (DST fall-back leaves a 23-hour gap otherwise).
    march10 = df2[df2.TRADE_DT .== Date(2023, 3, 10), :]
    march12 = df2[df2.TRADE_DT .== Date(2023, 3, 12), :]

    march11 = deepcopy(march12)
    march11.TRADE_DT .= Date(2023, 3, 11)

    for h in 1:24
        m10 = march10[march10.TRADE_HR .== h, :AVG_EM_RATE]
        m12 = march12[march12.TRADE_HR .== h, :AVG_EM_RATE]

        if !isempty(m10) && !isempty(m12)
            march11[march11.TRADE_HR .== h, :AVG_EM_RATE] .= (m10[1] + m12[1]) / 2
        end
    end

    y2023 = sort(vcat(df2, march11), [:TRADE_DT, :TRADE_HR])

    return clean_emissions(y2023)
end

"""
    load_emissions() -> Dict{Int, Vector{Float64}}

Loads whatever years of ERCOT average-emissions-rate data are available under `data/`.
Years with missing source files are skipped with a `@warn` rather than erroring, so
this only fails if no year could be loaded at all.
"""
function load_emissions()
    years_data = Dict{Int, Vector{Float64}}()

    em2022 = load_2022()
    em2022 !== nothing && (years_data[2022] = em2022)

    em2023 = load_2023()
    em2023 !== nothing && (years_data[2023] = em2023)

    df2024 = load_monthly_year(2024, MONTH_FILES_2024)
    df2024 !== nothing && (years_data[2024] = clean_emissions(df2024))

    df2025 = load_monthly_year(2025, month_files_2025())
    df2025 !== nothing && (years_data[2025] = clean_emissions(df2025))

    isempty(years_data) && error(
        "No emissions data could be loaded from $DATA_DIR — populate data/<year>/ with AER reports."
    )

    return years_data
end
