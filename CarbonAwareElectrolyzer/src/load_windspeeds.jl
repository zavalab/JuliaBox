#read in wind data
using DataFrames
using XLSX
using Dates
using Statistics

if !@isdefined(DATA_DIR)
    const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
end

function load_windspeeds()

    path = joinpath(DATA_DIR, "Data_weather.xlsx")

    xf = XLSX.readxlsx(path)
    sheet = xf["2024"]

    data = sheet["A1:L8761"]

    raw_headers = vec(data[1, :])
    headers = Symbol.([
        ismissing(h) ? "blank_$i" : String(h)
        for (i, h) in enumerate(raw_headers)])

    vals = data[2:end, :]
    data_2024= DataFrame(vals, headers)
    wind_data=data_2024[1:8760,7]

    return wind_data
end



