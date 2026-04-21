include("electrolyzer_struct.jl")
include("save_optimization.jl")
include("formulation.jl")
using CPUTime

# Purpose of this is to just run the sensitivity analysis of the MS-elec work

# Import the data from the ERCOT ISO
DAM_10yr = DataFrame(XLSX.readtable((@__DIR__)*"/../data/ERCOT_DAM_AVG_2014-2024.xlsx", "Sheet1"))
spp_arr_DAM = Float64.(DAM_10yr[:, "SettlementPointPrice"])
DAM_20yr_avg = vcat(spp_arr_DAM, spp_arr_DAM)

efficiencies = [0.60]
# efficiencies = [0.65, 0.70, 0.75]
# op_degrads = [0, 1.6, 3.2, 4.8, 6.4, 8.0, 10.0]
op_degrads = [1.6, 10.0]
 
for e in efficiencies 
    for d in op_degrads
        try
            # 1. Define the electrolyzer
            el = Electrolyzer(2.2, e, d/1000000, 1.9, 10, 0.0204/500, 3, .05, 1816*1000, 250*1000, 0.05)
            
            # 2. Run optimization
            # t_h_sb = @CPUelapsed h_sb = run_3st_opt(θ_h_sb, DAM_20yr_avg)
            t_e_d = @CPUelapsed m_e_d = run_3st_opt(el, spp_arr_DAM, gap=0.03)

            # 3. Check solver status (Crucial for "Intractable" results)
            # Replace 'termination_status' with the appropriate check for your specific 'run_3st_opt' output
            status = termination_status(m_e_d)
            
            if status == MOI.OPTIMAL
                filename = "../data/heatmap/eff_$(e)_deg_$(d)_results.JSON"
                # pop_and_save(m_e_d, el, spp_arr_DAM, filename)

                pop_and_save(m_e_d, t_e_d ,el, filename)
                println("Success: e=$e, d=$d")
            else
                println("Warning: Optimization failed for e=$e, d=$d with status: $status")
            end

        catch err
            # This catches code crashes, memory issues, or solver-level exceptions
            @warn "Skip: e=$e, d=$d produced an error" exception=(err, catch_backtrace())
            continue # Move to the next iteration
        end
    end
end
