using ColorSchemes
using Gurobi
using Statistics
using XLSX

##############################################
# Initialize Cases
##############################################

time_limit = 216000                 # seconds
af = 0.117
afcap_cost = 6.014735319879001e9    # from MIP solution to Basis (180 hr run)
cap_cost = afcap_cost / af          # required for full processing

# budgets = ["Basis", "Cap 90%", "Cap 75%", "Cap 50%", "Cap 25%"]   # all case study options
budgets = ["Basis"]

##############################################
# Begin Looping
##############################################

base_path = dirname(Base.source_path())

include("National_Upcycling_MIP_Function.jl")

for i=1:length(budgets)

    path = string(base_path, "/", budgets[i])

    if !isdir(path)
        mkdir(path)
    end

    National_Upcycling_MIP(path, (budgets[i] * cap_cost), time_limit)

end