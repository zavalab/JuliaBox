Example 1:
  Run two_objectives_v2.jl in Julia
  Use 'flagAltNadir' to toggle between original and alternate nadir point
  
Example 2:
  Run three_objectives2.jl in Julia
  Use 'flagAltNadir' to toggle between original and alternate nadir point
  Use stakeholder_hist.jl to visualize distribution of stakeholder dissatisfactions
  
Case Study:
  Run facility_location_problem3.jl in Julia. This script with create text files
  Use 'flagAltNadir' to toggle between original and alternate nadir point
  Use 'nStake = ' to adjust the number of stakeholders
  Run plot_facility_results in Python to create plots. 'iSaveAlt = ' and 'iSaveTrad = ' may need to be updated if the problem input data is changes (and the values of alpha for which the solution changes also changes)
  Use stakeholder_weights.jl to visualize distribution of stakeholder weights

Note:
  Julia Version 0.4.4-pre+24 was for the examples in the CVaR manuscript