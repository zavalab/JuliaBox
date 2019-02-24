Vbin = 3:25                               # mean wind speed set
seedSet = 0:9                             # seed set for each mean wind speed
nseed = length(seedSet)                   # number of seed per bin
S = length(Vbin)*length(seedSet)          # number of scenarios
SCEN = 1:S                                # scenario set

Vave = 10
deltaV = 1
threshold = 6                            # load threshold [1e7 N m]
Euler = 0.577216
#probability of each mean wind speed
deltaP = exp(-pi*  ((Vbin-deltaV/2)/2/Vave).^2) - exp(-pi*((Vbin+deltaV/2)/2/Vave).^2)
alpha = 0.5;           			  # probability level
firsttime = 10                                                                             #seconds
firstNt = round(Int, firsttime/dt)+1
firstTIMEG = 1:(firstNt)                                                                   #set of temporal grid points