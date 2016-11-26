# data for CHP problem
# Yankai Cao, Siyu Chen, Luis Fuentes
# UW-Madison, 2016

GD=readdlm("GDdata.csv",',')
WD=readdlm("WDdata.csv",',')
TAMB=readdlm("Temperature.csv",',')
CostE=readdlm("CostE.dat",',')

DT=1/12  # data frequency per hour(5 mins)
T=288    # number of time steps per day
TST0=100 # temperature of storage tank at time of 0
CkWe=400 # fixed cost for CHP system (kWe)
CST=125  # fixed cost for storage tank (stere)
kF=0.23  # depreciation rate (years)
PT=1/365 # scaling factor

CEm=40               # cost coefficient
fBio=0.000044*DT     # emission factor of biogas
fNG=0.00038*DT       # emission factor of natural gas
efWCHP=0.3725        # maximum system efficiency
efQCHP=0.475         # thermal efficiency
CostNG=0.01326296*DT # cost of natual gas
CostBio=0.034095*DT  # cost of biogas
CpW=0.00116/DT       # heat capacity of water
U=0.001*DT           # heat transfer coefficient
rho=1000             # density of water
TCHP=100             # discharge temperature of CHP system
TS=50                # temperature of supply water
CostW=0.00069        # cost of water
PSE=0.085*DT         # sales price of electricity
PHW=0.016*DT         # sales price of hot water
CMCHP=0.015*DT       # operational cost of CHP system