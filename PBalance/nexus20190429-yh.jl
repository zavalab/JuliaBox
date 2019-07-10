#=
Balance of phosphorous in the dairy nexus
Julia 0.6.4
Created by Dulce Lopez Diaz, January 2018
Modified by Yicheng Hu, Apr 2019
=#


#=Sectors
ES-environment
MIS-mining sector
CS-chemical sector
S-cultivated soil
RC-crop production
HC-harvested crop (crop processing)
F-food distribution sector
US-urban sector
WS-wastewater treatment sector
M -manure management sector
DC-Dairy Castle sector
DL-Dairy waste sector (deleted this one, combined with ES)
DS-Dairy production
B-Dairy sector losses (? balance stream)
E-Exports
=#

# Define parameters for the mass balance

##Mining sector
alphaMIS=0.045 #P removal efficiency (tonne loss P/tonne extracted P from phosphate)
#alphaMIS=0.0315


##Chemical Sector
alphaCS=0.09      # fraction of loss P in the chemical sector
#alphaCS=0.045
betaCSDC=2.5e-6  # requirements of additives by cow population (tonne P/cow/year)
#*** in P mass balance the additives are not significant comparating with other flows but in money balances the additives are important.

##Agriculture Sector
alphaSES=0.18    #fraction that represents the losses in the soil stage
#alphaSES=0.126
alphaRCS=0.15      # fraction biochar recycle to cultivated soil (conversion factor in the biochar process)
alphaRCES=0.08     # fraction losses of crops activities
#alphaRCES=0.04
alphaHCES=0.13     # loses from crops
#alphaHCES=0.075
betapc=15.966451*0.001 # corn consumption rate (tonne/person/year) fron USDA https://www.ers.usda.gov/data-products/chart-gallery/gallery/chart-detail/?chartId=85059

##Food Distribution Sector
alphaFES=0.4      # fraction for organic residues in food commodaties (losses)

##Urban Sector
alphaUSES=0.18     # fraction of losses by the human sector

##WWTP Sector
alphaWSS=0.75      # yield of residues that are used as fertilizers

##Dairy Sector
xcorn=3.5e-3      # P content in the corn (tonne P/ton corn)
betarc=3.45       # requirements of corn by the cows (tonne corn/cow/year)
alphaDCES=0.05    # fraction of dairy cow losses
betamc=0.9        # fraction of milk that is transformed to cheese
xcheese=7.809e-3  # P content in cheese (tonne P/tonne cheese)
betaprc=1.002e-4  # cheese yield (tonne cheese/L milk)
xmilk=9.578e-7    # P content in milk (tonne P/L milk)
betamilk=10367    # annual production of milk per cow (L milk/cow/year)
betamp=67.682     # milk consumption rate (L milk/person/year)
betachp=0.0174    # cheese consumption rate (tonne cheese/person/year)


##Manure Management Sector
alphaMES= 0.18    # factor of losses in the production of manure for transportation and disposal treatments


# Define parameters for economic analysis

##Mining Sector
crock=80          # price per ton of p rock (USD/tonne phosphate rock)
xrock=0.123        # P content in phosphate rock (tonne P/tonne phosphate rock) 11200/163000*(31*2/(31*2+5*16))
cwaste=0
#cwaste=99119     # cost to removal p from the waste

##Chemical Sector
cadd=25      # price of additives (USD/tonne addictives)
cfert=350        # price of fertilizer (USD/tonne fertilizer)
xadd=0.04        # P content in the addictives (tonne P/tonne addictives)
##TF JUBB 1988 by taking average
xfert=0.092       # P content in the DAP fertilzer (tonne P/tonne DAP)

##Agriculture Sector
ccorn=148.52      # price of corn (USD/tonne corn)
csilage=52        # price of silage (USD/tonne silage)
Xcorn=0.321       # fraction of corn in agricultural products
Xsilage=0.679     # fraction of silage in agricultural products
xsilage=2e-3      # P content in silage (tonne P/ton silage)
cmanuret=7.85     # price of treated manure (USD/tonne treated manure)
csluge=15         # price of sluge (USD/tonne sludge)
xsludge=0.064     # P content in sludge (tonne P/tonne sludge) from Maddie's report
xmanuret=1.13e-3    # P content in manure (tonne P/tonne manure) from price paper
betars=betarc/Xcorn*Xsilage  # requirements of silage from cows (tonne silage/cow/year)

##Urban Sector
cdesh=0.1567      # price of wastewater gate fee (USD/tonne wastewater)
xdesh=1.455e-4     # P content in wastewater inflow (tonne P/tonne wastewater) from Maddie's report

##Dairy Sector
xmanurer=1.3e-3   # P content in manure (tonne P/tonne manure) from price paper
cmanurer=7.56     # price of raw manure (USD/tonne raw manure)
Xmilk=0.1         # fraction of P appearing in milk
Xcheese=0.9       # fraction of P appearing in cheese
ccheese=3083      # price of cheese (USD/tonne cheese)
cmilk=0.386       # price of milk (USD/L milk)


# Some unused parameters, please double check this part
# crop yields
Yc=12.241         # yield corn production ton corn grain/ha
Yf=45.679         # yield forage production ton forage/ha

#parameters for harvested crops
betacd=0.223      # percentage of corn that is required by the dairy sector
betafd=0.463      # percentage of forage that is required by the dairy sector

# parameter for manure production
betaDCM=1.3e-3    # factor manuere production ton manure/cow head

avail=70000       # P reserves
betacfood=0.0154  # dairy products requirements per person by the urban sector

cdiaprod=20.64    #price for dairy products
#ccheese=2.94     #cheese milk
#cdairylos=1.663  #dairy cost losses  * COST OF DISCHARGE WASTE EFFLUENT IN THE DAIRY SECTOR (cost of treatment1)
cdairylos=0.5     #dairy cost losses, cost of treatment 2
cbiochar=2580     #price of the biochar as a phosphorus source


# Modeling

using JuMP
using Ipopt
using Gurobi

m = Model(solver=GurobiSolver());
#@variable(m, crock>=34500)
#@variable(m, 0<=cmilk<=100)

#@variable(m, cwaste>=0)
#@variable(m, cmilk>=0)
#@variable(m, ccheese>=0)
#@variable(m, ccorn)
#@varibale(m, csilage)
#@variable(m, crock>=0)
#@variable(m, cfert>=0)
#Global variables
pcow=1280000 # cow population, heads
phuman=5819000  #human population
fDSE=9455 # p exported in dairy products
fMIS=38000*1.3  #Phosphorus flowrate in the mining sector

#Mining sector(MS)
# mass balance part
@variable(m,fMISES)
@variable(m,fMISCS)
@constraint(m, fMISCS==fMIS-fMISES)
@constraint(m, fMISES==alphaMIS*fMIS)

# value balance part
@variable(m,cfMISES)                     #cost of mining treatment
#=is this wrong?
@constraint(m, cfMISES==crock*fMISES)
=#
@constraint(m, cfMISES==cwaste*fMISES)
@variable(m,cfMISCS)                     #market price of phosphate rock
@constraint(m, cfMISCS==crock*fMISCS/xrock)

@variable(m, BMIS)
@constraint(m, BMIS == cfMISCS - cfMISES)

#Chemical Sector(CS)
# mass balance part
@variable(m,fCSES)      # P disfribution losses
@variable(m,fCSS)          # P fertilizer from chemical sector
@variable(m,fCSDC)         # P aditives to dairy cows
@constraint(m, fCSS==fMISCS-fCSES-fCSDC)     # mass balance in the chemical sector
@constraint(m, fCSES==alphaCS*fMISCS)         # distribution and processing losses in the chemical sector
@constraint(m, fCSDC==betaCSDC*pcow)            # p from the chemical sector to dairy sector as additives

# value balance part
@variable(m,cfCSS)         # cost in the system by fertilizers
@variable(m,cfCSES)        # cost to removal p in the waste
@variable(m,cfCSDC)      # additive costs
@constraint(m, cfCSS==1/xfert*cfert*fCSS)            # cost of the phosphorus in the fertilizer
@constraint(m, cfCSES==cwaste*fCSES)         # cost of the phoshorus that is discharged to the enviroment by the chemical sector
@constraint(m, cfCSDC==1/xadd*cadd*fCSDC)      #cost of additives that are required by the dairy sector
@variable(m, BCS)
@constraint(m, BCS == cfCSDC + cfCSS - cfMISCS - cfCSES)

#Agriculture Sector(AS)
# mass balance part
# soil balance (S)
@variable(m,fSES)          # p erossion losses
@variable(m,fSRC)          # p crop content
@variable(m,fRCS)          # biochar
@variable(m,fWSS)          # p flow from sludge
@variable(m,fMS)           # p flow from manure

@constraint(m, fSRC==fCSS+fMS+fWSS-fSES+fRCS)      # mass balance in the cultivated soil sector
@constraint(m, fSES==alphaSES*fCSS)                 # erosion losses
@constraint(m, fRCS==alphaRCS*fSRC)                 # biochar recycle

# crop sector (RC)
@variable(m,fRCES)        # p crop losses
@variable(m,fRCHC)        # p crop production to sites of distribution-harvested crop sector
@constraint(m, fRCHC==fSRC-fRCS-fRCES)     # mass balance in the crop production sector
@constraint(m, fRCES==alphaRCES*fSRC)       # crop losses

# harvested crops (HC)
@variable(m,fHCF)        # p flow to food commodoties
@variable(m,fHCES)    # p harvested losses
@variable(m,fHCDC)       # p crop products to dairy sector
@variable(m,fHCB)        # a balance P stream indicating export/import

@constraint(m, fHCF==fRCHC-fHCDC-fHCES-fHCB)     # balance in the harvested crop sector
@constraint(m, fHCES==alphaHCES*fRCHC)       # p losses by harvested crops
@constraint(m, fHCF==xcorn*betapc*phuman)   # corn demand by urban (food center)

# Value balance part
@variable(m,cfSES)
@variable(m,cfMS)
@variable(m,cfHCF)
@variable(m,cfHCES)
@variable(m,cfRCES)
@variable(m,cfHCDC)
@variable(m,cfWSS)
@variable(m,cfHCB)
#@variable(m,cfSRC)
#@variable(m,cfRCS)
#@variable(m,cfRCHC)

@constraint(m, cfHCF==ccorn*fHCF/xcorn)
@constraint(m, cfHCDC==((ccorn*Xcorn)/xcorn+(csilage*Xsilage)/xsilage)*fHCDC)
@constraint(m, cfHCB==ccorn*fHCB/xcorn)
@constraint(m, cfMS==1/xmanuret*cmanuret*fMS)
@constraint(m, cfWSS==1/xsludge*csluge*fWSS)
@constraint(m, cfSES==cwaste*fSES)
@constraint(m, cfRCES==cwaste*fRCES)
@constraint(m, cfHCES==cwaste*fHCES)

#@constraint(m, cfRCHC==((ccorn*Xcorn/xcorn)+(csilage*Xsilage/xsilage))*fRCHC)
#@constraint(m, cfSRC==((ccorn*Xcorn/xcorn)+(csilage*Xsilage/xsilage))*fSRC)
#@constraint(m, cfRCS==cbiochar*fRCS)

@variable(m,BAS)
@constraint(m,BAS == cfHCF+cfHCDC+cfHCB-cfMS-cfWSS-cfCSS-cfSES-cfRCES-cfHCES)

#Food Distribution Sector(F)
# mass balance only for this sector
@variable(m,fFUS)     # P demand urban sector
@variable(m,fFES)     # P organic losses residues- food commodities
@variable(m,fDSF)     # P from demand of food commodities

@constraint(m, fFUS==fHCF+fDSF-fFES)    # mass balance in food commodities sector
@constraint(m, fFES==alphaFES*fHCF)     # P organic waste

#Urban Sector(US)
# mass balance part
@variable(m,fUSES)        # human sector losses
@variable(m,fUSWS)           # p human excreta disposal

@constraint(m, fUSWS==fFUS-fUSES)         # p mass balance of human sector
@constraint(m, fUSES==alphaUSES*fFUS)      # p losses in human sector

# value balance part
@variable(m,cfUSES)
@variable(m,cfUSWS)
@variable(m,cfFUS)
@variable(m,cfDSF)      # cost of the dairy products to the food commodoties sector

@constraint(m, cfFUS == cfHCF+cfDSF)
@constraint(m, cfUSES==cwaste*fUSES)
@constraint(m, cfUSWS==1/xdesh*cdesh*fUSWS)

#=I removed these two equations
@variable(m,cfFES)
@constraint(m, cfFUS==(0.1*cmilk+0.9*ccheese)*fFUS)
@constraint(m, cfFES==(cwaste*fFES))
=#

@variable(m,BUS)
@constraint(m,BUS == -cfFUS-cfUSWS-cfUSES)

#WWTP Sector(WS)
# mass balance part
@variable(m,fWSES)  # P WWTP losses

@constraint(m, fWSS==fUSWS-fWSES)         # p mass balance in wwtp sector
@constraint(m, fWSS==alphaWSS*fUSWS)       # organic fertilizer from the wwtp

# value balance part
@variable(m,cfWSES)
@constraint(m, cfWSES==cwaste*fWSES)

@variable(m,BWS)
@constraint(m,BWS == cfWSS + cfUSWS - cfWSES)

#Dairy Sector(DST)
# mass balance part
# dairy cattle sector (DC)
@variable(m,fDCES)    # dairy cattle losses
@variable(m,fDCM)        # p to manure
@variable(m,fDCDS)      # p from the dairy sector

@constraint(m, fDCDS==fCSDC+fHCDC-fDCES-fDCM)   # p mass balance dairy cattle sector
@constraint(m, fHCDC==xcorn*betarc*pcow+xsilage*betars*pcow)    # p to the dairy sector by corn and forage
@constraint(m, fDCES==alphaDCES*fHCDC)           # p from the dairy sector as dairy waste

#dairy production (DS)

@variable(m,fDSB)       # p dairy production losses
@variable(m,fpm)        # p in milk production
@variable(m,fpch)       # p in cheese production
@variable(m,pmilk)      # p milk production

@constraint(m, pmilk==betamilk*pcow)                # total milk production
@constraint(m, fpch==betamc*xcheese*betaprc*pmilk)   # p in cheese production by dairy sector
@constraint(m, fpm==(1-betamc)*xmilk*pmilk)        # p in milk production by dairy sector
@constraint(m, fDSF==fDCDS-fDSE-fDSB)                 # mass balance in dairy
@constraint(m, fDCDS==fpch+fpm)                     # p for each dairy product
@constraint(m, fDSF==(betamp*xmilk+betachp*xcheese)*phuman) # p dairy products to satisfy the urban demand

# value balance part
@variable(m,cfDCES)      # balance of money in dairy losses
@variable(m,cfDCM)       # money in manure production
@variable(m,cfDSE)      # money generates by dairy product exports
@variable(m,cfDSB)      # cost of the dairy produts for regulation

@constraint(m, cfDCM==1/xmanurer*cmanurer*fDCM)       #cost of p from the dairy sector to the manure stage
@constraint(m, cfDSB==(cmilk*Xmilk/xmilk+Xcheese*ccheese/xcheese)*fDSB)   #cost of the P flowrate
@constraint(m, cfDSE==(cmilk*Xmilk/xmilk+Xcheese*ccheese/xcheese)*fDSE)   #cost of the P in the dairy products from the dairy sector to regulation restriction
@constraint(m, cfDSF==(betamp*cmilk+betachp*ccheese)*phuman)             # cost of the dairy produts demanded by the urban sector
@constraint(m, cfDCES==cwaste*fDCES)     #cost of disposal daiy waste to the environment

#=I do not understand this, why divide by xsilage
@constraint(m, cfDCES==(cdairylos*fDCES)/xsilage)     #cost of disposal daiy waste to the environment
=#
#@constraint(m, cfDCDS==ccheese*betamc*betaprc*pmilk+cmilk*(1-betamc)*pmilk)  # cost of the total dairy products

@variable(m,BDST)
@constraint(m, BDST == cfDCM + cfDSB + cfDSE + cfDSF - cfHCDC - cfDCES - cfCSDC)

#Waste Processing Manure Sector(M)
# mass balance part
@variable(m,fMES)       # p losses from manure

@constraint(m, fMS==fDCM-fMES)                     # mass balance in manure production
@constraint(m, fMES==alphaMES*fDCM)                 # losses in manure production

# value balance part
@variable(m,cfMES)         # cost of the p losses from manure
@constraint(m, cfMES==cwaste*fMES)      # cost of the manure that is lost as discharges to the environment

@variable(m, BM)
@constraint(m, BM == cfMS - cfDCM - cfMES)

#Overall Value Balance
@variable(m, BTOTAL)
@constraint(m, BTOTAL == BMIS+BCS+BAS+BUS+BWS+BDST+BM)

#objetive function
#@objective(m, Max, BTOTAL)


#Solve the model and output result

status = solve(m)

#println("---------Model Gloabal Parameters----------")

#println("pcow=",getvalue(pcow))
#println("phuman=",getvalue(phuman))

println("----------Mass Balance Results----------")

println("----Mining Sector----")
println("fMIS=",fMIS)
println("fMISES=",getvalue(fMISES))
println("fMISCS=",getvalue(fMISCS))

println("----Chemical Sector----")
println("fMISCS=",getvalue(fMISCS))
println("fCSS=",getvalue(fCSS))
println("fCSES=",getvalue(fCSES))
println("fCSDC=",getvalue(fCSDC))

println("----Agriculture Sector----")
println("fSRC=",getvalue(fSRC))
println("fSES=",getvalue(fSES))
println("fCSS=",getvalue(fCSS))
println("fMS=",getvalue(fMS))
println("fWSS=",getvalue(fWSS))
println("fRCS=",getvalue(fRCS))
println("fRCES=",getvalue(fRCES))
println("fRCHC=",getvalue(fRCHC))
println("fHCF=",getvalue(fHCF))
println("fHCES=",getvalue(fHCES))
println("fHCDC=",getvalue(fHCDC))
println("fHCB=",getvalue(fHCB))

println("----Food Distribution----")
println("fHCF=",getvalue(fHCF))
println("fDSF=",getvalue(fDSF))
println("fFUS=",getvalue(fFUS))
println("fFES=",getvalue(fFES))

println("----Urban Sector----")
println("fFUS=",getvalue(fFUS))
println("fUSES=",getvalue(fUSES))
println("fUSWS=",getvalue(fUSWS))

println("----WWTP Sector----")
println("fUSWS=",getvalue(fUSWS))
println("fWSS=",getvalue(fWSS))
println("fWSES=",getvalue(fWSES))

println("----Dairy Sector----")
println("fCSDC=",getvalue(fCSDC))
println("fHCDC=",getvalue(fHCDC))
println("fDCDS=",getvalue(fDCDS))
println("fDCES=",getvalue(fDCES))
println("fDCM=",getvalue(fDCM))
println("fDSE=",fDSE)
println("fDSB=",getvalue(fDSB))
println("fDSF=",getvalue(fDSF))
println("fpm=",getvalue(fpm))
println("fpch=",getvalue(fpch))
println("pmilk=",getvalue(pmilk))

println("----Manure Management----")
println("fDCM=",getvalue(fDCM))
println("fMES=",getvalue(fMES))
println("fMS=",getvalue(fMS))

println("/n")
println("----------Value Balance Results----------")

println("----Mining Sector----")
println("cfMISES=",getvalue(cfMISES))
println("cfMISCS=",getvalue(cfMISCS))
println("bMIS=",getvalue(BMIS))

println("----Chemical Sector----")
println("cfCSDC=",getvalue(cfCSDC))
println("cfCSS=",getvalue(cfCSS))
println("cfCSES=",getvalue(cfCSES))
println("cfMISCS=",getvalue(cfMISCS))
println("bCS=",getvalue(BCS))

println("----Agriculture Sector----")
println("cfHCF=",getvalue(cfHCF))
println("cfHCDC=",getvalue(cfHCDC))
println("cfHCB=",getvalue(cfHCB))
println("cfMS=",getvalue(cfMS))
println("cfWSS=",getvalue(cfWSS))
println("cfCSS=",getvalue(cfCSS))
println("cfSES=",getvalue(cfSES))
println("cfRCES=",getvalue(cfRCES))
println("cfHCES=",getvalue(cfHCES))
println("bAS=",getvalue(BAS))

println("----Urban Sector----")
println("cfFUS=",getvalue(cfFUS))
println("cfUSWS=",getvalue(cfUSWS))
println("cfUSES=",getvalue(cfUSES))
println("bUS=",getvalue(BUS))

println("----WWTP Sector----")
println("cfWSS=",getvalue(cfWSS))
println("cfUSWS=",getvalue(cfUSWS))
println("cfWSES=",getvalue(cfWSES))
println("bWS=",getvalue(BWS))

println("----Dairy Sector----")
println("cfDCM=",getvalue(cfDCM))
println("cfDSB=",getvalue(cfDSB))
println("cfDSE=",getvalue(cfDSE))
println("cfDSF=",getvalue(cfDSF))
println("cfHCDC=",getvalue(cfHCDC))
println("cfDCES=",getvalue(cfDCES))
println("cfCSDC=",getvalue(cfCSDC))
println("bDST=",getvalue(BDST))

println("----Manure Management----")
println("cfMS=",getvalue(cfMS))
println("cfDCM=",getvalue(cfDCM))
println("cfMES=",getvalue(cfMES))
println("bM=",getvalue(BM))

println("----Overall Value----")
println("Btotal=",getvalue(BTOTAL))

#= Variables not used
println("cfDCDS=",getvalue(cfDCDS))
println("cfDSE=",getvalue(cfDSE))
println("cfDSF=",getvalue(cfDSF))
println("cfSRC=",getvalue(cfSRC))
println("cfRCS=",getvalue(cfRCS))
println("cfRCHC=",getvalue(cfRCHC))
println("cfFES=",getvalue(cfFES))
println("cfDCDS=",getvalue(cfDCDS))
println("bRC=",getvalue(BRC))
println("bHC=",getvalue(BHC))
println("bF=",getvalue(BF))
=#
