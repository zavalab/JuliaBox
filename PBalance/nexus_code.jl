#=
Balance of Phosphorus in the Dairy Nexus
Julia 0.6.4
Jan 2018: Created by Dulce Lopez Diaz
Apr 2019: Modified by Yicheng Hu
May 2019: Modified by Winnie Chan
=#

#= Sectors
ES  -environment
MIS -mining sector
CS  -chemical sector
S   -cultivated soil
RC  -crop production
HC  -crop processing (harvested crop)
F   -food distribution sector
US  -urban sector
WS  -wastewater treatment sector
M   -manure management sector
DC  -dairy livestock breeding
DS  -dairy product production
E   -exports
=#

# Define parameters for the mass balance
##Crop Parameters
xcorn=2.7036e-3   # P content in corn (tonne P/tonne corn)
xsilage=5.55e-4   # P content in silage (tonne P/tonne silage)
xhay=2.3608e-3    # P content in hay, alfalfa, haylage (tonne P/tonne hay, alfalfa, haylage)
xsoybean=6.05e-3  # P content in soybean (tonne P/tonne soybean)
xpotato=5.63e-4   # P content in potato (tonne P/tonne potato)
xwinwheat=3.31e-3 # P content in winter wheat (tonne P/tonne winter wheat)
xbeans=4.784e-3   # P content in beans (tonne P/tonne beans)
xsweetcorn=3.5e-3 # P content in sweetcorn (tonne P/tonne sweetcorn)
betapc=15.966e-3  # corn consumption rate (tonne/person/year)
betaps=0          # silage consumption rate (tonne/person/year)
betaph=0          # hay, alfalfa, haylage consumption rate (tonne/person/year)
betapsb=37.152e-3 # soybean consumption rate (tonne/person/year)
betapp=52.526e-3  # potato consumption rate (tonne/person/year)
betapw=59.783e-3  # winter wheat consumption rate (tonne/person/year)
betapb=3.4019e-3  # beans consumption rate (tonne/person/year)
betapsc=9.707e-3  # sweetcorn consumption rate (tonne/person/year)

##Mining Sector
alphaMIS=0.09  #P removal efficiency (tonne loss P/tonne extracted P from phosphate)

##Chemical Sector
alphaCS=0.11        # fraction of loss P in the chemical sector (tonne loss P/tonne P input to CS)
betaCSDC=2.5e-6     # requirements of additives by cow population (tonne P/cow/year)
#*** in P mass balance, the additives are not significant compared to other flows; but in money balances, the additives are important.

##Agriculture Sector
alphaSES=0.25      # fraction that represents the losses in the soil stage
alphaRCS=0.00      # fraction biochar recycle to cultivated soil (conversion factor in the biochar process)
alphaRCES=0.13     # fraction losses in crop activities, e.g. absorption or fixation
alphaHCES=0.09     # fraction losses from crops during processing

##Food Distribution Sector
alphaFES=0.12      # fraction for organic residues in food commodities (losses)

##Urban Sector
alphaUSES=0.25     # fraction losses in the human sector

##WWTP Sector
alphaWSS=0.75      # yield of residues used as fertilisers

##Dairy Sector
alphaDCES=0.18    # fraction of dairy cow losses ### true waste, not include incontrollable output to pasture
betamc=0.9        # fraction of milk that is transformed to cheese
xcheese=7.809e-3  # P content in cheese (tonne P/tonne cheese)
betaprc=1.03e-4   # cheese yield (tonne cheese/L milk)
xmilk=9.578e-7    # P content in milk (tonne P/L milk)
betamilk=10570    # annual production of milk per cow (L milk/cow/year)
betamp=67.682     # milk consumption rate (L milk/person/year)
betachp=0.0174    # cheese consumption rate (tonne cheese/person/year)
xbeef=6.10e-4     # P content in beef (tonne P/tonne beef)
betabfp=3.588e-2  # beef consumption rate (tonne beef/person/year)

##Manure Management Sector
#alphaMES=0.18     # factor of losses in the production of manure for transportation and disposal treatments



# Define parameters for economic analysis
##Mining Sector
crock=88        # price per tonne of phosphate rock (USD/tonne phosphate rock)
xrock=0.131     # P content in P rock (tonne P/tonne phosphate rock)
cwaste=0

##Chemical Sector
cadd=500        # price of additives (USD/tonne additive)
cfert=393.43    # price of fertiliser (USD/tonne fertiliser)
xadd=0.21       # P content in additives (tonne P/tonne additive)
xfert=0.20      # P content in DAP fertiliser (tonne P/tonne DAP)

##Agriculture Sector
ccorn=137.79        # price of corn (USD/tonne corn)
csilage=28          # price of silage (USD/tonne silgae)
chay=161.33         # price of hay, alfalfa, haylage (USD/tonne hay, alfalfa, haylage)
csoybean=316.18     # price of soybean (USD/tonne soybean)
cpotato=264.55      # price of potato (USD/tonne potato)
cwinwheat=165.35    # price of winter wheat (USD/tonne winter wheat)
cbeans=166.45       # price of beans (USD/tonne beans)
csweetcorn=92.37    # price of sweetcorn (USD/tonne sweetcorn)
cfeeds=96.75        # price of feeds and other feed grains (USD/tonne feeds)

gammahcdc=3.193e-3      # P per unit cattle feed (tonne P/tonne cattle feed)
gammahce=3.395e-3       # P per unit crop export (tonne P/tonne crop export)
kappahcecorn=0.295      # mass fraction of corn in crop exports
kappahcesoybean=0.234   # mass fraction of soybean in crop exports
kappahceothers=0.122    # mass fraction of other plant products in crop exports
kappahcefeeds=0.349     # mass fraction of feeds and other feed grains in crop exports

cmanuret=7.85     # price of treated manure (USD/tonne treated manure)
csludge=15        # price of sludge (USD/tonne sludge)
xmanuret=4.6e-3   # P content in manure (tonne P/tonne manure)
xsludge=9e-6      # P content in sludge (tonne P/tonne sludge)

##Urban Sector
cdesh=0.438         # price of wastewater gate fee (USD/tonne wastewater)
xdesh=1e-5          # P content in wastewater inflow (tonne P/tonne wastewater)

##Dairy Sector
xmanurer=9.2e-3   # P content in manure (tonne P/tonne manure)
cmanurer=7.56     # price of raw manure (USD/tonne raw manure)
cmilk=1.15        # price of milk (USD/L milk)
ccheese=724.4     # price of cheese (USD/tonne cheese)
cbeef=11786       # price of beef (USD/tonne beef)



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
#@variable(m, csilage)
#@variable(m, crock>=0)
#@variable(m, cfert>=0)

# Global variables
pcow=1280000    #cow population, heads
phuman=5819000  #human population
#fDSE=9455       #P exported in dairy products

## MASS BALANCE SECTION
# Mining Sector (MIS)
@variable(m,fMIS)       #total P extracted
@variable(m,fMISES)
@variable(m,fMISCS)
@constraint(m, fMISCS==fMIS-fMISES)
@constraint(m, fMISES==alphaMIS*fMIS)


# Chemical Sector (CS)
@variable(m,fCSES)  #P distribution losses
@variable(m,fCSS)   #P fertiliser from chemical sector
@variable(m,fCSDC)  #P additives to dairy cow feed
@constraint(m, fCSS==100100)           #set fixed number
@constraint(m, fCSS==fMISCS-fCSES-fCSDC) #mass balance in the chemical Sector
@constraint(m, fCSES==alphaCS*fMISCS)    #distribution and processing losses in the chemical sector
@constraint(m, fCSDC==betaCSDC*pcow)     #P from chemical sector to dairy sector as additives


# Agriculture Sector (AS)
# soil balance (S)
@variable(m,fSES)          #P erossion losses
@variable(m,fSRC)          #P crop content
@variable(m,fRCS)          #biochar
@variable(m,fWSS)          #P flow from sludge
@variable(m,fMS)           #P flow from manure

@constraint(m, fSRC==fCSS+fMS+fWSS-fSES+fRCS)      #mass balance in the cultivated soil sector
#@constraint(m, fSRC==fCSS+fMS+fWSS-fSES)            #mass balance in the cultivated soil sector
@constraint(m, fSES==alphaSES*(fCSS+fMS+fWSS+fRCS))                 #erosion losses
@constraint(m, fRCS==alphaRCS*fSRC)                #biochar recycle

# crop sector (RC)
@variable(m,fRCES)        #P crop losses
@variable(m,fRCHC)        #P crop production to sites of distribution-harvested crop sector
#@constraint(m, fRCHC==fSRC-fRCS-fRCES)     #mass balance in the crop production sector
@constraint(m, fRCHC==fSRC-fRCES)           #mass balance in the crop production sector
@constraint(m, fRCHC==79968.86)             # from excel "harvest sheet"
@constraint(m, fRCES==alphaRCES*fSRC)       #crop losses

# harvested crops (HC)
@variable(m,fHCF)        #P flow to food commodities
@variable(m,fHCES)       #P harvested losses
@variable(m,fHCDC)       #P crop products to dairy sector
@variable(m,fHCE)        #a balance P stream indicating export/import

cropDd = xcorn*betapc+xsilage*betaps+xhay*betaph+xsoybean*betapsb+xpotato*betapp+xwinwheat*betapw+xbeans*betapb+xsweetcorn*betapsc #P content of per capita food demand

@constraint(m,fHCF==fRCHC-fHCDC-fHCES-fHCE)     #balance in the harvested crop sector
@constraint(m,fHCES==alphaHCES*fRCHC)           #P losses by harvested crops
@constraint(m,fHCF==cropDd*phuman)              #crop demand by urban (food center)


# Food Distribution Sector(F)
# mass balance only for this sector
@variable(m,fFUS)     #P demand urban sector
@variable(m,fFES)     #P organic losses residues - food commodities
@variable(m,fDSF)     #P from demand of food commodities

@constraint(m, fFUS==fHCF+fDSF-fFES)    #mass balance in food commodities sector
@constraint(m, fFES==alphaFES*fHCF)     #P organic waste


# Urban Sector(US)
@variable(m,fUSES)        #human sector losses
@variable(m,fUSWS)        #P human excreta disposal

@constraint(m, fUSWS==fFUS-fUSES)         #P mass balance of human sector
@constraint(m, fUSES==alphaUSES*fFUS)     #P losses in human sector


# WWTP Sector (WS)
@variable(m,fWSES)  # P WWTP losses

@constraint(m, fWSS==fUSWS-fWSES)         #P mass balance in WWTP sector
@constraint(m, fWSS==alphaWSS*fUSWS)      #organic fertiliser from the WWTP


# Dairy Sector (DST)
# dairy cattle sector (DC)
@variable(m,fDCES)       #dairy cattle losses
@variable(m,fDCM)        #P to manure
@variable(m,fDCDS)       #P from the dairy cattle sector

@constraint(m, fDCDS==fCSDC+fHCDC-fDCES-fDCM)    #P mass balance dairy cattle sector
@constraint(m, fDCM==31335.44)                   #fix manure output from different cows
#@constraint(m, fHCDC==xcorn*betarc*pcow+xsilage*betars*pcow)    #P to the dairy sector by corn and forage
@constraint(m, fDCES==alphaDCES*fHCDC)           #P from the dairy sector as dairy waste

# dairy production (DS)
@variable(m,fDSE)

pmilk = betamilk*pcow
fpch = betamc*xcheese*betaprc*pmilk
fpm = (1-betamc)*xmilk*pmilk

bfweight=835744             # beef cattle slaughter liveweight (tonne/year)
#bfweight can also be calculated as (population of beef cattle)*(mass of each beef cow)*xbeef
fpbf = bfweight*xbeef

@constraint(m, fDSF==fDCDS-fDSE)                            #mass balance in dairy
@constraint(m, fDCDS==fpch+fpm+fpbf)                        #P for each dairy product
@constraint(m, fDSF==(betamp*xmilk+betachp*xcheese+betabfp*xbeef)*phuman) #P dairy products to satisfy urban demand


# Waste Processing Manure Sector (M)
@variable(m,fMES)       #P losses from manure
@constraint(m, fMS==fDCM-fMES)                       #mass balance in manure production
#@constraint(m, fMES==alphaMES*fDCM)                 #losses in manure production



##  VALUE BALANCE SECTION
# Mining Sector (MIS)
@variable(m,cfMISES)            #cost of mining treatment
@variable(m,cfMISCS)            #market price of phosphate rock
@variable(m,BMIS)
@constraint(m, BMIS == cfMISCS - cfMISES)
@constraint(m, cfMISCS == crock*fMISCS/xrock)
@constraint(m, cfMISES == cwaste*fMISES)


# Chemical Sector (CS)
@variable(m,cfCSES)     #cost to remove P in chemical waste
@variable(m,cfCSDC)     #feed additive cost
@variable(m,cfCSS)      #fertiliser cost
@variable(m,BCS)
@constraint(m, BCS == cfCSDC + cfCSS - cfMISCS - cfCSES)
@constraint(m, cfCSDC == (1/xadd)*cadd*fCSDC)
@constraint(m, cfCSS == (1/xfert)*cfert*fCSS)
@constraint(m, cfCSES == cwaste*fCSES)


# Agriculture Sector (AS)
@variable(m,cfMS)
@variable(m,cfWSS)
@variable(m,cfSES)
@variable(m,cfRCES)
@variable(m,cfHCDC)
@variable(m,cfHCF)
@variable(m,cfHCE)
@variable(m,cfHCES)
@variable(m,BAS)
@constraint(m,BAS == cfHCDC+cfHCF+cfHCE-cfMS-cfWSS-cfCSS-cfSES-cfRCES-cfHCES)
@constraint(m, cfHCDC == fHCDC/gammahcdc*cfeeds)
@constraint(m, cfHCF == phuman*(ccorn*betapc+csilage*betaps+chay*betaph+csoybean*betapsb+cpotato*betapp+cwinwheat*betapw+cbeans*betapb+csweetcorn*betapsc))
exportMass = fHCE/gammahce
@constraint(m, cfHCE == exportMass*(kappahcecorn*ccorn+kappahcesoybean*csoybean+kappahceothers*cpotato+kappahcefeeds*cfeeds))
@constraint(m, cfMS == 1/xmanuret*cmanuret*fMS)
@constraint(m, cfWSS == 1/xsludge*csludge*fWSS)
@constraint(m, cfSES == cwaste*fSES)
@constraint(m, cfRCES == cwaste*fRCES)
@constraint(m, cfHCES == cwaste*fHCES)


# Urban Sector (US)
@variable(m,cfUSES)
@variable(m,cfUSWS)
@variable(m,cfFUS)
@variable(m,cfDSF)      #cost of dairy products as food commodities
@variable(m,BUS)
@constraint(m, BUS == -cfFUS-cfUSWS-cfUSES)
@constraint(m, cfFUS == cfHCF+cfDSF)
@constraint(m, cfUSWS == 1/xdesh*cdesh*fUSWS)
@constraint(m, cfUSES == cwaste*fUSES)


# WWTP Sector (WS)
@variable(m,cfWSES)
@variable(m,BWS)
@constraint(m, BWS == cfWSS + cfUSWS - cfWSES)
@constraint(m, cfWSES == cwaste*fWSES)


# Dairy Sector (DST)
@variable(m,cfDCES)         #balance of money in dairy losses
@variable(m,cfDCM)          #revenue from manure production
@variable(m,cfDSE)          #revenue from dairy product exports
@variable(m,BDST)

fpdp = fpch+fpm+fpbf        #P flow in all dairy products
Kcheese = fpch/fpdp         #fraction of P in dairy products present in cheese
Kmilk = fpm/fpdp            #fraction of P in dairy products present in milk
Kbeef = fpbf/fpdp           #fraction of P in dairy products present in beef

@constraint(m, BDST == cfDCM + cfDSE + cfDSF - cfHCDC - cfDCES - cfCSDC)
@constraint(m, cfDCM == 1/xmanurer*cmanurer*fDCM)       #revenue/cost of manure from the dairy sector to the manure management sector
@constraint(m, cfDSE == (Kmilk*cmilk/xmilk+Kcheese*ccheese/xcheese+Kbeef*cbeef/xbeef)*fDSE)   #value of dairy sector exports
@constraint(m, cfDSF == (Kmilk*cmilk/xmilk+Kcheese*ccheese/xcheese+Kbeef*cbeef/xbeef)*fDSF)   #value of dairy sector produce for domestic market
@constraint(m, cfDCES == cwaste*fDCES)     #cost of dairy waste disposal to the environment


# Waste Processing Manure Sector (M)
@variable(m,cfMES)         #cost of P losses from manure
@variable(m, BM)
@constraint(m, BM == cfMS - cfDCM - cfMES)
@constraint(m, cfMES == cwaste*fMES)


# Overall Value Balance
@variable(m, BTOTAL)
@constraint(m, BTOTAL == BMIS+BCS+BAS+BUS+BWS+BDST+BM)



status = solve(m)


println("----------Mass Balance Results----------")

println("----Mining Sector----")
println("fMIS=",getvalue(fMIS))
println("fMISCS=",getvalue(fMISCS))
println("fMISES=",getvalue(fMISES))

println("----Chemical Sector----")
println("fMISCS=",getvalue(fMISCS))
println("fCSS=",getvalue(fCSS))
println("fCSDC=",getvalue(fCSDC))
println("fCSES=",getvalue(fCSES))

println("----Agriculture Sector----")
println("fSRC=",getvalue(fSRC))
println("fSES=",getvalue(fSES))
println("fCSS=",getvalue(fCSS))
println("fMS=",getvalue(fMS))
println("fWSS=",getvalue(fWSS))
#println("fRCS=",getvalue(fRCS))
println("fRCES=",getvalue(fRCES))
println("fRCHC=",getvalue(fRCHC))
println("fHCES=",getvalue(fHCES))
println("fHCDC=",getvalue(fHCDC))
println("fHCF=",getvalue(fHCF))
println("fHCE=",getvalue(fHCE))

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
println("fDSF=",getvalue(fDSF))
println("fDSE=",getvalue(fDSE))
#println("fpm=",getvalue(fpm))
#println("fpch=",getvalue(fpch))
#println("pmilk=",getvalue(pmilk))

println("----Manure Management----")
println("fDCM=",getvalue(fDCM))
println("fMS=",getvalue(fMS))
println("fMES=",getvalue(fMES))

println("/n")
println("----------Value Balance Results----------")

println("----Mining Sector----")
println("cfMISCS=",getvalue(cfMISCS))
println("cfMISES=",getvalue(cfMISES))
println("bMIS=",getvalue(BMIS))

println("----Chemical Sector----")
println("cfCSDC=",getvalue(cfCSDC))
println("cfCSS=",getvalue(cfCSS))
println("cfCSES=",getvalue(cfCSES))
println("cfMISCS=",getvalue(cfMISCS))
println("bCS=",getvalue(BCS))

println("----Agriculture Sector----")
println("cfHCDC=",getvalue(cfHCDC))
println("cfHCF=",getvalue(cfHCF))
println("cfHCE=",getvalue(cfHCE))
println("cfCSS=",getvalue(cfCSS))
println("cfMS=",getvalue(cfMS))
println("cfWSS=",getvalue(cfWSS))
println("cfSES=",getvalue(cfSES))
println("cfRCES=",getvalue(cfRCES))
println("cfHCES=",getvalue(cfHCES))
println("bAS=",getvalue(BAS))

println("----Urban Sector----")
println("cfFUS=",getvalue(cfFUS))
println("cfUSES=",getvalue(cfUSES))
println("cfUSWS=",getvalue(cfUSWS))
println("bUS=",getvalue(BUS))

println("----WWTP Sector----")
println("cfUSWS=",getvalue(cfUSWS))
println("cfWSS=",getvalue(cfWSS))
println("cfWSES=",getvalue(cfWSES))
println("bWS=",getvalue(BWS))

println("----Dairy Sector----")
println("cfCSDC=",getvalue(cfCSDC))
println("cfHCDC=",getvalue(cfHCDC))
println("cfDCES=",getvalue(cfDCES))
println("cfDCM=",getvalue(cfDCM))
println("cfDSF=",getvalue(cfDSF))
println("cfDSE=",getvalue(cfDSE))
println("bDST=",getvalue(BDST))

println("----Manure Management----")
println("cfDCM=",getvalue(cfDCM))
println("cfMS=",getvalue(cfMS))
println("cfMES=",getvalue(cfMES))
println("bM=",getvalue(BM))

println("----Overall Value----")
println("Btotal=",getvalue(BTOTAL))
