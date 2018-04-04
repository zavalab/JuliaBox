#C8 production model (whole, 444 nodes)
#consider mixing technology: yield factor linear combination
#C8 yiled factor based on 2.5% COD conversion, C6 7.5%

using JuMP
using Gurobi

m = Model(solver=GurobiSolver(Threads = 1,MIPGap = 5e-2, NodefileStart=0.25, TimeLimit = 86400));
#MIPGap = 1e-2; #Default value is 1e-6. #Not larger than 3e-2

#Importing Data
#technology_matrix = readdlm("technology_matrix.csv",','); # all 3 kinds of tech, each with several capacities available
node_matrix = readdlm("node_matrix.csv",',');              # all 5 kinds of nodes:county,CAFO,WWTP,LF,collection site
#product_matrix = readdlm("product_matrix.csv",',');        # all 6 kinds of products
demand_matrix = readdlm("demand_matrix.csv",',');
#alpha_matrix = readdlm("alpha_matrix.csv",',');      
supply_matrix = readdlm("supply_matrix.csv",',');

##Data pretreatment

#Yield_Factor_Part
C6_COD = 0.163
C8_COD = 0.017
CH4_vol= 0.51
VS_w     = [0.1115;0.0520;0.0800;0.10835;0.0548];
TS_w     = [0.1316;0.0630;0.1000;0.12844;0.0667];
CH4_yd_w = [0.01491038;0.01761488;0.02805254;0.1*0.02805254 + 0.9*0.01491038;0.1*0.02805254 + 0.9*0.01761488];
FD_mois  = 0.3

CH4_m  = (CH4_vol*16)/(CH4_vol*16 + (1-CH4_vol)*44)
bg_yd_w  = CH4_yd_w/CH4_m
FD_yd_w  = TS_w-bg_yd_w+FD_mois

COD_w  = bg_yd_w/0.7712*1000
C6_yd  = COD_w*C6_COD/8/32*120/1000
C8_yd  = COD_w*C8_COD/11/32*144/1000

COD_ID   = COD_w*(1-C6_COD-C8_COD)
bg_yd_ID = COD_ID*0.7712/1000
VS_ID    = VS_w*(1-C6_COD-C8_COD)
TS_ID    = VS_ID + (TS_w-VS_w)
FD_yd_ID = TS_ID - bg_yd_ID + FD_mois


alpha_matrix = [
    -1 0 0 C8_yd[1] C6_yd[1] 1-C8_yd[1]-C6_yd[1] 0 0 0 0 0 0 0 0 0 0;
    0 -1 0 C8_yd[2] C6_yd[2] 0 1-C8_yd[2]-C6_yd[2] 0 0 0 0 0 0 0 0 0;
    0 0 -1 C8_yd[3] C6_yd[3] 0 0 1-C8_yd[3]-C6_yd[3] 0 0 0 0 0 0 0 0;
    -0.9 0 -0.1 C8_yd[4] C6_yd[4] 0 0 0 1-C8_yd[4]-C6_yd[4] 0 0 0 0 0 0 0;
    0 -0.9 -0.1 C8_yd[5] C6_yd[5] 0 0 0 0 1-C8_yd[5]-C6_yd[5] 0 0 0 0 0 0;
    -1 0 0 0 0 0 0 0 0 0 bg_yd_w[1] FD_yd_w[1] 0 0 0 0;
    0 -1 0 0 0 0 0 0 0 0 bg_yd_w[2] 0 FD_yd_w[2] 0 0 0;
    0 0 -1 0 0 0 0 0 0 0 bg_yd_w[3] 0 0 FD_yd_w[3] 0 0;
    -0.9 0 -0.1 0 0 0 0 0 0 0 bg_yd_w[4] 0 0 0 FD_yd_w[4] 0;
    0 -0.9 -0.1 0 0 0 0 0 0 0 bg_yd_w[5] 0 0 0 0 FD_yd_w[5];
    0 0 0 0 0 -1 0 0 0 0 bg_yd_ID[1] FD_yd_ID[1] 0 0 0 0;
    0 0 0 0 0 0 -1 0 0 0 bg_yd_ID[2] 0 FD_yd_ID[2] 0 0 0;
    0 0 0 0 0 0 0 -1 0 0 bg_yd_ID[3] 0 0 FD_yd_ID[3] 0 0;
    0 0 0 0 0 0 0 0 -1 0 bg_yd_ID[4] 0 0 0 FD_yd_ID[4] 0;
    0 0 0 0 0 0 0 0 0 -1 bg_yd_ID[5] 0 0 0 0 FD_yd_ID[5]
    ];


#Tech_Cost_Part
#waste_biogas part
caps_bg  = [30000; 50000; 150000; 250000; 15000; 35000; 70000; 100000; 5000; 10000; 20000; 30000; 30000; 50000; 150000; 250000; 15000; 35000; 70000; 100000];
AD_bg = 937.12*caps_bg.^0.6 + 75355;
SP_bg = 17869*log(caps_bg) - 95066;
EG_bg = AD_bg*0.67851070;
OM_bg = AD_bg*0.09649075;
cln_bg = [caps_bg[1:4]*bg_yd_w[1]*0.08/1.15*1000;caps_bg[5:8]*bg_yd_w[2]*0.08/1.15*1000;caps_bg[9:12]*bg_yd_w[3]*0.08/1.15*1000;caps_bg[13:16]*bg_yd_w[4]*0.08/1.15*1000;caps_bg[17:20]*bg_yd_w[5]*0.08/1.15*1000];
inv_bg = 1.231*(AD_bg + EG_bg) + SP_bg;
opr_bg = 1.231*(OM_bg + cln_bg)+ 0.048*SP_bg;


#ID_biogas part 
caps_ID  = caps_bg
AD_ID = 937.12*caps_ID.^0.6 + 75355;
SP_ID = 17869*log(caps_ID) - 95066;
EG_ID = AD_ID*0.67851070;
OM_ID = AD_ID*0.09649075;
cln_ID = [caps_ID[1:4]*bg_yd_ID[1]*0.08/1.15*1000;caps_ID[5:8]*bg_yd_ID[2]*0.08/1.15*1000;caps_ID[9:12]*bg_yd_ID[3]*0.08/1.15*1000;caps_ID[13:16]*bg_yd_ID[4]*0.08/1.15*1000;caps_ID[17:20]*bg_yd_ID[5]*0.08/1.15*1000];
inv_ID = 1.231*(AD_ID + EG_ID)+SP_ID;
opr_ID = 1.231*(OM_ID + cln_ID)+0.048*SP_ID;

#C6C8 part
caps_CC  = caps_bg;
C6C8_CC  = [caps_CC[1:4]*(C8_yd[1]+C6_yd[1]); caps_CC[5:8]*(C8_yd[2]+C6_yd[2]); caps_CC[9:12]*(C8_yd[3]+C6_yd[3]); caps_CC[13:16]*(C8_yd[4]+C6_yd[4]); caps_CC[17:20]*(C8_yd[5]+C6_yd[5])];
AD_CC    = (caps_CC/3629793.6).^0.6*7538395
Pump1_CC = (caps_CC/3629793.6).^0.6*224943
Pump2_CC = (caps_CC/3629793.6).^0.6*29540
Sepex_CC = (C6C8_CC/18044.8992).^0.6*859554.23
Pump3_CC = (C6C8_CC/18044.8992).^0.6*29230
SepC1_CC = (C6C8_CC/18044.8992).^0.6*1344717
SepC2_CC = (C6C8_CC/18044.8992).^0.6*796836
Inv_CC   = 1.202*(AD_CC + Pump1_CC + Pump2_CC + Sepex_CC + Pump3_CC + SepC1_CC + SepC2_CC);
OM_CC    = 0.1*Inv_CC;
heat_CC0 = 1.1*caps_CC*1000*2.8*(35+273-280.5)/1000/55.6/0.7*35.314667*10/1000/0.8
heat_CC  = (C6C8_CC/18044.8992)*(3.041940863+0.323383717)*24*365*4.184/55.6/0.7*35.314667*10/0.8
base_CC  = (caps_CC/3629793.6)*16840000
opr_CC   = OM_CC + heat_CC + heat_CC0+ 1.202*base_CC;

#Tech matrix
technology_matrix = Matrix(60, 6);
technology_matrix[:,1] = [
    "tA1";"tA2";"tA3";"tA4";"tB1";"tB2";"tB3";"tB4";"tC1";"tC2";"tC3";"tC4";"tD1";"tD2";"tD3";"tD4";"tE1";"tE2";"tE3";"tE4";
    "tF1";"tF2";"tF3";"tF4";"tG1";"tG2";"tG3";"tG4";"tH1";"tH2";"tH3";"tH4";"tI1";"tI2";"tI3";"tI4";"tJ1";"tJ2";"tJ3";"tJ4";
    "tK1";"tK2";"tK3";"tK4";"tL1";"tL2";"tL3";"tL4";"tM1";"tM2";"tM3";"tM4";"tN1";"tN2";"tN3";"tN4";"tO1";"tO2";"tO3";"tO4"]
technology_matrix[:,2] = [caps_CC;caps_bg;caps_ID];
technology_matrix[:,3] = [Inv_CC;inv_bg;inv_ID];
technology_matrix[:,4] = [opr_CC;opr_bg;opr_ID];
technology_matrix[:,5] = [
    "p1";"p1";"p1";"p1";"p2";"p2";"p2";"p2";"p3";"p3";"p3";"p3";"p1";"p1";"p1";"p1";"p2";"p2";"p2";"p2";
    "p1";"p1";"p1";"p1";"p2";"p2";"p2";"p2";"p3";"p3";"p3";"p3";"p1";"p1";"p1";"p1";"p2";"p2";"p2";"p2";
    "p6";"p6";"p6";"p6";"p7";"p7";"p7";"p7";"p8";"p8";"p8";"p8";"p9";"p9";"p9";"p9";"p10";"p10";"p10";"p10"];
technology_matrix[:,6] = [(heat_CC+heat_CC0)/183.4528156;zeros(20,1);zeros(20,1)];


#Carbon_Emission_Part
function emission(VS_frac,Bo,frac,t)
    temp = readdlm("temperatures.csv",',');
    temp = 273.15 + (temp-32)/1.8;
    
    if t == 1
        VS_s = VS_frac*0.37;
        VS_l = VS_frac*(1-0.37)
        
        VS_tot_l = zeros(365);
        VS_loss_l = zeros(365);
        CH4_emit_l = zeros(365);
        VS_tot_s = zeros(365);
        VS_loss_s = zeros(365);
        CH4_emit_s = zeros(365);
        
        for i in [1,184]
            VS_tot_s[i] = VS_s
            VS_tot_l[i] = VS_l
            CH4_emit_l[i] = (24*VS_tot_l[i]*frac/1000*exp(43.33-112700/8.314/temp[i])) + 0.01*(24*VS_tot_l[1]*(1-frac)/1000*exp(43.33-112700/8.314/temp[1]));
            CH4_emit_s[i] = frac * max(VS_tot_s[i]*Bo*0.67*(0.201*(temp[i]-273.15)-0.29)/100/100,0);
            VS_loss_l[i] = 3*CH4_emit_l[i];
            VS_loss_s[i] = 3*CH4_emit_s[i];
        end
        
        for i in 2:183
            VS_tot_s[i] = VS_s * i - VS_loss_s[i-1];
            VS_tot_l[i] = VS_l * i - VS_loss_l[i-1];
            
            CH4_emit_l[i] = (24*VS_tot_l[i]*frac/1000*exp(43.33-112700/8.314/temp[i])) + 0.01*(24*VS_tot_l[i]*(1-frac)/1000*exp(43.33-112700/8.314/temp[i]));
            CH4_emit_s[i] = frac * max(VS_tot_s[i]*Bo*0.67*(0.201*(temp[i]-273)-0.29)/100/100,0);

            VS_loss_l[i] = VS_loss_l[i-1] + 3*CH4_emit_l[i];
            VS_loss_s[i] = VS_loss_s[i-1] + 3*CH4_emit_s[i];
        end
        
        for i in 185:365
            VS_tot_s[i] = VS_s * (i-183) - VS_loss_s[i-1];
            VS_tot_l[i] = VS_l * (i-183) - VS_loss_l[i-1];
            
            CH4_emit_l[i] = (24*VS_tot_l[i]*frac/1000*exp(43.33-112700/8.314/temp[i])) + 0.01*(24*VS_tot_l[i]*(1-frac)/1000*exp(43.33-112700/8.314/temp[i]));
            CH4_emit_s[i] = frac * max(VS_tot_s[i]*Bo*0.67*(0.201*(temp[i]-273)-0.29)/100/100,0);

            VS_loss_l[i] = VS_loss_l[i-1] + 3*CH4_emit_l[i];
            VS_loss_s[i] = VS_loss_s[i-1] + 3*CH4_emit_s[i];
        end
        
        factor_l = 0.5*sum(CH4_emit_l[1:183])/183 + sum(CH4_emit_l[184:365])/182;
        factor_s = 0.5*sum(CH4_emit_s[1:183])/183 + sum(CH4_emit_s[184:365])/182;
        if VS_s <= 0.07
            factor_s = 1.4*factor_s;
        end
        if VS_l <= 0.07
            factor_l = 1.4*factor_l;
        end
        return factor_s+factor_l
    end

    
    if t == 0
        VS_tot = zeros(365);
        VS_loss = zeros(365);
        CH4_emit = zeros(365);
        
        for i in [1,184]
            VS_tot[i] = VS_frac
            CH4_emit[i] = (24*VS_tot[i]*frac/1000*exp(43.33-112700/8.314/temp[i])) + 0.01*(24*VS_tot[1]*(1-frac)/1000*exp(43.33-112700/8.314/temp[1]));
            VS_loss[i] = 3*CH4_emit[i];
        end
        
        for i in 2:183
            VS_tot[i] = VS_frac * i - VS_loss[i-1];
            CH4_emit[i] = (24*VS_tot[i]*frac/1000*exp(43.33-112700/8.314/temp[i])) + 0.01*(24*VS_tot[i]*(1-frac)/1000*exp(43.33-112700/8.314/temp[i]));
            VS_loss[i] = VS_loss[i-1] + 3*CH4_emit[i];
        end
        
        for i in 185:365
            VS_tot[i] = VS_frac * (i-183) - VS_loss[i-1];
            CH4_emit[i] = (24*VS_tot[i]*frac/1000*exp(43.33-112700/8.314/temp[i])) + 0.01*(24*VS_tot[i]*(1-frac)/1000*exp(43.33-112700/8.314/temp[i])); 
            VS_loss[i] = VS_loss[i-1] + 3*CH4_emit[i];
        end
        factor = 0.5*sum(CH4_emit[1:183])/183 + sum(CH4_emit[184:365])/182;
        if VS_frac <= 0.07
            factor = 1.4*factor;
        end
        return factor
    end
end

VS_matrix = [0.1115;0.0520;0.0800;0;0;VS_ID;0;VS_w*0.52]
product_alias = [
    "p1";"p2";"p3";"p4";"p5";"p6";"p7";"p8";"p9";"p10";"p11";"p12";"p13";"p14";"p15";"p16"]    
product_names  = [
    "Manure";"Sludge";"Food Waste";"C8";"C6";"ID1";"ID2";"ID3";"ID4";"ID5";"Biogas";"FD1";"FD2";"FD3";"FD4";"FD5"]
trans_cost    = 0.16*ones(16,1);
product_matrix = Matrix(16,5);
product_matrix[:,1] = product_alias;
product_matrix[:,2] = product_names;
product_matrix[:,3] = trans_cost;
product_matrix[:,4] = [
    emission(VS_matrix[1],0.24,0.5,0);emission(VS_matrix[2],0.24,0.5,0);emission(VS_matrix[3],0.6,0.5,0);0;0;
    emission(VS_matrix[6],0.24,0.5-C6_COD-C8_COD,0);emission(VS_matrix[7],0.24,0.5-C6_COD-C8_COD,0);
    emission(VS_matrix[8],0.6,0.5-C6_COD-C8_COD,0);emission(VS_matrix[9],0.276,0.5-C6_COD-C8_COD,0);
    emission(VS_matrix[10],0.276,0.5-C6_COD-C8_COD,0);
    0;emission(VS_matrix[12],0.24,5/55,1);emission(VS_matrix[13],0.24,5/55,1);emission(VS_matrix[14],0.6,5/55,1);
    emission(VS_matrix[15],0.276,5/55,1);emission(VS_matrix[16],0.276,5/55,1)]
product_matrix[:,5] = product_matrix[:,4]/16*5*44;
product_matrix[11,5] = 44/(CH4_vol*16+(1-CH4_vol)*44);
    
#Saved Social Cost of Carbon Constant
η = 0.3 # efficiency to generate electricity
SCC_bg = 55.6*CH4_m*1000*1000/3600*η*0.84/1000


#Define Sets

#Numbers about intervals

n1=72;
n2=317;
n3=379;

TECHS = technology_matrix[:,1];                   # set of all technologies
TECH_C8ALL = technology_matrix[1:20,1];           # set of all C8 technologies
TECH_BgALL = technology_matrix[21:end,1];         # set of all Biogas technologies
TECH = Matrix(15,1);
for i in 1:15
    TECH[i] = technology_matrix[Int((i-1)*4+1):Int(4*i),1];
end


NODES = node_matrix[:,1];                         # set of all nodes
NODE1 = node_matrix[1:n1,1];                      # set of nodes1-county
NODE2 = node_matrix[n1+1:n2,1];                   # set of nodes2-CAFO
NODE3 = node_matrix[n2+1:n3,1];                   # set of nodes3-WWTP
NODE4 = node_matrix[n3+1:end-1,1];                # set of nodes4-LF
NODE5 = node_matrix[end,1];                       # set of nodes5-collection site

PRODS = product_matrix[:,1];                      # set of products
DEMS  = demand_matrix[:,1];                       # set of demands
SUPS  = supply_matrix[:,1];                       # set of supplies



#Define Dictionaries
tech_cap    =   Dict(zip(TECHS, technology_matrix[:, 2]));                     # technology capacity tonne/year
tech_alias  =   Dict(zip(TECHS, technology_matrix[:, 1]));                     # technology name alias
tech_invcost    =   Dict(zip(TECHS, technology_matrix[:, 3]));                 # technology investment cost $/year
tech_opcost =   Dict(zip(TECHS, technology_matrix[:, 4]));                     # technology operating cost $/year
tech_refprod    =   Dict(zip(TECHS, technology_matrix[:, 5]));                 # technology products 
tech_co2    =   Dict(zip(TECHS, technology_matrix[:, 6]));                     # technology CO2 emission tonne CO2/year

node_lat    =   Dict(zip(NODES, node_matrix[:, 3]));                           # node latitude
node_long   =   Dict(zip(NODES, node_matrix[:, 4]));                           # node longitude
node_alias  =   Dict(zip(NODES, node_matrix[:, 2]));                           # node alias name
prod_alias  =   Dict(zip(PRODS, product_matrix[:, 2]));                        # product alias name

prod_transcost  =   Dict(zip(PRODS, product_matrix[:, 3]));                    # product transportation cost
prod_ch4 = Dict(zip(PRODS, product_matrix[:, 4]));                           # product CH4 emission coefficient
prod_co2 = Dict(zip(PRODS, product_matrix[:, 5]));                           # product CO2 emission coefficient (not contribute to GW)

dem_node    =   Dict(zip(DEMS, demand_matrix[:, 2]));                          # demand node
dem_prod    =   Dict(zip(DEMS, demand_matrix[:, 3]));                          # demand product
dem_cap     =   Dict(zip(DEMS, demand_matrix[:, 4]));                          # demand flow capacity
dem_price   =   Dict(zip(DEMS, demand_matrix[:, 5]));                          # demand price

sup_node    =   Dict(zip(SUPS, supply_matrix[:, 2]));                          # supply node
sup_prod    =   Dict(zip(SUPS, supply_matrix[:, 3]));                          # supply product
sup_value   =   Dict(zip(SUPS, supply_matrix[:, 4]));                          # supply value
sup_price   =   Dict(zip(SUPS, supply_matrix[:, 5]));                          # supply price



#Define each prod from and go sets
NOTFROM = Dict("p1"=>[NODE1;NODE3;NODE4;NODE5], "p2"=>[NODE2;NODE4;NODE5], "p3"=>[NODE2;NODE3;NODE4;NODE5] , "p4"=>[NODE1;NODE5],"p5"=>[NODE1;NODE5], "p6"=>[NODE1;NODE3;NODE4;NODE5], "p7"=>[NODE1;NODE2;NODE5], "p8"=>[NODE1;NODE2;NODE3;NODE5],"p9"=>[NODE1;NODE3;NODE4;NODE5],"p10"=>[NODE1;NODE2;NODE5],"p11"=>[NODE1;NODE5],"p12"=>[NODE1;NODE5],"p13"=>[NODE1;NODE5],"p14"=>[NODE1;NODE5],"p15"=>[NODE1;NODE5],"p16"=>[NODE1;NODE5]);
NOTTO   = Dict("p1"=>[NODE1;NODE3;NODE4;NODE5],"p2"=>[NODE1;NODE2;NODE5],"p3"=>[NODE1;NODE5],"p4"=>[NODE1;NODE2;NODE3;NODE4], "p5"=>[NODE1;NODE2;NODE3;NODE4], "p6"=>[NODE3;NODE5],"p7"=>[NODE2;NODE5],"p8"=>[NODE2;NODE3;NODE5], "p9"=>[NODE3;NODE5],"p10"=>[NODE2;NODE5],"p11"=>[NODE1;NODE2;NODE3;NODE4;NODE5], "p12"=>[NODE3;NODE5], "p13"=>[NODE3;NODE5], "p14"=>[NODE3;NODE5], "p15"=>[NODE3;NODE5], "p16"=>[NODE3;NODE5]);


# Emissions Metric
const co2_per_km  =   0.2e-3                         # 0.2e-3ton of CO2 emitted per km per ton of freight



# TradeOff Analysis Parameters
##cost_min    =   0
##cost_max    =   5000.000000000002        
#budget      =   15e5*365

epsilon     =   1                           # Varied between 0 and 1 to manipulate the budget

# Haversine Formula Parameters. Used to estimate distance from latitude and longitude data
R = 6335.439

## Defining Two Variable Dictionaries ##
transfer = Dict(("tA1","p1") => 0.5)         #Just used as an initiator to set up the dictionary with two keys
transfer = Dict(("tA1","p1") => 0.5)         #Just used as an initiator to set up the dictionary with two keys
for i in 1:15
    for j in 1:4
      for k in 1: length(PRODS)
        transfer[(TECH[i][j], PRODS[k])] = alpha_matrix[i, k]
      end
    end
end

distance = Dict(("n1", "n2") => 1.1)

# Using the Haversine formula
for i in NODES
  for j in NODES
    distance[(i, j)] = 2*R*asin(sqrt(sin((node_lat[j] - node_lat[i])*pi/2/180)^2 + cos(node_lat[j]*pi/180)*cos(node_lat[i]*pi/180)*sin((node_long[j] - node_long[i])*pi/2/180)^2))
  end
end

M=1e15; #big M
#Rmax=200; #Define largest transportation distance

## Define variables
#flows

@variable(m,f[NODES,NODES,PRODS]>=0);

for p in PRODS
@constraint(m, [i in NODES, j in NOTTO[p]], f[i,j,p] ==0);
@constraint(m, [i in NOTFROM[p], j in NODES],f[i,j,p] ==0);             
end
#for i in NODES
    #for j in NODES
        #if distance[i,j] >= Rmax
            #for pr in ["p1";"p2";"p3"]
                #@constraint(m, toofar[pr in ["p1";"p2";"p3"]], f[i,j,pr] == 0);
            #end
        #end
    #end
#end

#logic variable--if tech installed
@variable(m, y[TECHS,NODES],Bin);
@constraint(m, notech1[i in [NODE1;NODE5], j in TECHS], y[j,i] == 0);
@constraint(m, notech2[i in NODE2, j in [TECH[2];TECH[3];TECH[5];TECH[7];TECH[8];TECH[10];TECH[12];TECH[13];TECH[15]]], y[j,i] == 0);
@constraint(m, notech3[i in NODE3, j in [TECH[1];TECH[3];TECH[4];TECH[6];TECH[8];TECH[9];TECH[11];TECH[13];TECH[14]]], y[j,i] == 0);
@constraint(m, notech4[i in NODE4, j in [TECH[1];TECH[2];TECH[6];TECH[7];TECH[11];TECH[12]]], y[j,i] ==0);

#demand
@variable(m, dem[DEMS] >=0);
@variable(m, d[NODES,PRODS] >=0);
@variable(m, sup[SUPS] >=0);
@variable(m, s[NODES,PRODS] >=0);
#generated/consumed amount by tech t in node i
@variable(m, x[NODES,PRODS,TECHS]);
#total x (sum over t)

@variable(m, p[NODES, PRODS]);

@variable(m,transcost>=0);
@variable(m,invcost>=0);
@variable(m,opcost>=0);

#social welfare
@variable(m,swf);
#unit: Million $/year $/tonne * Million tonne/year = Million $/year
@constraint(m, swf == sum(sup_price[s]*sup_value[s] for s in SUPS) - sum(dem_price[d]*dem[d] for d in DEMS)); 

#demand and supply
#unit: Million tonne/year
@constraint(m, demeq[n in NODES, pr in PRODS], d[n,pr] == sum(dem[dd] for dd in DEMS if dem_prod[dd]==pr && dem_node[dd]==n));
@constraint(m, supeq[n in NODES, pr in PRODS], s[n,pr] == sum(sup[ss] for ss in SUPS if sup_prod[ss]==pr && sup_node[ss]==n));

## Balance Constraints
@constraint(m,balance[i in NODES, pr in PRODS],s[i,pr]+p[i,pr]+sum(f[j,i,pr] for j in NODES) == sum(f[i,j,pr]  for j in NODES)+d[i,pr]);
@constraint(m,process[i in NODES, pr in PRODS], p[i,pr] == sum(x[i,pr,t] for t in TECHS));    
@constraint(m, transfer_pr[i in NODES, t in TECHS, pr in PRODS], x[i,pr,t] == transfer[t,pr]/transfer[t,tech_refprod[t]]*x[i,tech_refprod[t],t]);

##Logic & capacity constriants
#tech capacity & x=0 if no tech installed
@constraint(m, techonofflb[i in NODES, t in TECHS, pr in PRODS], x[i,pr,t] >= -y[t,i]*tech_cap[t]*1e-6);
@constraint(m, techonoffub[i in NODES, t in TECHS, pr in PRODS], x[i,pr,t] <= +y[t,i]*tech_cap[t]*1e-6);

#at most one capacity in each tech can be installed
#@constraint(m, onetech[i in [NODE2;NODE3;NODE4],pr in ["p1";"p2";"p3"]], sum{y[t,i], t in TECH_PROD[pr]} <= 1);
#at most one tech can be installed
@constraint(m, notech11[i in NODES], sum(y[t,i] for t in TECH_BgALL) <=0);
@constraint(m, notech33[i in NODES], sum(y[t,i] for t in TECH_C8ALL) <=1);
#@constraint(m, notech22[i in NODES], sum(y[t,i] for t in [TECH2;TECH3;TECH4;TECH5;TECH7;TECH8;TECH10;TECH11;TECH12;TECH13] i in NODES) ==0 );
                                    
#eliminate self-flow
@constraint(m, selfflow[i in NODES,pr in PRODS], f[i,i,pr] == 0);

#demand capacity 
# unit: Million tonne/year                                
@constraint(m, demand_capacity[i in DEMS], dem[i] <= dem_cap[i]*1e-6);
#fix supply
@constraint(m, supply_fix[i in SUPS], sup[i] == sup_value[i]*1e-6);

## Budget constraints
#unit: $*binary -> *1e-6 = Million $/year
@constraint(m, invcost == 1e-6*sum(tech_invcost[t]*y[t,n] for t in TECHS for n in NODES));
#unit: Million $/Million ton/km*km*Million ton/year =Million $/year                                            
@constraint(m, transcost == sum(prod_transcost[p]*f[n,nn,p]*distance[n,nn] for n in NODES for nn in NODES for p in PRODS));
#unit: $/year -> Million $/year                                                                    
@constraint(m, opcost == 1e-6*sum(tech_opcost[t]*y[t, n] for t in TECHS for n in NODES));
                                                                                    
@variable(m, totalcost>=0)
@constraint(m, transcost + opcost + invcost/20 == totalcost);
#@constraint(m, totalcost<= epsilon*budget);

##Carbon Cost
@variable(m, totalprod[PRODS]>=0);
@variable(m, ch4prod[PRODS] >=0);                                                                                    
@variable(m, co2prod[PRODS] >=0);                                                                                    
@constraint(m, totald[p in PRODS], totalprod[p] == sum(d[i,p] for i in NODES)); # million tonne/year
@constraint(m, cal_ch4[p in PRODS], ch4prod[p] == totalprod[p]*prod_ch4[p]);  # million tonne CH4/year
#@constraint(m, cal_co2[p in PRODS], co2prod[p] == totalprod[p]*prod_co2[p]);                                                                                    

@variable(m, CO2trans>=0);
@constraint(m, CO2trans == co2_per_km*sum(f[n,nn,p]*distance[n,nn] for n in NODES for nn in NODES for p in PRODS));
                       #million ton CO2/km/million ton * km * million ton/year = million ton CO2/year                                                                                 

@variable(m, unamp>=0);                                                                                    
@variable(m, unsp >=0);                                                                                    
@variable(m, unfwp>=0);                                                                                    
@constraint(m, unamp==totalprod["p1"]/11217129.219*1e6);
@constraint(m, unsp ==totalprod["p2"]/2512392.157*1e6);
@constraint(m, unfwp==totalprod["p3"]/634847.07*1e6);
                                                                                                        
@variable(m, totalch4>=0);
@variable(m, totalco2>=0);                                                                                                        
@variable(m, SCC_ch4>=0);
@variable(m, SCC_co2>=0);                                                                                                        
@constraint(m, totalch4 == sum(ch4prod[p] for p in PRODS));
@constraint(m, totalco2 == CO2trans + 1e-6*sum(y[t,i]*tech_co2[t] for t in TECHS for i in NODES));                                                                                                        
@constraint(m, SCC_ch4 == 28*152*totalch4); # Million $/Million tonne CH4 * Million ton CH4
@constraint(m, SCC_co2 == 152*totalco2);                                                                                                        

@variable(m, TSCC>=0);
@variable(m, SCO2>=0)                                                                                                        
@variable(m, SSCC>=0);                                                                                                        
@variable(m, NSCC);
@constraint(m, TSCC == SCC_ch4+SCC_co2);
@constraint(m, SCO2 == SCC_bg*totalprod["p11"]); # million tonne CO2/year
@constraint(m, SSCC == 152*SCO2); # Million $/year
@constraint(m, NSCC == TSCC-SSCC)                                                                                                        

##Objective                                                                                                         
@variable(m,profit);
@constraint(m, profit == -swf-totalcost);
@variable(m, Bgtech>=0);                                                                                        
@variable(m, C8tech>=0);    
@constraint(m, Bgtech == sum(y[t,i] for t in TECH_BgALL for i in NODES));
@constraint(m, C8tech == sum(y[t,i] for t in TECH_C8ALL for i in NODES));                                                                                                                                                                                                                                       
#@constraint(m, profit-NSCC >= 0);                                                                                                                                                                                                                                                                                                                                                                                      
@objective(m, Max, profit-NSCC);
