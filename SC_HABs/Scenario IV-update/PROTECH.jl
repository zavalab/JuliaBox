## PROTECH/PROTBAS Model to predict algae growth in vertical dimension
## Coded by Yicheng Hu 2018-07

## Read data
algae_matrix  = readdlm("algae_information.csv",',');
temperature   = readdlm("temperature_profile.csv",',');
mixed_layer   = readdlm("mixed_layer.csv",',');
weather       = readdlm("weather.csv",',');
inflow        = readdlm("inflow.csv",',');
initial_value = readdlm("initial_composition.csv",',');

V = 12.5*39.4*10^6; # volume of the waterbody (m3)

global Xed          # edible amount
global lm           # the last layer in mixed layer
## Define sets and dictionaries
algae = algae_matrix[:,1];
day   = 1:length(temperature[1,:]);
layer = 1:length(temperature[:,1]);

s = Dict(zip(algae, algae_matrix[:,2])); # surface area of each algae [μm2]
v = Dict(zip(algae, algae_matrix[:,3])); # cell volume of each algae [μm2]
m = Dict(zip(algae, algae_matrix[:,4])); # maximum cell dimension of each algae [μm]
feature = Dict(zip(algae, algae_matrix[:,5])) # feature of each species (if N fixer)
feature2= Dict(zip(algae, algae_matrix[:,6])) # feature of each species (if edible)

θ = Dict((day[1],layer[1]) => 20.0);  # water temperature in a given day in a given layer [C]
for d in day
    for l in layer
        θ[(d,l)] = temperature[l,d];
    end
end

z = Dict(zip(layer, 0.1*layer - 0.05))   # the depth in the middle of each layer [m]
T = Dict(zip(day, weather[:,2]));        # daylight period [hr]
E = Dict(zip(day, weather[:,3]));        # total electromagnetic radiative flux arriving atmosphere [mol/(m2·s)]
fc= Dict(zip(day, weather[:,4]));        # portion of clear sky, 1 if completely clear
mix = Dict(zip(day, mixed_layer[:,1]));  # the depth of mixed layer [m]

Qinflow = Dict(zip(day, inflow[:,1])); # daily water inflow to the lake [m3]
Pinflow = Dict(zip(day, inflow[:,2])); # daily N inflow to the lake [mg]
Ninflow = Dict(zip(day, inflow[:,3])); # daily P inflow to the lake [mg]

## Algae growth modeling

# for record purpose
X = Dict((algae[1],day[1],layer[1]) => 0.0); # concentration of chl-a in each species in each layer in each day [mg/m3]
P = Dict(zip([0;day], [45; zeros(length(day),1)]));  # concentration of P in lake [mg/m3] with initial value (at end of day, or available the next day)
N = Dict(zip([0;day], [360;zeros(length(day),1)]));  # concentration of N in lake [mg/m3] with initial value
ε = Dict(zip((day),zeros(length(day),1)));     #vertical light extinction coefficient [m-1]
for i in 1:length(algae)
    for l in layer
        X[(algae[i],0,l)] = initial_value[l,i];  # assign initial value of concentration of chl-a
        for d in day
            X[(algae[i],d,l)] = 0.0;
        end
    end
end

# calculate algae related parameters
r20 = Dict(zip(algae, zeros(length(algae),1))); # ideal growth rate at 20 C [day-1]
b   = Dict(zip(algae, zeros(length(algae),1))); # parameter for temperature impact
αr  = Dict(zip(algae, zeros(length(algae),1))); # slope of Ik-rθ
for a in algae
    r20[a] = 1.142*(s[a]/v[a])^0.325;
    b[a]   = 3.378-2.505*log10(s[a]/v[a]);
    αr[a]  = 0.257*(m[a]*s[a]/v[a])^0.236;
end

Glast = Dict((algae[1],layer[1]) => 0.0); #rate of loss due to grazing [mg/m3 day-1] in each layer for each species
glast = Dict((algae[1],layer[1]) => 0.0);
for a in algae
    for l in layer
        Glast[(a,l)] = 0;
        glast[(a,l)] = 0;
    end
end

for d in day
    # calculate daily parameters
    TSS = 5 #[mg/L]
    ε[d] = 0.01*sum(X[a,d-1,l] for a in algae for l in layer) + 1.16 + 0.145*TSS #vertical light extinction coefficient [m-1]
    I0 = 0.47*0.8*E[d]*(0.3+0.7*fc[d]);
    lm = maximum(z[l]*(z[l] <= mix[d]) for l in layer); # the mixed layer index
    #lm = 125;
    # calculate daily variables in differet depth for each algae
    r   = Dict((algae[1],layer[1]) => 0.1);   # final growth rate [day-1]
    rθ  = Dict((algae[1],layer[1]) => 0.1);   # ideal growth rate at θ C [day-1]
    rθI = Dict((algae[1],layer[1]) => 0.1);   # modified growth rate considering light intensity I [day-1]
    rcorθI = Dict((algae[1],layer[1]) => 0.1);# modified growth rate considering respiration [day-1]
    hp  = Dict((algae[1],layer[1]) => 0.1);   # light compensate depth for each algae [m]
    Ik  = Dict((algae[1],layer[1]) => 0.1);   # photon flux necessary to saturate the instantaneous growth rate [mol/(m2·s)]
    for l in layer
        for a in algae
            rθ[(a,l)] = r20[a]*10^(b[a]*(1000/293 - 1000/(273+θ[d,l])));
            Ik[(a,l)]  = rθ[a,l]/αr[a]/3600/24;
            hp[(a,l)]  = log(2*I0/Ik[a,l])/ε[d];
            if z[l] <= hp[a,l]
                rθI[(a,l)] = rθ[a,l]*T[d]/24;
            else
                rθI[(a,l)] = 0.223*αr[a]*I0*3600*24*exp(-ε[d]*z[l]);
            end
            rcorθI[(a,l)] = 1.055*rθI[a,l] - 0.07*rθ[a,l];
        end
    end
    # consider the impact of nutrients
    Pdem  = Dict((layer[1]) => 0.1);
    Ndem  = Dict((layer[1]) => 0.1);
    ratioP  = Dict((layer[1]) => 0.1);
    ratioN= Dict((layer[1]) => 0.1);
    for l in layer
        if l <= lm
            Pdem[l] = 1.2*sum((exp(rcorθI[a,ll])-1)*X[a,d-1,ll]*(rcorθI[a,ll]>=0) for a in algae for ll in 1:lm)/lm; # P demand for algae to finish growth
            Ndem[l] = 8.3*sum((exp(rcorθI[a,ll])-1)*X[a,d-1,ll]*(rcorθI[a,ll]>=0)*(feature[a]!="Yes") for a in algae for ll in 1:lm)/lm; # N demand for algae to finish growth
        else #algae can only use local nutrients
            Pdem[l] = 1.2*sum((exp(rcorθI[a,l])-1)*X[a,d-1,l]*(rcorθI[a,l]>=0) for a in algae); # P demand for algae to finish growth
            Ndem[l] = 8.3*sum((exp(rcorθI[a,l])-1)*X[a,d-1,l]*(rcorθI[a,l]>=0)*(feature[a]!="Yes") for a in algae); # N demand for algae to finish growth
        end
        #print(Pdem[l])
        ratioP[l] = 1; # the ratio of available P over demanded P
        ratioN[l] = 1; # the ratio of available N over demanded N
        if Pdem[l] > (P[d-1]-3)
            if P[d-1] >= 3
                ratioP[l] = (P[d-1]-3)/Pdem[l];
            else
                ratioP[l] = 0;
            end
        end
        if Ndem[l] > (N[d-1]-80)
            if N[d-1] >= 80
                    ratioN[l] = (N[d-1]-80)/Ndem[l];
            else
                ratioN[l] = 0;
            end
        end
        #print(ratioP[l]);
    end
    #println(sum(Pdem[ll] for ll in layer)/length(layer));
    rP   = Dict((algae[1],layer[1]) => 1e20); # growth rate limited by P [day-1]
    rN   = Dict((algae[1],layer[1]) => 1e20); # growth rate limited by N [day-1]
    for l in layer
        for a in algae
             if rcorθI[a,l] >= 0;
                 rP[(a,l)] = ratioP[l] * rcorθI[a,l];
                 if feature[a] != "Yes"
                     rN[(a,l)] = ratioN[l] * rcorθI[a,l];
                 else
                     rN[(a,l)] = 1e20;
                 end
             else
                 rP[(a,l)] = 1e20;
                 rN[(a,l)] = 1e20;
             end
             r[(a,l)] = min(rcorθI[a,l], rP[a,l], rN[a,l]);
        end
    end
    # consider the impact of grazing
    G = Dict((algae[1],layer[1]) => 0.0); #rate of loss due to grazing [mg/m3 day-1] in each layer for each species
    g = Dict((algae[1],layer[1]) => 0.0);
    Xed = sum(X[a,d-1,l]*(feature2[a] == "Yes") for a in algae for l in layer)/length(layer);
    for a in algae
        if feature2[a] == "Yes"
            if Xed <= 1.6
                for l in layer
                    #G[(a,l)] = Glast[a,l] + 0.005*exp(0.2297*θ[d,l]);
                    G[(a,l)] = 0;
                    g[(a,l)] = 0;
                end
            else
                for l in layer
                    if θ[d,l] >= 7
                        if Xed >= 10
                            g[(a,l)] = glast[a,l] + 0.1911*(θ[d,l]/13);
                        else
                            g[(a,l)] = glast[a,l] + 0.1911*(θ[d,l]/13)*(sum(X[a,d-1,l] for l in layer)/8.4);
                        end
                    else
                        g[(a,l)] = glast[a,l]
                    end
                    G[(a,l)] = X[a,d-1,l] * (1-exp(-g[a,l]));
                end
            end
        else
            for l in layer
                G[(a,l)] = 0;
                g[(a,l)] = 0;
            end
        end
    end
    Glast = G;
    glast = g;
    # consider the impact of dilution (assume steady state)
    w = Qinflow[d]/V # dilution ratio
    D = Dict((algae[1],layer[1]) => 0.0);
    for a in algae
        for l in layer
            D[(a,l)] = X[a,d-1,l]*(1-exp(-w));
        end
    end
    # grow
    Xtemp = Dict((algae[1],layer[1]) => 0.0); # temporary record
    R = Dict((algae[1],layer[1]) => 0.0); # grow
    for a in algae
        for l in layer
            R[(a,l)] = (exp(r[a,l])-1)*X[a,d-1,l];
            Xtemp[(a,l)] = X[a,d-1,l] + R[a,l] - G[a,l] - D[a,l];
        end
    end
    # distribute evenly in mixed layer
    for a in algae
        for l in 1:lm
            Xtemp[a,l] = sum(Xtemp[a,l] for l in 1:lm)/lm;
        end
    end
    # algae movement
    I = Dict(zip(layer, zeros(length(layer),1)));
    for l in layer
        I[l] = I0*exp(-ε[d]*z[l]);
    end

    for a in algae
        if a == "Chlorella" || a == "Oocystis"
            for l in layer
                if l == 1;
                    X[a,d,l] = 0;
                elseif l == length(layer)
                    X[a,d,l] = Xtemp[a,l] + Xtemp[a,l-1];
                else
                    X[a,d,l] = Xtemp[a,l-1];
                end
            end
        elseif a == "Asterionella"
            for l in layer
                in = 0;
                out = 0;
                if 1 <= l < length(layer)
                    out = Xtemp[a,l];
                    if l-2 >= 1
                        in = in + Xtemp[a,l-2]*(I[l-2] <= 500);
                    end
                    if  l-10 >= 1
                        in = in + Xtemp[a,l-10]*(I[l-10] > 500);
                    end
                else
                    out = 0;
                    in = Xtemp[a,l-1] + Xtemp[a,l-2] + sum(Xtemp[a,l-i]*(I[l-i] > 500) for i in 3:10);
                end
                X[a,d,l] = Xtemp[a,l] + in - out;
            end
        elseif a == "Anabaena" || a == "Planktothrix"
            for l in layer
                in = 0;
                out = 0;
                if l == 1
                    out = Xtemp[a,l]*(I[l]>30);
                    in = Xtemp[a,l+1]*(I[l+1] <= 10);
                elseif  1 < l < length(layer)
                    out = Xtemp[a,l]*(I[l]<=10 || I[l] >30);
                    in = Xtemp[a,l-1]*(30 < I[l-1] <= 100);
                    if l-3 >= 1
                        in = in + Xtemp[a,l-3]*(I[l-3] > 100);
                    end
                    in = in + Xtemp[a,l+1]*(I[l+1] <= 10);
                else
                    out = Xtemp[a,l]*(I[l] <= 10);
                    in = Xtemp[a,l-1]*(I[l-1] >30) + sum(Xtemp[a,l-i]*(I[l-i] > 100) for i in 2:3);
                end
                X[a,d,l] = Xtemp[a,l] + in - out;
            end
        elseif a == "Microcystis"
            for l in layer
                in = 0;
                out = 0;
                if l == 1
                    out = Xtemp[a,l]*(I[l] > 30);
                    in = sum(Xtemp[a,l+i]*(I[l+i]<=30)*(I[l+i]>5) for i in 1:2) + sum(Xtemp[a,l+i]*(I[l+i]<=5) for i in 3:30);
                elseif 1 < l < length(layer)
                    out = Xtemp[a,l];
                    if l-2 >= 1
                        in = in + Xtemp[a,l-2]*(30 < I[l-2] <= 100);
                    end
                    if  l-6 >= 1
                        in = in + Xtemp[a,l-6]*(I[l-6] > 100);
                    end
                    if l+2 <= length(layer)
                        in = in + Xtemp[a,l+2]*(5 < I[l+2] <= 30);
                    end
                    if  l+30 <= length(layer)
                        in = in + Xtemp[a,l+30]*(I[l+30] <= 100);
                    end
                else
                    out = Xtemp[a,l]*(I[l] <= 30);
                    in = sum(Xtemp[a,l-i]*(I[l-i]<=100)*(I[l-i]>30) for i in 1:2) + sum(Xtemp[a,l-i]*(I[l-i]>100) for i in 3:6);
                end
                X[a,d,l] = Xtemp[a,l] + in - out;
            end
        elseif a == "Ceratium"
            for l in layer
                in = 0;
                out = 0;
                if l == 1;
                    out = Xtemp[a,l]*(I[l] > 100);
                    in = sum(Xtemp[a,l+i]*(I[l+i]<=100) for i in 1:10);
                elseif 1 < l < length(layer)
                    out = Xtemp[a,l];
                    if  l-10 >= 1
                        in = in + Xtemp[a,l-10]*(I[l-10] > 100);
                    end
                    if l+10 <= length(layer)
                        in = in + Xtemp[a,l+10]*(I[l+10] <= 100);
                    end
                else
                    out = Xtemp[a,l]*(I[l] <= 100);
                    in = sum(Xtemp[a,l-i]*(I[l-i]>100) for i in 1:10);
                end
                X[a,d,l] = Xtemp[a,l] + in - out;
            end
        elseif a == "Rhodomonas"
            for l in layer
                in = 0;
                out = 0;
                if l == 1
                    out = Xtemp[a,l]*(I[l] > 100);
                    in = sum(Xtemp[a,l+i]*(I[l+i] <= 30) for i in 1:5);
                elseif 1 < l < length(layer)
                    out = Xtemp[a,l]*(I[l]<=30 || I[l] >100);
                    if  l-5 >= 1
                        in = in + Xtemp[a,l-5]*(I[l-5] > 100);
                    end
                    if l+5 <= length(layer)
                        in = in + Xtemp[a,l+5]*(I[l+5] <= 30);
                    end
                else
                    out = Xtemp[a,l]*(I[l] <= 30);
                    in = sum(Xtemp[a,l-i]*(I[l-i] > 100) for i in 1:5);
                end
                X[a,d,l] = Xtemp[a,l] + in - out;
            end
        end
    end
    # prevent extinction
    for a in algae
        for l in layer
            if X[a,d,l] <= 0
                X[a,d,l] = 1e-10
            end
        end
    end
    # update nutrient environment: nutrient in last period - consumption by algae + released by grazing + inflow
    P[d] = (V*(P[d-1]*exp(-w)  - 1.2*sum((X[a,d,l]-X[a,d-1,l]) for a in algae for l in layer)/length(layer)) + Pinflow[d])/V;
    N[d] = (V*(N[d-1]*exp(-w)  - 8.3*sum((X[a,d,l]-X[a,d-1,l]) for a in algae for l in layer)/length(layer)) + Ninflow[d])/V;
end
y = Dict(zip(day, zeros(length(day),1)));
for d in day
    y[d] = sum(X[a,d,l] for a in algae for l in layer)/length(layer);
end
a0 = Dict(zip(algae, zeros(length(algae),1)));
af = Dict(zip(algae, zeros(length(algae),1)));
for a in algae
    a0[a] = sum(X[a,1,l] for l in layer);
    af[a] = sum(X[a,day[end],l] for l in layer);
end

open("daily_prediction.csv","w") do pr
    for i in 1:length(day)
       println(pr,day[i],",",y[i],",",P[i],",",N[i]);
    end
end

open("daily_algae_composition.csv","w") do al
    print(al,"#day",",")
    for a in algae
        print(al,a,",");
    end
    println(al);
    for d in day
        print(al,d,",");
        for a in algae
            print(al, sum(X[a,d,l] for l in layer)/length(layer), ",");
        end
        println(al);
    end
end

open("daily_chla_distribution.csv","w") do cl
    print(cl,"0",",")
    for l in layer
        print(cl,l,",");
    end
    println(cl);
    for d in day
        print(cl,d,",");
        for l in layer
            print(cl, sum(X[a,d,l] for a in algae), ",");
        end
        println(cl);
    end
end

open("biweek_prediction.csv","w") do pr
    for t in 1:Int(length(day)/14)
        dayin = (t-1)*14+1:(t-1)*14+14;
        println(pr,t,",",1/14*sum(y[d] for d in dayin),",",1/14*sum(P[d] for d in dayin),",",1/14*sum(N[d] for d in dayin));
    end
end

open("biweek_algae_composition.csv","w") do al
    print(al,"#",",")
    for a in algae
        print(al,a,",");
    end
    println(al);
    for t in 1:Int(length(day)/14)
        print(al,t,",");
        for a in algae
            dayin = (t-1)*14+1:(t-1)*14+14;
            print(al, 1/14*sum(X[a,d,l] for l in layer for d in dayin)/length(layer), ",");
        end
        println(al);
    end
end
