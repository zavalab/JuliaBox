# generate inflow file for protech

open("inflow.csv","w") do pr
    println(pr, "#daily flow (m3)",",","P (mg)",",", "N (mg)");

    nutrient_results = readdlm("NEWS2_results.csv",',');
    runoff      = readdlm("daily_runoff.csv",',');
    TIME    = nutrient_results[:,1];
    N       = nutrient_results[:,2];
    P       = nutrient_results[:,3];

    rainfall_matrix = readdlm("rainfall.csv",',');      # rainfall matrix [mmH2O]
    rainfall = rainfall_matrix[:,2];

    point_emit = 370000/170000*130917*120*0.003785      # population * 120 gallon/cap/day * [gallon -> m3]

    A_ws    = 529.81*10^6                   # watershed area m2
    A_lake  = 39.4*10^6                     # lake area m2
for d in 1:14*length(TIME)
    Q_day = point_emit + A_lake*rainfall[d]*0.001 + runoff[d]*A_ws*0.001;
    res = d%7;
    if res != 0
        weekly_runoff = sum(runoff[(d-(res-1)):(d+(7-res))]);
    else
        weekly_runoff = sum(runoff[(d-6):d]);
    end
    P_day = P[Int(ceil(d/14))]*runoff[d]/weekly_runoff*10^6;
    N_day = N[Int(ceil(d/14))]*runoff[d]/weekly_runoff*10^6;
    println(pr, Q_day, ",", P_day,",", N_day);
end
end
