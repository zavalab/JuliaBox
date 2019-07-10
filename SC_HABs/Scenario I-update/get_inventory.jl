# used for get inventory levels

Level = zeros(15,1);
prod = "p3";

inv_results = readdlm("inventory_results.csv",',');
item = inv_results[:,1];

for i in 1:length(item)
    if inv_results[i,2] == prod
        T = inv_results[i,4];
        value = inv_results[i,3];
        Level[T] = Level[T] + value;
    end
end

open("p3_inventory.csv", "w") do pp
for t in 1:15
    println(pp, Level[t])
end
end
