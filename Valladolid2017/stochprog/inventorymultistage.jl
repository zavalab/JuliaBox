
using JuMP 
using Ipopt 

# Model parameters
NST =3;
STAGE = collect(1:NST);
SCEN0 = collect(1:1);
SCENK = collect(1:2);

Demd = Dict()

# construct scenario tree
Dem=zeros(length(STAGE),length(SCEN0),length(SCENK),length(SCENK))
  p=zeros(length(SCEN0),length(SCENK),length(SCENK))
for st in STAGE, i in SCEN0, j in SCENK, k in SCENK
    p[i,j,k]=1/4;
    if st == 1 && i == 1
    Dem[st,i,j,k]=0.0
    end
    if st == 2 && i == 1 && j == 1
    Dem[st,i,j,k]=5.0
    end
    if st == 2 && i == 1 && j == 2
    Dem[st,i,j,k]=25.0
    end
    if st == 3 && i == 1 && j == 1 && k == 1
    Dem[st,i,j,k]=0.0
    end
    if st == 3 && i == 1 && j == 1 && k == 2
    Dem[st,i,j,k]=25.0
    end
    if st == 3 && i == 1 && j == 2 && k == 1
    Dem[st,i,j,k]=50.0
    end   
    if st == 3 && i == 1 && j == 2 && k == 2
    Dem[st,i,j,k]=25.0
    end       
end

# Create Model (standard)
m = Model(solver=IpoptSolver(print_level=0))

@variable(m, 0<=Inv[STAGE,SCEN0,SCENK,SCENK] <=1000)    
@variable(m, 0<=Buy[STAGE,SCEN0,SCENK,SCENK] <=25)

@objective(m, Min,  sum{p[i,j,k]*Buy[st,i,j,k], st in STAGE, i in SCEN0, j in SCENK, k in SCENK}
                  + sum{ p[i,j,k]*Inv[1,i,j,k],              i in SCEN0, j in SCENK, k in SCENK})

@constraint(m, consinv[st in 1:NST-1,i in SCEN0,j in SCENK,k in SCENK], 
               Inv[st+1,i,j,k] == Inv[st,i,j,k] + Buy[st,i,j,k] - Dem[st,i,j,k])

# non-anticipativity stage 1
@constraint(m,  nonantInv1[st in STAGE,i in SCEN0,j in SCENK,k in SCENK; st==1 && i==1],Inv[st,i,j,k] == Inv[st,1,1,1])
@constraint(m,  nonantBuy1[st in STAGE,i in SCEN0,j in SCENK,k in SCENK; st==1 && i==1],Buy[st,i,j,k] == Buy[st,1,1,1])

# non-anticipativity stage 2
@constraint(m, nonantBuy21[st in STAGE,i in SCEN0,j in SCENK,k in SCENK; st==2 && i==1 && j==1],Buy[st,i,j,k] == Buy[st,1,1,1])
@constraint(m, nonantBuy22[st in STAGE,i in SCEN0,j in SCENK,k in SCENK; st==2 && i==1 && j==2],Buy[st,i,j,k] == Buy[st,1,2,1])

# terminal constraint
@constraint(m,    finalInv[ i in SCEN0,j in SCENK,k in SCENK],Inv[NST,i,j,k]==100)

solve(m)

# Results
println(getvalue(Buy))
println(getvalue(Inv))
println("obj ",getobjectivevalue(m))


