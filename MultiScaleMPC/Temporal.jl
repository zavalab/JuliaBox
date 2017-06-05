using JuMP, Gurobi , DataFrames

mkpath("output/Temporal")

env=Gurobi.Env();
mySolver=GurobiSolver(env,OutputFlag=0);
disturbanceFunction(t)=4*sin(2*2*pi*t/T)+sin(12*2*pi*t/T);
disturbanceFunction1(t)=4*sin(2*2*pi*t/T);
disturbanceFunction2(t)=sin(12*2*pi*t/T);
maxIteration=30;
coarseIteration=4;
tolerance=1e-4;
plotFrequency=10;

N=10;
MCoarses=[4];

M=100;
zBar=0;
delta=.1;
T=N*M*delta;
timeGrid=(delta:delta:T);
d=disturbanceFunction(timeGrid);
d1=disturbanceFunction1(timeGrid);
d2=disturbanceFunction2(timeGrid);

# Full Solution
fullProblem = Model(solver=mySolver)
@variable(fullProblem, u[i=1:N*M])
@variable(fullProblem, z[i=0:N*M])

@objective(fullProblem, Min, sum{z[i]^2+u[i]^2,i=1:N*M})
@constraint(fullProblem,fullI, z[0]==zBar)
@constraint(fullProblem, fullD[i=0:N*M-1], z[i+1]==z[i]+(u[i+1]+d[i+1])*delta);


statusFull = solve(fullProblem);
zFull=getvalue(z[1:N*M]);
uFull=getvalue(u);

lambda=zeros(N+1,1);

# GS iteration
zGS=zeros(N*M)
uGS=zeros(N*M)
uCoarse=zeros(N*M)
zCoarse=zeros(N*M)
zInitial=zBar;

PIup=zeros(1,M+1);
PIup[1]=1;
PIdown=zeros(1,M+1);
PIdown[end]=1;

errorSave=[]
errorSaveWC=[]
zSave=[]


for MCoarse in MCoarses
  numPerCoarse=div(M,MCoarse)
  zInitial=zBar
  for k=1:N
    GSProblem=Model(solver=mySolver)
    @variable(GSProblem, zC[j=0:MCoarse])
    @variable(GSProblem, uC[j=1:MCoarse])
    @expression(GSProblem, dC[j=1:MCoarse], sum{d[M*(k-1)+i], i in numPerCoarse*(j-1)+1:numPerCoarse*j})
    @expression(GSProblem, z[i=0:M],zC[fld(i-1,numPerCoarse)+1])
    @expression(GSProblem, u[i=1:M],uC[fld(i-1,numPerCoarse)+1])

    @constraint(GSProblem,GSD[i=0:MCoarse-1], zC[i+1]==zC[i]+uC[i+1]*numPerCoarse*delta+dC[i+1]*delta);
    @constraint(GSProblem,GSI, (z[0]-zInitial)==0)

    @objective(GSProblem, Min, sum{z[i]^2+u[i]^2,i=1:M}+lambda[k+1]*z[M])


    statusGS = solve(GSProblem);
    zTemp=getvalue(z[1:M]);
    uTemp=getvalue(u);
    zGS[M*(k-1)+1:M*k]=copy(zTemp)
    uGS[M*(k-1)+1:M*k]=copy(uTemp)
    lambda[k]=getdual(GSI)
    zInitial=zTemp[end];
  end

  error=norm(zFull-zGS);
  uCoarse=copy(uGS)
  zCoarse=copy(zGS)
  push!(errorSaveWC,error)
  println("error=$error")
end


countIteration=0;
while error>tolerance && countIteration<maxIteration
  countIteration+=1
  zInitial=zBar;
  for k=1:N
      GSProblem=Model(solver=mySolver)
      @variable(GSProblem, u[i=1:M])
      @variable(GSProblem, z[i=0:M])

      @objective(GSProblem, Min, sum{z[i]^2+u[i]^2,i=1:M}+lambda[k+1]*z[M])
      @constraint(GSProblem,GSI, z[0]==zInitial)
      @constraint(GSProblem,GSD[i=0:M-1], z[i+1]==z[i]+(u[i+1]+d[M*(k-1)+i+1])*delta);


    statusGS = solve(GSProblem);
    zTemp=getvalue(z[1:M]);
    uTemp=getvalue(u);
    zGS[M*(k-1)+1:M*k]=copy(zTemp)
    uGS[M*(k-1)+1:M*k]=copy(uTemp)
    lambda[k]=getdual(GSI)
    zInitial=zTemp[end];
  end
  error=norm(zFull-zGS);
  println("error=$error")
  push!(errorSave,error)
  push!(errorSaveWC,error)
  push!(zSave,copy(zGS))
end
DataFrames.writetable("output/Temporal/timeGrid.csv",DataFrame(timeGrid[:,:]),header=false)
DataFrames.writetable("output/Temporal/dFull.csv",DataFrame(d[:,:]),header=false)
DataFrames.writetable("output/Temporal/d1.csv",DataFrame(d1[:,:]),header=false)
DataFrames.writetable("output/Temporal/d2.csv",DataFrame(d2[:,:]),header=false)
DataFrames.writetable("output/Temporal/zFull.csv",DataFrame(zFull[:,:]),header=false)
DataFrames.writetable("output/Temporal/zCoarse.csv",DataFrame(zCoarse[:,:]),header=false)
DataFrames.writetable("output/Temporal/zGS.csv",DataFrame(zGS[:,:]),header=false)
DataFrames.writetable("output/Temporal/errorSave.csv",DataFrame(errorSave[:,:]),header=false,quotemark='\t')
DataFrames.writetable("output/Temporal/errorSaveWC.csv",DataFrame(errorSaveWC[:,:]),header=false,quotemark='\t')

for i=1:countIteration
  DataFrames.writetable("output/Temporal/zSave$(i).csv",DataFrame(zSave[i][:,:]),header=false)
end
for i=countIteration+1:100
  if ispath("output/Temporal/zSave$(i).csv")
    rm("output/Temporal/zSave$(i).csv")
  end
end
