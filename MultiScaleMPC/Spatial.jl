using JuMP, Gurobi, DataFrames
include("Spatial_sub.jl")
mkpath("output/Spatial")

# Choose Solver
env=Gurobi.Env()
mySolver=GurobiSolver(env,OutputFlag=0)

# User Defined Variables
errorTolerance=1e-10
gridSize=100
partitionSize=10
heatCoef=1 #Heat transfer coefficient
maxIteration=20
weightFactor=1
coarseSize=5

# Set parameters
partitionPerSide=div(gridSize,partitionSize)
totNumGrid=gridSize^2;
totNumPartitions=div(totNumGrid,partitionSize^2);
partitionNumGrid=partitionSize^2

# Set index
totalIndex,neighborIndex,innerIndex,neighbor=SetIndexSet(gridSize)

# Set Partition Index
totalPartitionIndex,totalIndexPart,neighborIndexPart,boundaryIndexPart,innerIndexPart,partitionTransfer=SetPartitionIndexSet(partitionPerSide,partitionSize)

# Set Coarse index
coarseIndexPart,coarseTotalIndex,coarseInnerIndex,coarseNeighborIndex,coarseTransfer=SetCoarseIndexSet(coarseSize,partitionSize,totalPartitionIndex,totalIndexPart)


# Set disturbance and setpoint
d1=Dict(ind=>8*sin(4*pi*ind[1]/(gridSize+1))*sin(-4*pi*ind[2]/(gridSize+1)) for ind in innerIndex)
d2=Dict(ind=>0 for ind in innerIndex)
d3=Dict(ind=>sin(-32*pi*ind[1]/(gridSize+1))*sin(32*pi*ind[2]/(gridSize+1)) for ind in innerIndex)

d4=Dict(ind=>60*exp(-((ind[1]-(gridSize)*0.6)^2+(ind[2]-(gridSize)*0.25)^2)/10) for ind in innerIndex)
d5=Dict(ind=>-0*exp(-((ind[1]-(gridSize)*0.75)^2+(ind[2]-(gridSize)*0.1)^2)/5) for ind in innerIndex)
d6=Dict(ind=>-50*exp(-((ind[1]-(gridSize)*0.25)^2+(ind[2]-(gridSize)*0.5)^2)/5) for ind in innerIndex)
d7=Dict(ind=>40*exp(-((ind[1]-(gridSize)*0.7)^2+(ind[2]-(gridSize)*0.7)^2)/10) for ind in innerIndex)
d8=Dict(ind=>-40*exp(-((ind[1]-(gridSize)*0.2)^2+(ind[2]-(gridSize)*0.1)^2)/3) for ind in innerIndex)

dWave=Dict(ind=>d1[ind]+d2[ind]+d3[ind] for ind in innerIndex)
dPeak=Dict(ind=>d4[ind]+d5[ind]+d6[ind]+d7[ind]+d8[ind] for ind in innerIndex)

d=Dict(ind=>d1[ind]+d2[ind]+d3[ind]+d4[ind]+d5[ind]+d6[ind]+d7[ind]+d8[ind] for ind in innerIndex)
zset=Dict(ind=>0 for ind in innerIndex)

# Solve Full Problem
fullProblem = Model(solver=mySolver)

@variable(fullProblem, z[ind in totalIndex]) # State (Potential, Temperature)
@variable(fullProblem, u[ind in innerIndex]) # Input (Forcing, heat)

@objective(fullProblem, Min, sum(weightFactor*(z[ind]-zset[ind])^2+u[ind]^2 for ind in innerIndex))
@constraint(fullProblem, fullBalance[ind in innerIndex], u[ind]+d[ind]-heatCoef*(4*z[ind]-sum(z[ind2] for ind2 in neighbor[ind]))==0)
@constraint(fullProblem, fullBoundary[ind in neighborIndex], z[ind]==0)

statusFull=solve(fullProblem)

zFull=zeros(gridSize,gridSize)
dFull=zeros(gridSize,gridSize)
dWaveFull=zeros(gridSize,gridSize)
dPeakFull=zeros(gridSize,gridSize)
zsetFull=zeros(gridSize,gridSize)
for ind in innerIndex
  zFull[ind[1],ind[2]]=getvalue(z[ind])
  dFull[ind[1],ind[2]]=d[ind]
  dWaveFull[ind[1],ind[2]]=dWave[ind]
  dPeakFull[ind[1],ind[2]]=dPeak[ind]
  zsetFull[ind[1],ind[2]]=zset[ind]
end



zOld=Dict(ind=>0.0 for ind in totalIndex)

dualBoundary=Dict(ind=>0.0 for ind in totalIndex)
dualBoundarySave=Dict(ind=>0.0 for ind in totalIndex)

error=Inf
errorSave=[]
zGS=zeros(gridSize,gridSize)
zSave=[]
countIteration=0
iterOrder=SetOrder1()

# Solve Coarse Problem
for part in iterOrder
    CoarseProblem = Model(solver=mySolver)
    @variable(CoarseProblem, zCoarse[coarse in coarseTotalIndex])

    @expression(CoarseProblem, z[ind in totalIndexPart[part]], zCoarse[coarseTransfer[part][ind]]/sqrt(coarseSize))
    @expression(CoarseProblem, u[ind in innerIndexPart[part]], -d[ind]+heatCoef*(4*z[ind]-sum(z[ind2] for ind2 in neighbor[ind]))) # Input (Forcing, heat)

    @objective(CoarseProblem, Min, sum(weightFactor*(z[ind]-zset[ind])^2+u[ind]^2 for ind in innerIndexPart[part])+sum(z[ind]*dualBoundary[ind] for ind in boundaryIndexPart[part]))
    @constraint(CoarseProblem , CoarseBoundary[coarse in coarseNeighborIndex], sum(z[ind]-zOld[ind] for ind in coarseIndexPart[part][coarse])==0)

    statusCoarse=solve(CoarseProblem)

    for ind in innerIndexPart[part]
      zOld[ind]=getvalue(z[ind])
      zGS[ind[1],ind[2]]=zOld[ind]
      dualBoundarySave[ind]=dualBoundary[ind]
      dualBoundary[ind]=0
    end
    for coarse in coarseNeighborIndex
        for ind in coarseIndexPart[part][coarse]
          dualBoundary[ind]=dualBoundary[ind]+getdual(CoarseBoundary[coarse])
        end
      end
end

for ind in totalIndex
  if dualBoundary==0
    dualBoundary[ind]=dualBoundarySave[ind]
  end
  dualBoundary[ind]=dualBoundary[ind]
end
zCoarse=copy(zGS)

# Solve GS Problem
while error>errorTolerance && countIteration<maxIteration
  for part in iterOrder
      GSProblem = Model(solver=mySolver)
      @variable(GSProblem, z[ind in totalIndexPart[part]]) # State (Potential, Temperature)
      @variable(GSProblem, u[ind in innerIndexPart[part]]) # Input (Forcing, heat)

      @objective(GSProblem, Min, sum(weightFactor*(z[ind]-zset[ind])^2+u[ind]^2 for ind in innerIndexPart[part])+sum(z[ind]*dualBoundary[ind] for ind in boundaryIndexPart[part]))
      @constraint(GSProblem, GSBalance[ind in innerIndexPart[part]], u[ind]+d[ind]-heatCoef*(4*z[ind]-sum(z[ind2] for ind2 in neighbor[ind]))==0)
      @constraint(GSProblem , GSBoundary[ind in neighborIndexPart[part]], z[ind]==zOld[ind])

      statusGS=solve(GSProblem)

      for ind in innerIndexPart[part]
        zOld[ind]=getvalue(z[ind])
        zGS[ind[1],ind[2]]=zOld[ind]
        dualBoundary[ind]=0
      end
      for ind in neighborIndexPart[part]
        dualBoundary[ind]=dualBoundary[ind]+getdual(GSBoundary[ind])
      end
  end
  error=sum((zGS-zFull).^2)
  push!(errorSave,error)
  countIteration+=1
  push!(zSave,copy(zGS))
  println("errror=$error")
end

DataFrames.writetable("output/Spatial/dFull.csv",DataFrame(dFull),header=false)
DataFrames.writetable("output/Spatial/dWaveFull.csv",DataFrame(dWaveFull),header=false)
DataFrames.writetable("output/Spatial/dPeakFull.csv",DataFrame(dPeakFull),header=false)
DataFrames.writetable("output/Spatial/zFull.csv",DataFrame(zFull),header=false)
DataFrames.writetable("output/Spatial/zCoarse.csv",DataFrame(zCoarse),header=false)
DataFrames.writetable("output/Spatial/zGS.csv",DataFrame(zGS),header=false)
DataFrames.writetable("output/Spatial/zsetFull.csv",DataFrame(zsetFull),header=false)
DataFrames.writetable("output/Spatial/errorSave.csv",DataFrame(reshape(errorSave,countIteration,1)),header=false,quotemark='\t')

for i=1:countIteration
  DataFrames.writetable("output/Spatial/zSave$(i).csv",DataFrame(zSave[i]),header=false)
end
for i=countIteration+1:maxIteration
  if ispath("output/Spatial/zSave$(i).csv")
    rm("output/Spatial/zSave$(i).csv")
  end
end
