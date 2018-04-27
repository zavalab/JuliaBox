using DataFrames

# Define type ExpData
type ExpData
  nameExp::AbstractString #String # Name of experiment

  numMicrobe::Int # Number of microbiobe that is involved in this experiment
  species # Species number (e.g. if BH and BU is involved, it saves [1 3])

  timeLength::Int # The length of time grid
  timeLengthExtended::Int # The length of time grid with further discretization

  timePoints # Time grid
  timePointsExtended # Extended ime grid (with further discretization)
  timeDiff # Time difference between extended timegrid

  abundance # Absolute abundance data
  IC # Initial Condition

  numOfDisc # Number of discretization that is used
  weight # Weight for each point
end

# This function sets object for ExpData
function SetExpData(nameExp,species,timePoints,abundance,IC,numOfDisc,weightComp)
  timeLength=length(timePoints)
  timePointsExtended=zeros((timeLength-1)*numOfDisc+1)
  for t in 1:timeLength-1
    timePointsExtended[(t-1)*numOfDisc+1:t*numOfDisc+1]=linspace(timePoints[t],timePoints[t+1],numOfDisc+1)
  end
  weight=ones(timeLength)*weightComp
  if IC=="noIC"
    weight[1]=0
  end

  return ExpData(nameExp,length(species),species,
    timeLength,length(timePointsExtended),
    timePoints,timePointsExtended,diff(timePointsExtended),
    abundance,IC,
    numOfDisc,weight)
end

# This function returns the species number from given species name
function GetSpecies(speciesName,speciesOrder)
  lengthSpeciesOrder=size(speciesOrder,1)
  lengthSpeciesName=size(speciesName,1)

  output=Array{Vector{Int64}}(lengthSpeciesName)

  for i=1:lengthSpeciesName
    temp=speciesName[i,:]
    count=0
    for p=1:length(temp)
      if temp[p]!=""
        count+=1
      end
    end
    numOfMicrobe=fld(count,2)
    output[i]=Vector{Int64}(numOfMicrobe)

    for j=1:lengthSpeciesOrder
      for k=1:numOfMicrobe
        if speciesOrder[j,:]==speciesName[i,2*k-1:2*k]
          output[i][k]=j
        end
      end
    end
  end

  return output
end

# This function returns absolute abundance from the given total abundance and
# relative abundance, time sequence

function GetAbsAbundance(timeseq,relabund,growth0,growth1,growth2,delta)

t0 = timeseq[1:3]
t1 = timeseq[3:5]-timeseq[3]
t2 = timeseq[5:end]-timeseq[5]

numMicrobe=size(relabund,2)

if isnan(t2[end])
  t2=t2[1:end-1]
  timeseq=timeseq[1:end-1]
  relabund=relabund[1:end-1,:]
end
totBiomass0 = Interpolation([0:(length(growth0)-1);]*delta,growth0,t0[1:3]);
totBiomass1 = Interpolation([0:(length(growth1)-1);]*delta,growth1,t1[2:end]);
totBiomass2 = Interpolation([0:(length(growth2)-1);]*delta,growth2,t2[2:end]);

totBiomass1 = [totBiomass0[end]/20;totBiomass1]
totBiomass2 = [totBiomass1[end]/20;totBiomass2]

absabund0 = Array{Vector}(numMicrobe)
absabund1 = Array{Vector}(numMicrobe)
absabund2 = Array{Vector}(numMicrobe)

for i=1:numMicrobe
  absabund0[i] = relabund[1:3,i].*totBiomass0
  absabund1[i] = relabund[3:5,i].*totBiomass1
  absabund2[i] = relabund[5:end,i].*totBiomass2
end

return absabund0, absabund1, absabund2, t0, t1, t2

end

# This function interpolates
function Interpolation(x,y,x0)
  y0=zeros(length(x0))
  for j=1:length(x0)
    for i=1:length(x)-1
      if x[i]<= x0[j] <x[i+1]
        y0[j]= y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(x0[j]-x[i])
      elseif x0[j]>=x[i+1]
        y0[j]=y[i+1]
      end
    end
  end
  return y0
end

# Here we start data processing ------------------------------------------------
# Sets empty array/sets that will be used to store data
data=[]
link=Set() # This will be used to link phase 1-phase 2-phase 3
dataIC=Set() # Saves index for phase 1 data
dataNoIC=Set() # Saves index for phase 2 data

# Experiment "Mono" ------------------------------------------------------------
# Species name
singlen0=readcsv("input/Mono/singlen0.csv")
singlen1=readcsv("input/Mono/singlen1.csv")
singlen2=readcsv("input/Mono/singlen2.csv")

# Time sequence
time0=readcsv("input/Mono/time0.csv")
time1=readcsv("input/Mono/time1.csv")
time2=readcsv("input/Mono/time2.csv")

# Total abundance
singlegrowthn0=readcsv("input/Mono/singlegrowthn0.csv")
singlegrowthn1=readcsv("input/Mono/singlegrowthn1.csv")
singlegrowthn2=readcsv("input/Mono/singlegrowthn2.csv")

numDataMono0=size(singlegrowthn0,1)

if includeMono==true
  speciesMono0=GetSpecies(singlen0,speciesOrder)
  speciesMono1=GetSpecies(singlen1,speciesOrder)
  speciesMono2=GetSpecies(singlen2,speciesOrder)

  ICmono0=zeros(numOfSpecies)
  for i=1:numOfSpecies
    k=1
    while singlegrowthn0[i,k]==0
      k+=1
    end
    ICmono0[i]=singlegrowthn0[i,k]
  end

  for i=1:numDataMono0
    abund0=[]
    push!(abund0,singlegrowthn0[i,:])
    push!(data,SetExpData("mono",speciesMono0[i],time0,
      abund0,ICmono0[i],numOfDiscMono,1))
    push!(dataIC,length(data))
  end
end

# Experiment "Pair"-------------------------------------------------------------
# Species name
uniquename=readcsv("input/Pair/uniquename.csv")

# Time sequence
TM=readcsv("input/Pair/TM.csv")

# Total abundance
pairgrowthn0=readcsv("input/Pair/pairgrowthn0.csv")
pairgrowthn1=readcsv("input/Pair/pairgrowthn1.csv")
pairgrowthn2=readcsv("input/Pair/pairgrowthn2.csv")

# Relative abundance
ABM=readcsv("input/Pair/ABM.csv")

# Time gap between total abundance data points
deltaPair=0.5

numDataPair=size(pairgrowthn0,1)
if includePair==true

  # Setting Initial condition
  ICpair=[0.02*ABM[:,1] 0.02*(1-ABM[:,1])]

  ABM=ABM'
  speciesPair=GetSpecies(uniquename,speciesOrder)
  for i=1:numDataPair
    (absAbundance0, absAbundance1, absAbundance2, t0, t1, t2)=
      GetAbsAbundance(TM[i,:],[ABM[:,i] 1-ABM[:,i]],
      pairgrowthn0[i,:],pairgrowthn1[i,:],pairgrowthn2[i,:],deltaPair
      )
    l=length(t0)+length(t1)+length(t2)-2
    push!(data,SetExpData("Pair",speciesPair[i],t0,
      copy(absAbundance0),ICpair[i,:],numOfDiscMulti,1))
    ind0=length(data)
    push!(dataIC,ind0)

    push!(data,SetExpData("Pair",speciesPair[i],t1,
      copy(absAbundance1),"noIC",numOfDiscMulti,1))
    ind1=length(data)
    push!(dataNoIC,ind1)

    push!(data,SetExpData("Pair",speciesPair[i],t2,
      copy(absAbundance2),"noIC",numOfDiscMulti,1))
    ind2=length(data)
    push!(dataNoIC,ind2)

    push!(link,(ind0,ind1,0.05))
    push!(link,(ind1,ind2,0.05))
  end
end

# PW22 Data---------------------------------------------------------------------
# Species name
relnames=readcsv("input/PW22/relnames.csv")

# Time sequence
PW22time=readcsv("input/PW22/PW22time.csv")

# Total abundance
absdata0=readcsv("input/PW22/absdata0.csv")
absdata1=readcsv("input/PW22/absdata1.csv")
absdata2=readcsv("input/PW22/absdata2.csv")

# Relative abundance
relabm=readcsv("input/PW22/relabm.csv")

# Time gap between total abundance data points
deltaPW22=0.5


numDataPW22=size(absdata0,1)
if includePW22==true
  ICPW22=zeros(numDataPW22,2)

  # Setting Initial condition
  for i=1:numDataPW22
    ICPW22[i,1]=0.0158*relabm[2*i-1,1]
    ICPW22[i,2]=0.0158*relabm[2*i,1]
  end

  speciesPW22=GetSpecies(relnames,speciesOrder)
  for i=1:numDataPW22
    (absAbundance0, absAbundance1, absAbundance2, t0, t1, t2)=
      GetAbsAbundance(PW22time,relabm[2*i-1:2*i,:]',
      absdata0[i,:],absdata1[i,:],absdata2[i,:],deltaPW22)
      l=length(t0)+length(t1)+length(t2)-2

    push!(data,SetExpData("PW22",speciesPW22[i],t0,
      copy(absAbundance0),ICPW22[i,:],numOfDiscMulti,1))
    ind0=length(data)
    push!(dataIC,ind0)

    push!(data,SetExpData("PW22",speciesPW22[i],t1,
      copy(absAbundance1),"noIC",numOfDiscMulti,1))
    ind1=length(data)
    push!(dataNoIC,ind1)

    push!(data,SetExpData("PW22",speciesPW22[i],t2,
      copy(absAbundance2),"noIC",numOfDiscMulti,1))
    ind2=length(data)
    push!(dataNoIC,ind2)

    push!(link,(ind0,ind1,0.05))
    push!(link,(ind1,ind2,0.05))
  end
end

# TM1 Data----------------------------------------------------------------------
# Species name
TMname=readcsv("input/TM1/TMname.csv")

# Time sequence
TM1time=readcsv("input/TM1/TM1time.csv")

# Total abundance
means_0=readcsv("input/TM1/means_0.csv")
means_1=readcsv("input/TM1/means_1.csv")
means_2=readcsv("input/TM1/means_2.csv")

# Relative abundance
TM1rel=readcsv("input/TM1/TM1rel.csv")

# Time gap between total abundance data points
deltaTM1=0.5

numDataTM1=size(means_0,1)
if includeTM1==true
  ICTM1=zeros(numDataTM1,3)

  # Setting Initial condition
  for i=1:numDataTM1
    ICTM1[i,1]=0.03*TM1rel[1,3*i-2]
    ICTM1[i,2]=0.03*TM1rel[1,3*i-1]
    ICTM1[i,3]=0.03*TM1rel[1,3*i]
  end

  speciesTM1=GetSpecies(TMname,speciesOrder)
  for i=1:numDataTM1
    (absAbundance0, absAbundance1, absAbundance2, t0, t1, t2)=
      GetAbsAbundance(TM1time,TM1rel[:,3*i-2:3*i],
      means_0[i,:],means_1[i,:],means_2[i,:],deltaTM1)
      l=length(t0)+length(t1)+length(t2)-2

    push!(data,SetExpData("TM1",speciesTM1[i],t0,
      copy(absAbundance0),ICTM1[i,:],numOfDiscMulti,1))
    ind0=length(data)
    push!(dataIC,ind0)

    push!(data,SetExpData("TM1",speciesTM1[i],t1,
      copy(absAbundance1),"noIC",numOfDiscMulti,1))
    ind1=length(data)
    push!(dataNoIC,ind1)

    push!(data,SetExpData("TM1",speciesTM1[i],t2,
      copy(absAbundance2),"noIC",numOfDiscMulti,1))
    ind2=length(data)
    push!(dataNoIC,ind2)

    push!(link,(ind0,ind1,0.05))
    push!(link,(ind1,ind2,0.05))
  end
end

# L1O6 Data---------------------------------------------------------------------
# Species name
L1O6comm=readcsv("input/L1O6/L1O6comm.csv")

# Time sequence
tL1O6 = [0 780 1425 2280 2910 3720 4310]./60

# Total abundance
hgc0=readcsv("input/L1O6/hgc0.csv")
hgc1=readcsv("input/L1O6/hgc0.csv")
hgc2=readcsv("input/L1O6/hgc0.csv")

# Relative abundance
L1O6r=readcsv("input/L1O6/L1O6r.csv")

# Time gap between total abundance data points
deltaL1O6=1


numDataL1O6=size(hgc0,1)
if includeL1O6==true
  speciesL1O6=GetSpecies(L1O6comm,speciesOrder)

  ICL1O6=Array{Vector}(numDataL1O6)

  # Setting Initial condition
  for i=1:numDataL1O6
    ICL1O6[i]=0.01*ones(length(speciesL1O6[i]))
  end

  for i=1:numDataL1O6
    numMicrobe=length(speciesL1O6[i])
    (absAbundance0, absAbundance1, absAbundance2, t0, t1, t2)=
      GetAbsAbundance(tL1O6,L1O6r[:,numMicrobe*(i-1)+1:numMicrobe*i],
      hgc0[i,:],hgc1[i,:],hgc2[i,:],deltaL1O6)
    l=length(t0)+length(t1)+length(t2)-2

    push!(data,SetExpData("L1O6",speciesL1O6[i],t0,
      copy(absAbundance0),ICL1O6[i],numOfDiscMulti,1))
    ind0=length(data)
    push!(dataIC,ind0)

    push!(data,SetExpData("L1O6",speciesL1O6[i],t1,
      copy(absAbundance1),"noIC",numOfDiscMulti,1))
    ind1=length(data)
    push!(dataNoIC,ind1)

    push!(data,SetExpData("L1O6",speciesL1O6[i],t2,
      copy(absAbundance2),"noIC",numOfDiscMulti,1))
    ind2=length(data)
    push!(dataNoIC,ind2)

    push!(link,(ind0,ind1,0.05))
    push!(link,(ind1,ind2,0.05))
  end
end


numOfData=length(data) # Total number of objects "data"
numExp=includeMono*1+
  includePair*numDataPair+
  includePW22*numDataPW22+
  includeTM1*numDataTM1+
  includeL1O6*numDataL1O6 # This will be used in objective function

ubMat=ones(numOfSpecies,numOfSpecies)*upperBoundA
lbMat=ones(numOfSpecies,numOfSpecies)*lowerBoundA
for i=1:numOfSpecies
  ubMat[i,i]=upperBoundR2
  lbMat[i,i]=lowerBoundR2
end

# param0=readcsv("output/param6.csv")
