# Facility Location Case Study
# Created by Alex Dowling (adowling2@wisc.edu)
# University of Wisconsin-Madison
# Department of Chemical and Biological Engineering
# Last Modified Feb. 17th, 2016

using JuMP
using PyPlot
using Gurobi

function calcDistanceMatrix(landmarks,candidates)

	n = size(landmarks,1)
	m = size(candidates,1)
	
	D = zeros(n,m)
	
	for i = 1:n
		for j = 1:m
			D[i,j] = norm(landmarks[i,:] - candidates[j,:])
			
		end
	end
	
	return D

end

# Required by CVaR library
type ProbData
	farmLoc
	cityLoc
	waterLoc
	candLoc
	dFarmCand
	dCityCand
	dWaterCand
	
	function ProbData(nFarms,nUrban,nWater,nCandidate)
		this = new()
		farmLoc = rand(nFarms,2)
		cityLoc = rand(nUrban,2)
		waterLoc = rand(nWater,2)
		candLoc = rand(nCandidate,2)
		
		this.dFarmCand = calcDistanceMatrix(farmLoc,candLoc)
		this.dCityCand = calcDistanceMatrix(cityLoc,candLoc)
		this.dWaterCand = calcDistanceMatrix(waterLoc,candLoc)
	
		this.farmLoc = farmLoc
		this.cityLoc = cityLoc
		this.waterLoc = waterLoc
		this.candLoc = candLoc
	
		return this
	
	end
end

# Required by CVaR library
type Solution

	f
	ySmall
	yLarge
	farmToFacilityNetwork
	minDistUtoC
	minDistWtoC
	z # Either a single float or an array of floats
	alpha
	cvar
end

# This function is used to write data to a file for plotting using Python
function writeDataToFile(d,filename)
	file = open(filename,"w+")
	
	@printf(file,"Farms: %d\n",size(d.farmLoc,1))
	for i = 1:size(d.farmLoc,1)
		@printf(file,"%f,%f\n",d.farmLoc[i,1],d.farmLoc[i,2])
	end
	@printf(file,"\n")
	
	@printf(file,"Cities: %d\n",size(d.cityLoc,1))
	for i = 1:size(d.cityLoc,1)
		@printf(file,"%f,%f\n",d.cityLoc[i,1],d.cityLoc[i,2])
	end
	@printf(file,"\n")
	
	@printf(file,"Water: %d\n",size(d.waterLoc,1))
	for i = 1:size(d.waterLoc,1)
		@printf(file,"%f,%f\n",d.waterLoc[i,1],d.waterLoc[i,2])
	end
	@printf(file,"\n")
	
	@printf(file,"Facilities: %d\n",size(d.candLoc,1))
	for i = 1:size(d.candLoc,1)
		@printf(file,"%f,%f\n",d.candLoc[i,1],d.candLoc[i,2])
	end
	
	close(file)
	
	return 0

end

# This function is used to write data to a file for plotting using Python
function writeSolutionToFile(sln, d, filename)
	file = open(filename,"w+")
	
	nFarm = size(d.farmLoc,1)
	nCand = size(d.candLoc,1)
		
	@printf(file,"Small Facilities:\n")
	for i = 1:nCand
		if(sln.ySmall[i] >= 1 - 1E-6)
			@printf(file,"%d\n",i)
		end
	end
	
	@printf(file,"Large Facilities:\n")
	for i = 1:nCand
		if(sln.yLarge[i] >= 1 - 1E-6)
			@printf(file,"%d\n",i)
		end
	end
	
	@printf(file,"Farm-Facility Matrix (Sparse Form):\n")
	for i = 1:nFarm
		for j = 1:nCand
			if(sln.farmToFacilityNetwork[i,j] > 1E-6)
				@printf(file,"(%d,%d,%f)\n",i,j,sln.farmToFacilityNetwork[i,j])
			end
		end
	end
	
	# Print additional information not used for plotting
	@printf(file,"f = [%f, %f, %f, %f]\n",sln.f[1], sln.f[2], sln.f[3], sln.f[4]);
	
	if(length(sln.z) > 1)
		n = length(sln.z)
		@printf(file,"z = [");
		for i = 1:(n-1)
			@printf(file,"%f, ",sln.z[i])
		end
		@printf(file,"%f]\n",sln.z[end])
	
		@printf(file,"alpha = %f\n",sln.alpha)
		@printf(file,"CVaR = %f\n",sln.cvar)
	else
		@printf(file,"z = %f\n",sln.z)
	end
	
	close(file)
	
	return 0 
end

# Required by CVaR library
function generateBaseModel(d)
# Solve either for a single weight or CVaR formulation

	m = Model(solver=GurobiSolver())
	
	nF = size(d.dFarmCand,1)
	nU = size(d.dCityCand,1)
	nW = size(d.dWaterCand,1)
	nC = size(d.dFarmCand,2)
	
	F = 1:nF
	U = 1:nU
	W = 1:nW
	C = 1:nC
	
	capacityW = [2, 10]
	costW = [2, 8.5]
	waterW = [2, 9.5]
	safetyW = [2, 2]
	
	# Define objective functions
	@defVar(m,f[1:4])
	
	# Existence of facility at each candidate site
	@defVar(m,ySmall[C],Bin)
	@defVar(m,yLarge[C],Bin)
	
	# Distance to nearest facility for each farm
	@defVar(m, 0 <= transCost[F,C] <= 1)
	
	# Smallest distance to facility
	@defVar(m, minDistUtoC[U] >= 0)
	@defVar(m, minDistWtoC[W] >= 0)
	
	#### Compute objectives
	
	# Total distance from each farm to nearest facility
	@addConstraint(m, CalcTotalDistanceToFarms, f[1] == sum{sum{d.dFarmCand[i,j] * transCost[i,j],i=F},j=C})

	
	# Calculate total distance of constructed facilities to each urban area/water source
	# This formulation prefers building all of the facilities
	# @addConstraint(m, CalcTotalDistanceToCities, f[2] == -sum{sum{d.dCityCand[i,j]*(ySmall[j]*safetyW[1] + yLarge[j]*safetyW[2]),i=U},j=C})
	# @addConstraint(m, CalcTotalDistanceToWater, f[3] == -sum{sum{d.dWaterCand[i,j]*(ySmall[j]*waterW[1] + yLarge[j]*waterW[2]),i=W},j=C})
	
	@addConstraint(m, CalcTotalDistanceToCities, f[2] == -sum{minDistUtoC[i],i=U})
	@addConstraint(m, CalcTotalDistanceToWater, f[3] == -sum{minDistWtoC[i],i=W})
	
	# Calculate total cost
	@addConstraint(m, CalcTotalCost, f[4] == sum{ySmall[i]*costW[1] + yLarge[i]*costW[2],i=C})

	
	### Compute transportation costs/network topology
	@addConstraint(m, NetworkIn[i=F], sum{transCost[i,j],j=C} == 1)
	@addConstraint(m, NetworkOut[j=C], sum{transCost[i,j],i=F} <= capacityW[1]*ySmall[j] + capacityW[2]*yLarge[j])
	
	### Compute min distances
	M = sqrt(2)
	
	# Should I use a convex-hull relaxation here?
	@addConstraint(m, CalcMinDistUtoC[i=U, j=C], minDistUtoC[i] <= d.dCityCand[i,j]*(safetyW[1]*ySmall[j] + safetyW[2]*yLarge[j]) + M*safetyW[2]*(1 - ySmall[j] - yLarge[j]))
	@addConstraint(m, CalcMinDistUtoC[i=W, j=C], minDistWtoC[i] <= d.dWaterCand[i,j]*(waterW[1]*ySmall[j] + waterW[2]*yLarge[j]) + M*waterW[2]*(1 - ySmall[j] - yLarge[j]))
	
	### Additional constraints
	
	# Total capacity
	# @addConstraint(m, TotalCapacity, 12 <= sum{ySmall[i]*capacityW[1] + yLarge[i]*capacityW[2],i=C})
		
	# One facility size per location
	@addConstraint(m, OnePerLocation[i=C], 1 >= ySmall[i] + yLarge[i])
	
	return m
end

function assembleSolution(m,alpha)

	f = getVar(m, :f)
	ySmall = getVar(m, :ySmall)
	yLarge = getVar(m, :yLarge)
	transCost = getVar(m, :transCost)
	minDistUtoC = getVar(m, :minDistUtoC)
	minDistWtoC = getVar(m, :minDistWtoC)
	z = getVar(m, :z)
	
	if(0 <= alpha && alpha <= 1)
		cvar = getVar(m, :cvar)
		s = Solution(getValue(f),getValue(ySmall),getValue(yLarge),getValue(transCost),getValue(minDistUtoC),getValue(minDistWtoC),getValue(z),alpha,getValue(cvar))
	else
		s = Solution(getValue(f),getValue(ySmall),getValue(yLarge),getValue(transCost),getValue(minDistUtoC),getValue(minDistWtoC),getValue(z),[],[])
	end

end

include("generic_CVaR_library.jl")

##### Script Part

# Adjust this flag to toggle between the original (pessimistic) and alternate nadir point
flagAltNadir = true

# Seed random number generator to keep results consistent
srand(2015)

nFarms = 12
nUrban = 2
nWater = 5
nCandidate = 30


d = ProbData(nFarms,nUrban,nWater,nCandidate)

# Calculate Utopia and Nadir points
mooData = computeNadirUtopiaPoints(d, flagAltNadir)

# Print solutions to file for plotting with Python script?
flagPlotSingleObjSolutions = false
flagSaveSingleObjectiveSolutions = false

for i = 1:length(mooData.nObj)

	if(flagSaveSingleObjectiveSolutions)
		writeSolutionToFile(mooData.slnSingleObjTrad[i],d,string("Solution_Obj_",i,"_only_trad.txt"));
	end
	
	if(flagSaveSingleObjectiveSolutions)
		writeSolutionToFile(mooData.slnSingleObjAlt[i],d,string("Solution_Obj_",i,"_only_alt.txt"));
	end
	
	
end


nStake = 50
W = rand(nStake,4)

for i = 1:size(W,1)
	W[i,:] = W[i,:] / sum(W[i,:])
end

if(flagAltNadir)
	# Use alternate nadir point
	fileNameBase = "Alt_Nadir"
else
	# Use traditional nadir point
	fileNameBase = "Trad_Nadir"
end

stakeholderData  = solveForPreferredPoints(d,W,mooData)

flagPlotStakeholderPreferredPoints = true

if(flagPlotStakeholderPreferredPoints)
	for i = 1:nStake
		writeSolutionToFile(stakeholderData.slnIdeal[i],d,string("Solution_Stakeholder_",i,"_",fileNameBase,".txt"))
	end
end


alpha = 0:0.01:1

sCVaR = Array(Solution,length(alpha))
cvar = zeros(length(alpha))

for i = 1:length(alpha)
	sCVaR[i] = solveCVaRProblem(d,stakeholderData,alpha[i])
	cvar[i] = sCVaR[i].cvar
	
end

zPrev = zeros(nStake)
for i = 1:length(alpha)
# If any of the stakeholder dissatisfactions change, write the solution to a file
	if(sum(abs(zPrev - sCVaR[i].z) .< 1E-8) < nStake)
		writeSolutionToFile(sCVaR[i],d,string("Solution_CVaR_",i,"_",fileNameBase,".txt"));
		println("Change in solution when alpha = ",alpha[i],", i = ",i)
	end
	
	zPrev = sCVaR[i].z
end


figure()
plot(alpha,cvar)
xlabel(L"\alpha")
ylabel("CVaR")
