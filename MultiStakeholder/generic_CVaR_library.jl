###### GENERIC CVAR LIBRARY
# Created by Alex Dowling
# University of Wisconsin-Madison
# Department of Chemical and Biological Engineering
# Last Feb. 16th, 2016
#
# This file contains functions to solve multi-stakeholder
# multi-objective problems using the CVaR framework
#
#
# This library requests two data types and two functions to be declared elsewhere
#
# type ProbData
#
# # Place input data for the optimization problem here
#
# end
#
#
# type Solution
#
#	f # Array of objective function values
#	z # Either a single float or an array of floats
#	alpha # Either a float or empty	
#	cvar # Either a float or empty
#
# # Place additional solution information here
#  	
# end
#
#
#

# function generateBaseModel(d)
# # Create instance of base module using ProbData object "d"
#
#	m = Model()
#	
#	return m
# end

# function assembleSolution(m,alpha)
#
#	f = getVar(m, :f)
#	z = getVar(m, :z)
#	
# # Extract additional variable values from optimal solution here
#	
#	if(0 <= alpha && alpha <= 1)
#		cvar = getVar(m, :cvar)
#
# # Don't forget to add any additional variables to the Solution object here
#		s = Solution(getValue(f),getValue(z),alpha,getValue(cvar))
#
#	else
#
# # Don't forget to add any additional variables to the Solution object here
#		s = Solution(getValue(f),getValue(z),[],[])
#	end
#
# end

type MOOProbData
	nObj
	fSingleObjTrad
	fSingleObjAlt
	fUtopia
	fNadirTrad
	fNadirAlt
	flagUseAltNadir
end

type StakeholderData
	nStakeholders
	nObj
	W
	zIdeal
	slnIdeal
	fUtopia
	fNadir
end

function solveWithIndividualWeight(d, w, fUtopia, fNadir)
	
	m = generateBaseModel(d)
	
	f = getVar(m, :f)
	nObjectives = length(f)
	O = 1:nObjectives
	
	@defVar(m, z)
	
	if(isempty(fUtopia) && isempty(fNadir))
		# Do not scale objective with Utopia and Nadir points
		@addConstraint(m, OBJECTIVE, z == sum{f[i]*w[i],i=O})
		
	else
		# Scale objective with Utopia and Nadir points
		@addConstraint(m, OBJECTIVE_SCALED, z == sum{w[i]*(f[i] - fUtopia[i])/(fNadir[i] - fUtopia[i]),i=O})
	end
	
	@setObjective(m, Min, z)
	solve(m)
	
	return assembleSolution(m, -1)
end

function solveForNadirPoint(d,iRestrict,zRestrict)

	m = generateBaseModel(d)

	f = getVar(m, :f)
	nObjectives = length(f)
	O = 1:nObjectives

	w = 1/(nObjectives-1)*ones(nObjectives)
	w[iRestrict] = 0

	@defVar(m, z)
	@addConstraint(m, OBJECTIVE, z == sum{f[i]*w[i],i=O})
	@addConstraint(m, RESTRICT_F, f[iRestrict] <= zRestrict)
	
	@setObjective(m, Min, z)
	solve(m)
	
	return assembleSolution(m, -1)
end

function determineUtopiaNadir(fArray)

	# Columns: Objectives
	# Runs: Solutions
	n = size(fArray,1)
	fUtopia = zeros(n)
	fNadir = zeros(n)
	
	for i = 1:n
		fUtopia[i] = minimum(fArray[:,i])
		fNadir[i] = maximum(fArray[:,i])
	end

	return (fUtopia, fNadir)
end

function solveForPreferredPoints(d,W,mooData)

	nS = size(W,1)
	zPreferred = zeros(nS)
	slnVector = Array(Solution,nS)
	
	fUtopia = mooData.fUtopia
	
	if(mooData.flagUseAltNadir)
		fNadir = mooData.fNadirAlt
	else
		fNadir = mooData.fNadirTrad
	end
		
	for i = 1:nS
		w = W[i,:]
		slnVector[i] = solveWithIndividualWeight(d, w, fUtopia, fNadir)
		zPreferred[i] = slnVector[i].z
	end

	return StakeholderData(nS,size(W,2),W,zPreferred,slnVector,fUtopia,fNadir)

end

function plotStakeholderDisatisfaction(sln, zPreferred)

	d = zeros(length(sln.z))

	barh(1:length(d),d)

end

function solveCVaRProblem(d,stakeData,alpha)

	W = stakeData.W
	zStakeholder = stakeData.zIdeal
	fUtopia = stakeData.fUtopia
	fNadir = stakeData.fNadir
	nStake = stakeData.nStakeholders
	nObjectives = stakeData.nObj


	m = generateBaseModel(d)
	
	S = 1:nStake
	O = 1:nObjectives
	
	f = getVar(m, :f)
	@defVar(m, z[1:nStake])
	@addConstraint(m, STAKE_OBJ[i=S], z[i] == sum{W[i,j]*(f[j] - fUtopia[j])/(fNadir[j] - fUtopia[j]), j = O})
	
	@defVar(m, zSlack[i=S] >= 0)
	@defVar(m, nu)
	@addConstraint(m, SLACK_CALCULATION[i=S], zSlack[i] >= z[i] - zStakeholder[i] - nu)
	
	@defVar(m, cvar)
	
	if(alpha < 1)
		@addConstraint(m, CVAR_CALCULATION, cvar == sum{nu + 1/(1-alpha)*zSlack[i], i=S}/nStake)
	else
		@addConstraint(m, CVAR_CALCULATION2[i=S], cvar >= zSlack[i])
		setUpper(nu,0)
	end
	
	@setObjective(m, Min, cvar)
	solve(m)

	return assembleSolution(m, alpha)
end

function computeNadirUtopiaPoints(d, flagUseAltNadir)

	# Determine number of objectives 
	m = generateBaseModel(d)
	f = getVar(m, :f)
	nObj = length(f)

	fSingleObjTrad = zeros(nObj,nObj)
	fSingleObjAlt = zeros(nObj,nObj)

	for i = 1:nObj
		w = zeros(nObj)
		w[i] = 1
	
		s = solveWithIndividualWeight(d,w,[],[])
	
		fSingleObjTrad[i,:] = s.f
	
		sStar = solveForNadirPoint(d,i,s.f[i])

		fSingleObjAlt[i,:] = sStar.f
	
	end

	(fUtopiaTrad, fNadirTrad) = determineUtopiaNadir(fSingleObjTrad)
	(fUtopiaAlt, fNadirAlt) = determineUtopiaNadir(fSingleObjAlt)

	return MOOProbData(nObj, fSingleObjTrad, fSingleObjAlt, fUtopiaTrad, fNadirTrad, fNadirAlt, flagUseAltNadir)

end