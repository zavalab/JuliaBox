# Example 1
# Created by Alex Dowling (adowling2@wisc.edu)
# University of Wisconsin-Madison
# Department of Chemical and Biological Engineering
# Last Modified February 17th, 2016

using JuMP
using PyPlot
using Gurobi


type ProbData
	xUpperBound
	A
	b
	
	ProbData() = new(5,[1 2],5)
end

type Solution

	f # Array of objective function values
	z # Either a single float or an array of floats
	alpha # Either a float or empty	
	cvar # Either a float or empty

	x

end

function generateBaseModel(d)
# Create instance of base module using ProbData object "d"

	m = Model()

	@defVar(m, 0 <= x[1:2] <= d.xUpperBound)
	@defVar(m, f[1:2])
	@addConstraint(m, SINGLE_CONSTRAINT, d.A[1]*x[1] + d.A[2]*x[2] >= d.b)
	@addConstraint(m, OBJECTIVES[i=1:2], x[i] == f[i])
	
	return m
	
 end

function assembleSolution(m,alpha)

	f = getVar(m, :f)
	z = getVar(m, :z)

	x = getVar(m, :x)
		
	if(0 <= alpha && alpha <= 1)
		cvar = getVar(m, :cvar)
		s = Solution(getValue(f),getValue(z),alpha,getValue(cvar),getValue(x))
	else
		s = Solution(getValue(f),getValue(z),[],[],getValue(x))
	end

end

include("generic_CVaR_library.jl")

function solveWithRandW(d,mooData)
	nObj = mooData.nObj
	
	# Adjust the number of stakeholders here
	nStake = 10
	W = rand(nStake,nObj)

	for i = 1:size(W,1)
		W[i,:] = W[i,:] / sum(W[i,:])
	end

	stakeholderData = solveForPreferredPoints(d,W,mooData)

	# Iterate over alpha values and save results
	alpha = 0:0.01:1
	X = zeros(length(alpha),2)
	Z = X

	for i = 1:length(alpha)
		s = solveCVaRProblem(d,W,stakeholderData,alpha[i])
		for j = 1:2
			X[i,j] = s.x[j]
		end
	end

	return (alpha, X, W)
end

function printSolutionChanges(alpha,X)
	indDiff = find(sum(abs(diff(X)),2) .>= 1E-6)
	println("**********")
	println("alpha = ",alpha[1], " x = ",X[1,:])
	for i = 1:length(indDiff)

		println("...")
		println("alpha = ",alpha[indDiff[i]], " x = ",X[indDiff[i],:])
		println("alpha = ",alpha[indDiff[i]+1], " x = ",X[indDiff[i]+1,:])

	end

	println("**********")
end

### Script part

# Determine Utopia and Nadir points
d = ProbData()

# Input 2: Number of objectives
# Input 3: Use alternate nadir point (true/false)
mooData = computeNadirUtopiaPoints(d, true)

srand(2015)


W = eye(2)*0.8 + ones(2,2)*0.1
m = size(W,2)

stakeholderData = solveForPreferredPoints(d,W,mooData)

# Iterate over alpha values and save results
alpha = 0:0.01:1
X = zeros(length(alpha),2)
D = zeros(length(alpha),m)
cvar = zeros(length(alpha))
for i = 1:length(alpha)
	s = solveCVaRProblem(d,stakeholderData,alpha[i])
	for j = 1:2
		X[i,j] = s.x[j]
	end
	for j = 1:size(W,1)
		D[i,j] = s.z[j] - stakeholderData.zIdeal[j]
	end
	cvar[i] = s.cvar
end

plotFlag = true

if(plotFlag)
	figure()
	
	subplot(211)
	
	grid()
	plot(alpha,X[:,1],linestyle="-",color="blue",label=L"$x_0$",linewidth=3,drawstyle="steps-mid")
	plot(alpha,X[:,2],linestyle="--",color="red",label=L"$x_1$",linewidth=3,drawstyle="steps-mid")
	ylabel("Variables",fontsize=18)
	legend(fancybox="True", shadow="True",ncol=2)

	subplot(212)

	grid()
	plot(alpha,D[:,1],linestyle="-",color="blue",label="Stakeholder 0",linewidth=3,drawstyle="steps-mid")
	plot(alpha,D[:,2],linestyle="--",color="red",label="Stakeholder 1",linewidth=3,drawstyle="steps-mid")
	xlabel(L"$\alpha$",fontsize=18)
	ylabel("Dissatisfaction",fontsize=18)
	legend(fancybox="True", shadow="True",loc=4)

figure()
grid()
plot(alpha,cvar)
xlabel(L"$\alpha$",fontsize=18)
ylabel(L"$CVaR_{\alpha}$",fontsize=18)
end


# This section generates plots that were not included in the final version of the manuscript
#=
(alpha1, X1, W1) = solveWithRandW(d,fUtopiaUsed,fNadirUsed)
(alpha2, X2, W2) = solveWithRandW(d,fUtopiaUsed,fNadirUsed)
(alpha3, X3, W3) = solveWithRandW(d,fUtopiaUsed,fNadirUsed)
(alpha4, X4, W4) = solveWithRandW(d,fUtopiaUsed,fNadirUsed)

#println("Analyzing Set 1")
#printSolutionChanges(alpha1,X1)

#println("Analyzing Set 2")
#printSolutionChanges(alpha2,X2)

#println("Analyzing Set 3")
#printSolutionChanges(alpha3,X3)

#println("Analyzing Set 4")
#printSolutionChanges(alpha4,X4)

figure()

x = [5; 0]
y = [0; 5/2]

#scatter(x,y,color="black",marker="*",s=100.0)
plot(x,y,color="black",marker="*",markersize=10.0,label="Pareto Surface",linestyle="--")

plot(X1[:,1], X1[:,2],color="r", marker="o",label="Set A")
plot(X2[:,1], X2[:,2],color="b", marker="^",label="Set B")
plot(X3[:,1], X3[:,2],color="g", marker="v",label="Set C")
plot(X4[:,1], X4[:,2],color="m", marker="s",label="Set D")


legend(fancybox="True", shadow="True")

xlim([-0.5, 5.5])
ylim([-0.5, 3])

xlabel(L"x_1",fontsize=18)
ylabel(L"x_2",fontsize=18)

figure()
plot3D(X1[:,1], X1[:,2], alpha1, color="r", marker="o",label="Set A")
plot3D(X2[:,1], X2[:,2], alpha2, color="b", marker="^",label="Set B")
plot3D(X3[:,1], X3[:,2], alpha3, color="g", marker="v",label="Set C")
plot3D(X4[:,1], X4[:,2], alpha4, color="m", marker="s",label="Set D")

# scatter3D(x,y,[0,0],color="black",marker="*",s=100.0)

xlabel(L"x_1",fontsize=18)
ylabel(L"x_2",fontsize=18)
zlabel(L"\alpha",fontsize=18)
legend(bbox_to_anchor=(0.5, 1.0), loc="upper center", ncol=4, fancybox="True", shadow="True")

=#