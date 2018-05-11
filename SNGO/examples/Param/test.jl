using JuMP, Ipopt, DataFrames
push!(LOAD_PATH, ENV["SINGODIR"])
using Plasmo, Ipopt, SCIP
using JuMP
using MAT
mingap = 0.01

# User defined Variables ---------------------------------------
# Number of discretization: decide how fine your grid will be
numOfDiscMono=10 # Number Of Discretization for monospecies experiment
numOfDisc = numOfDiscMono
numOfDiscMulti=240 # Number of Discrteization for multispecies experiment

includeMono=true # Include data of experiment "Mono"
includePair=false # Include data of experiment "Pair"
includePW22=false # Include data of experiment "PW22"
includeTM1=false # Include data of experiment "TM1"
includeL1O6=false # Include data of experiment "L1O6"

equalWeight=true # Equal weight vs Original weight

# Setting bound, +Inf or -Inf indicates that there will be no bound
# Note that setting problem unbounded may make problem infeasible
upperBoundR1=2
lowerBoundR1=0.001

upperBoundR2=-0.01
lowerBoundR2=-10

upperBoundA=10
lowerBoundA=-10

# Species order. This wil decide the order of parameters
speciesOrder=[
"B" "H" # 1
"C" "A" # 2
"B" "U" # 3
"P" "C" # 4
"B" "O" # 5
"B" "V" # 6
"B" "T" # 7
"E" "L" # 8
"F" "P" # 9
"C" "H" # 10
"D" "P" # 11
"E" "R" # 12
]
numOfSpecies=size(speciesOrder,1)
if equalWeight==true
  include("data_processing_equal_weight.jl")
else
  include("data_processing.jl")
end





function createMono(j::Int)
    m = Model()
    @variable(m, lowerBoundR1<=r[j]<=upperBoundR1)
    @variable(m, lbMat[j,j]<=a[j,j]<=ubMat[j,j])

    lt0 = length(time0)
    @variable(m, 0<=y[i=0:numOfDisc*(length(time0)-1)]<=1)
    @variable(m, 0<=ysquare[j,i=1:numOfDisc*(length(time0)-1)]<=1)
    @constraint(m, square[i=1:numOfDisc*(length(time0)-1)], ysquare[i] == y[i]^2)

    # Dynamic Constraint
    @constraint(m, ODEMono[j,i=1:numOfDisc*(length(time0)-1)], y[j,i]==y[j,i-1]+(r[j]*y[j,i]+a[j,j]*ysquare[j,i])/numOfDisc*0.5)

    # Initial Conditions
    @constraint(m, y[j,0]==data[j].IC[1])

    # Objective Function
    @variable(m, objective)
    @constraint(m, objective>=10000*(sum{(y[j,numOfDisc*(i-1)]-data[j].abundance[1][i])^2, i=1:lt0}))
    @objective(m, Min, objective)
    return m
end

#=
m = createMono(1)
m = copyNLModel(m)
mcopy = copyNLModel(m)
mcopy.solver = IpoptSolver()
mcopy_status = solve(mcopy)
m.colVal = copy(mcopy.colVal)
#m.solver = SCIPSolver("limits/gap", mingap, "limits/absgap", mingap, "limits/time", 86400.0, "setobjlimit", mcopy.objVal)
m.solver = BaronSolver()
push!(m.solver.options, (:EpsA, mingap), (:EpsR, mingap), (:MaxTime, 86400.0)) #,(:PrLevel, 0))
solve(m)
=#



lt0 = length(time0)
function createMonoPiece(j::Int, t, piece, npiece)
    m = Model()
    @variable(m, lowerBoundR1<=r[j]<=upperBoundR1)
    @variable(m, lbMat[j,j]<=a[j,j]<=ubMat[j,j])
    #@variable(m, lowerBoundA<=a[j,j]<=upperBoundA)

    if piece == 1
       @variable(m, 0<=yc[piece]<=1)
    elseif piece == npiece
       @variable(m, 0<=yc[piece-1]<=1)
    else
       @variable(m, 0<=yc[(piece-1):piece]<=1)
    end

    @variable(m, 0<=y[i=0:numOfDisc*(length(t)-1)]<=1)
    @variable(m, 0<=ysquare[i=1:numOfDisc*(length(t)-1)]<=1)
    @constraint(m, square[i=1:numOfDisc*(length(t)-1)], ysquare[i] == y[i]^2)

    # Dynamic Constraint
    @constraint(m, ODEMono[i=1:numOfDisc*(length(t)-1)], y[i]==y[i-1]+(r[j]*y[i]+a[j,j]*ysquare[i])*0.5/numOfDisc)

    # Boundary Conditions
   if piece == 1
       @constraint(m, y[0]==data[j].IC[1])
       @constraint(m, y[numOfDisc*(length(t)-1)]==yc[piece])
    elseif piece == npiece
       @constraint(m, y[0]==yc[piece-1])
    else
       @constraint(m, y[0]==yc[piece-1])
       @constraint(m, y[numOfDisc*(length(t)-1)]==yc[piece])
    end

    # Objective Function
    @variable(m, objective)
    ind_start = round(Int, t[1]/0.5)

    if ind_start == 0
        @constraint(m, objective>=10000*(sum{(y[numOfDisc*(i-1)]-data[j].abundance[1][i+ind_start])^2, i=1:(length(t))}))
    else
	@constraint(m, objective>=10000*(sum{(y[numOfDisc*(i-1)]-data[j].abundance[1][i+ind_start])^2, i=2:(length(t))}))
    end
    @objective(m, Min, objective)
    return m
end


j = 1
println("j:   ", j)
lt0 = length(time0)
nstepPpiece = 1   
npiece = Int(ceil((lt0-1)/nstepPpiece))
println("npiece:   ",npiece)
m = NetModel()
@variable(m, lowerBoundR1<=r[j]<=upperBoundR1)
@variable(m, lbMat[j,j]<=a[j,j]<=ubMat[j,j])
@variable(m, 0<=yc[1:(npiece-1)]<=1)
for i = 1:npiece
           if i < npiece
              t = time0[(1+nstepPpiece*(i-1)):(nstepPpiece*i+1)]
           else
              t = time0[(1+nstepPpiece*(i-1)):lt0]
           end
           println(t)
           node = createMonoPiece(j, t, i, npiece)
           # add first-stage variables and constraints
           @addNode(m, node, "s$i")
           @constraint(m, getvariable(node, :r)[j]==r[j])
           @constraint(m, getvariable(node, :a)[j,j]==a[j,j])
           if i == 1
               @constraint(m, getvariable(node, :yc)[i]==yc[i])
           elseif i == npiece
               @constraint(m, getvariable(node, :yc)[i-1]==yc[i-1])
           else
               @constraint(m, getvariable(node, :yc)[i-1]==yc[i-1])
               @constraint(m, getvariable(node, :yc)[i]==yc[i])
           end
end


m.solver = SCIPSolver("limits/gap", mingap, "limits/absgap", mingap, "limits/time", 43200.0)
m = copyNLModel(m)
solve(m)
