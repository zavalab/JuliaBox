# Created by Alex Dowling at the University of Wisconsin-Madison in 2016

using JuMP, StructJuMP, PyPlot

function convertToArray2(A,n,keys=[];rnd=false)
	# Converts input A, which is a JuMPDict into an array
	
	if(length(n) != 2)
		warn("Expecting length(A) == 2")
	end
	
	AA = getvalue(A)
	
	B = zeros(n[1], n[2])
	
	if(isempty(keys))
		keys = 1:n[2]
	end
	
	for i in 1:n[1]
		for j in 1:n[2]
			if(rnd)
				B[i,j] = round(AA[i,keys[j]])
			else
				B[i,j] = AA[i,keys[j]]
			end
		end
	end
	
	return B
end

type ProblemData

	days::Array{Int64,1}
	nBlocks::Int64
	daysInBlocks::Array{Array{Int64,1},1}

	maxIn::Array{Float64,2}
	price::Array{Float64,2}

end

function ProblemData(days::Array{Int64,1},nBlocks::Int64, A::Array{Float64,2} = [0.95 0.05; 0.6 0.4])

	nD = length(days)
	
	### Allocate days evenly to blocks
	daysPerBlock = div(nD,nBlocks)
	nextDay = 1
	
	daysInBlocks = Array(Array{Int64,1}, 0)
	
	for i = 1:nBlocks
	
		if(i < nBlocks)
			j = nextDay:(nextDay + daysPerBlock - 1)
		else
			j = nextDay:nD
		end
		
		
		push!(daysInBlocks, days[j])
		
		nextDay = j[end] + 1
	end

	for i = 1:nBlocks
		println("Block ",i," ","contains ",length(daysInBlocks[i])," days")
		println(daysInBlocks[i])
		println(" ") 
	end


	nD = length(days)

	# Seed random number generator
	srand(9000)

	# Generate input data
	maxIn = rand(nH, nD)
	price = rand(nH, nD)

	maxIn_ = 0.5
	price_ = 0.5

	# Make the signals auto-correlated
	for i = 1:nD
		for j = 1:nH
			maxIn[j,i] = A[1,1]*maxIn_ + A[1,2]*maxIn[j,i]
			price[j,i] = A[2,1]*price_ + A[2,2]*price[j,i]
		
			maxIn_ = maxIn[j,i]
			price_ = price[j,i]
		end
	end

	return ProblemData(days, nBlocks, daysInBlocks, maxIn, price)

end

type Solution

	E::Array{Float64,2}
	Pin::Array{Float64,2}
	Pout::Array{Float64,2}

end

function Solution(m::JuMP.Model, pd::ProblemData)

	r = sum(getvalue(getvariable(m,:revenue)))
	println("Revenue = ",r)
	println("Revenue without storage = ",sum(pd.maxIn .* pd.price)*etaIn*etaOut)

	nD = length(pd.days)

	E_ = zeros(nH,nD)
	Pin_ = zeros(nH, nD)
	Pout_ = zeros(nH, nD)

	for i = 1:pd.nBlocks

		d = pd.daysInBlocks[i]

		dim = [nH, length(d)]

		E_[:,d] = convertToArray2(getvariable(getchildren(m)[i],:E), dim, d)
		Pin_[:,d] = convertToArray2(getvariable(getchildren(m)[i],:Pin), dim, d)
		Pout_[:,d] = convertToArray2(getvariable(getchildren(m)[i],:Pout), dim, d)

	end

	return Solution(E_, Pin_, Pout_)

end



# Declare CONSTANT model parameters
const etaIn = 0.9
const etaOut = 0.7
const maxE = 1.0
const startE = maxE/2
const endE = startE 
const hours = 1:24
const nH = length(hours)

function print_model_size(m::JuMP.Model)
	println(m.numCols, " variables ")
	println(MathProgBase.numconstr(m) ," constraints")
	println(" ")
end

function assemble_stage1!(m::JuMP.Model, pd::ProblemData)

	blocks = 1:pd.nBlocks
	blocks_ = 1:(pd.nBlocks + 1)

	# Compute maximum revenue for each block
	maxRevB = zeros(pd.nBlocks)

	for i = blocks
		maxRevB[i] = sum(sum(pd.price,1)[pd.daysInBlocks[i]])
	end

	# Add linking variables to model
	@variable(m, 0 <= E0[i = blocks_] <= maxE)
	@variable(m, 0 <= revenue[i = blocks] <= maxRevB[i])

	# Add objective
	@objective(m, :Min, -sum{revenue[i], i=blocks})

	# Set initial and final storage energy state bounds
	setupperbound(E0[1], startE)
	setlowerbound(E0[pd.nBlocks + 1], endE)

	println("Stage 1:")
	print_model_size(m)

	return nothing

end

function assemble_stage2_block!(m::JuMP.Model, pd::ProblemData, i::Int64, for_pips::Bool, flat::Bool)

	if(flat)
		bl = m
	else
		bl = StructuredModel(parent=m, id=i)
	end
	
	bdays = pd.daysInBlocks[i]

	# Declare state variable (storage energy level)
	@variable(bl, 0 <= E[t = hours, d = bdays] <= maxE)

	# Declare energy flow variables
	@variable(bl, 0 <= Pin[t = hours, d = bdays] <= pd.maxIn[t, d])
	@variable(bl, 0 <= Pout[t = hours, d = bdays] <= 1)

	# Load stage 1 variables
	E0 = getvariable(m, :E0)
	revenue = getvariable(m, :revenue)

	# Declare energy balance
	d_ = bdays[1]
	t_ = 0

	# Iterate over timesteps
	for d = bdays
		for t = hours
	
			# For first timestep, use linking variable
			if(t_ == 0)
				E_ = E0[i]
							
			# For other timesteps, use energy state from previous timestep
			else
				E_ = E[t_, d_]
			end
	
			if(for_pips && t_ == 0)
				@NLconstraint(bl, E[t,d] == E_ + Pin[t,d]*etaIn - Pout[t,d]/etaOut)
			else
				@constraint(bl, E[t,d] == E_ + Pin[t,d]*etaIn - Pout[t,d]/etaOut)
			end
	
			d_ = d
			t_ = t
		end
	end

	# Link energy state
	if(for_pips)
		@NLconstraint(bl, E0[i + 1] == E[hours[end],bdays[end]])
	else
		@constraint(bl, E0[i + 1] == E[hours[end],bdays[end]])
	end

	# Compute revenue for block
	if(for_pips)
		@NLconstraint(bl, revenue[i] == sum{Pout[t,d]*pd.price[t,d], t=hours, d=bdays})
	else
		@constraint(bl, revenue[i] == sum{Pout[t,d]*pd.price[t,d], t=hours, d=bdays})
	end
	
	if(flat)
		println("Flattened model with Stage 1 and Stage 2, Blocks 1 through ",i,":")
	else
		println("Stage 2, Block ",i,":")
	end
	print_model_size(bl)
	
	return nothing

end

function assemble_stage2!(m::JuMP.Model, pd::ProblemData, for_pips::Bool, flat::Bool)

	# Iterate over each block
	for i in 1:pd.nBlocks
		assemble_stage2_block!(m, pd, i, for_pips, flat)
	end

	return nothing

end

function create_struct_model(pd::ProblemData, for_pips::Bool=false, flat::Bool = false)

	# Create structured model
	m = StructuredModel(num_scenarios = pd.nBlocks)

	assemble_stage1!(m, pd)
	assemble_stage2!(m, pd, for_pips, flat)

	return m

end

function plot_solution(s::Solution, pd::ProblemData)

	nD = length(pd.days)

	t = (1:(nH*nD))/nH

	plot(t, reshape(s.E, nH*nD), color="green",linewidth=2)
	plot(t, reshape(s.Pin, nH*nD), color="red",linewidth=2)
	plot(t, reshape(s.Pout, nH*nD), color="blue",linewidth=2)

	return nothing 
end
