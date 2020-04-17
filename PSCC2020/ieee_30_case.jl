# This is a duplicate of the IEEE matpower 30-bus case

# raw data
# bus data
#	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
bus = [
	1	3	0	0	0	0	1	1	0	135	1	1.05	0.95;
	2	2	21.7	12.7	0	0	1	1	0	135	1	1.1	0.95;
	3	1	2.4	1.2	0	0	1	1	0	135	1	1.05	0.95;
	4	1	7.6	1.6	0	0	1	1	0	135	1	1.05	0.95;
	5	1	0	0	0	0.19	1	1	0	135	1	1.05	0.95;
	6	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	7	1	22.8	10.9	0	0	1	1	0	135	1	1.05	0.95;
	8	1	30	30	0	0	1	1	0	135	1	1.05	0.95;
	9	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	10	1	5.8	2	0	0	3	1	0	135	1	1.05	0.95;
	11	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	12	1	11.2	7.5	0	0	2	1	0	135	1	1.05	0.95;
	13	2	0	0	0	0	2	1	0	135	1	1.1	0.95;
	14	1	6.2	1.6	0	0	2	1	0	135	1	1.05	0.95;
	15	1	8.2	2.5	0	0	2	1	0	135	1	1.05	0.95;
	16	1	3.5	1.8	0	0	2	1	0	135	1	1.05	0.95;
	17	1	9	5.8	0	0	2	1	0	135	1	1.05	0.95;
	18	1	3.2	0.9	0	0	2	1	0	135	1	1.05	0.95;
	19	1	9.5	3.4	0	0	2	1	0	135	1	1.05	0.95;
	20	1	2.2	0.7	0	0	2	1	0	135	1	1.05	0.95;
	21	1	17.5	11.2	0	0	3	1	0	135	1	1.05	0.95;
	22	2	0	0	0	0	3	1	0	135	1	1.1	0.95;
	23	2	3.2	1.6	0	0	2	1	0	135	1	1.1	0.95;
	24	1	8.7	6.7	0	0.04	3	1	0	135	1	1.05	0.95;
	25	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
	26	1	3.5	2.3	0	0	3	1	0	135	1	1.05	0.95;
	27	2	0	0	0	0	3	1	0	135	1	1.1	0.95;
	28	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	29	1	2.4	0.9	0	0	3	1	0	135	1	1.05	0.95;
	30	1	10.6	1.9	0	0	3	1	0	135	1	1.05	0.95;
]

# generator data
#	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
gen = [
	1	23.54	0	150	-20	1	100	1	80	0	0	0	0	0	0	0	0	0	0	0	0;
	2	60.97	0	60	-20	1	100	1	80	0	0	0	0	0	0	0	0	0	0	0	0;
	22	21.59	0	62.5	-15	1	100	1	50	0	0	0	0	0	0	0	0	0	0	0	0;
	27	26.91	0	48.7	-15	1	100	1	55	0	0	0	0	0	0	0	0	0	0	0	0;
	23	19.2	0	40	-10	1	100	1	30	0	0	0	0	0	0	0	0	0	0	0	0;
	13	37	0	44.7	-15	1	100	1	40	0	0	0	0	0	0	0	0	0	0	0	0;
]

# branch data
#	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
branch = [
	1	2	0.02	0.06	0.03	130	130	130	0	0	1	-360	360;
	1	3	0.05	0.19	0.02	130	130	130	0	0	1	-360	360;
	2	4	0.06	0.17	0.02	65	65	65	0	0	1	-360	360;
	3	4	0.01	0.04	0	130	130	130	0	0	1	-360	360;
	2	5	0.05	0.2	0.02	130	130	130	0	0	1	-360	360;
	2	6	0.06	0.18	0.02	65	65	65	0	0	1	-360	360;
	4	6	0.01	0.04	0	90	90	90	0	0	1	-360	360;
	5	7	0.05	0.12	0.01	70	70	70	0	0	1	-360	360;
	6	7	0.03	0.08	0.01	130	130	130	0	0	1	-360	360;
	6	8	0.01	0.04	0	32	32	32	0	0	1	-360	360;
	6	9	0	0.21	0	65	65	65	0	0	1	-360	360;
	6	10	0	0.56	0	32	32	32	0	0	1	-360	360;
	9	11	0	0.21	0	65	65	65	0	0	1	-360	360;
	9	10	0	0.11	0	65	65	65	0	0	1	-360	360;
	4	12	0	0.26	0	65	65	65	0	0	1	-360	360;
	12	13	0	0.14	0	65	65	65	0	0	1	-360	360;
	12	14	0.12	0.26	0	32	32	32	0	0	1	-360	360;
	12	15	0.07	0.13	0	32	32	32	0	0	1	-360	360;
	12	16	0.09	0.2	0	32	32	32	0	0	1	-360	360;
	14	15	0.22	0.2	0	16	16	16	0	0	1	-360	360;
	16	17	0.08	0.19	0	16	16	16	0	0	1	-360	360;
	15	18	0.11	0.22	0	16	16	16	0	0	1	-360	360;
	18	19	0.06	0.13	0	16	16	16	0	0	1	-360	360;
	19	20	0.03	0.07	0	32	32	32	0	0	1	-360	360;
	10	20	0.09	0.21	0	32	32	32	0	0	1	-360	360;
	10	17	0.03	0.08	0	32	32	32	0	0	1	-360	360;
	10	21	0.03	0.07	0	32	32	32	0	0	1	-360	360;
	10	22	0.07	0.15	0	32	32	32	0	0	1	-360	360;
	21	22	0.01	0.02	0	32	32	32	0	0	1	-360	360;
	15	23	0.1	0.2	0	16	16	16	0	0	1	-360	360;
	22	24	0.12	0.18	0	16	16	16	0	0	1	-360	360;
	23	24	0.13	0.27	0	16	16	16	0	0	1	-360	360;
	24	25	0.19	0.33	0	16	16	16	0	0	1	-360	360;
	25	26	0.25	0.38	0	16	16	16	0	0	1	-360	360;
	25	27	0.11	0.21	0	16	16	16	0	0	1	-360	360;
	28	27	0	0.4	0	65	65	65	0	0	1	-360	360;
	27	29	0.22	0.42	0	16	16	16	0	0	1	-360	360;
	27	30	0.32	0.6	0	16	16	16	0	0	1	-360	360;
	29	30	0.24	0.45	0	16	16	16	0	0	1	-360	360;
	8	28	0.06	0.2	0.02	32	32	32	0	0	1	-360	360;
	6	28	0.02	0.06	0.01	32	32	32	0	0	1	-360	360;
]

## processed data for future use

# N = number of nodes
# A = adjacency matrix of physical transmission lines
# F = flow capacity of physical transmission lines
# V = adjacency matrix of spatial virtual links
# D = load capacity array
# S = power supply capacity array
# α_d = load bidding prices
# α_s = supply bidding prices

N = size(bus)[1]
A = zeros(Int8, N, N)
F = zeros(Float64, N, N)
V = zeros(Int8, N, N)
D = zeros(Float64, N)
S = zeros(Float64, N)
α_d = zeros(Float64, N)
α_s = zeros(Float64, N)
Del = zeros(Float64, N, N)
B = zeros(N,N)

for i in 1:N
	D[i] = sqrt(bus[i,3]^2 + bus[i,4]^2)
	α_d[i] = 20.0 + 40.0/N*i
end


num_supply = size(gen)[1]
for i in 1:num_supply
	node_supply = Int(gen[i,1])
	S[node_supply] = gen[i,9]
	α_s[node_supply] = 10 + 20/num_supply*i
end

for i in 1:size(branch)[1]
	node1 = Int(branch[i,1])
	node2 = Int(branch[i,2])
	temp_capacity = sqrt(branch[i,3]^2 + branch[i,4]^2)*90
	A[node1,node2] = 1
	A[node2,node1] = 1
	F[node1,node2] = temp_capacity
	F[node2,node1] = temp_capacity
	susceptance = branch[i,4]/(branch[i,3]^2+branch[i,4]^2)
	B[node1,node2] = susceptance
	B[node2,node1] = susceptance
end

α_f = 5 * copy(F .> 0)

## THE FOLLOWINGS ARE PARAMETERS FOR SPACE-TIME NETWORK ONLY
using Random
T = 3
rng = MersenneTwister(0)

#=
demand_fluc = 5 * randn(rng, Float64, (N,T))
supply_fluc = 5 * randn(rng, Float64, (N,T))
α_d_fluc = 5 * randn(rng, Float64, (N,T))
α_s_fluc = 5 * randn(rng, Float64, (N,T))
D_st = D .+ demand_fluc
S_st = S .+ supply_fluc
D_st = D_st .* (D_st .> 0)
S_st = S_st .* (S_st .> 0)
α_s_st = α_s .+ α_s_fluc
α_d_st = α_d .+ α_d_fluc
α_s_st = α_s_st .* (α_s_st .> 0)
α_d_st = α_d_st .* (α_d_st .> 0)
=#

D_st = zeros(N,T)
S_st = zeros(N,T)
α_d_st = zeros(N,T)
α_s_st = zeros(N,T)

for i in 1:T
	for j in 1:N
		if D[j] > 0
			temp_rndnum = 5 * randn(rng,Float64)
			D_st[j,i] = D[j] + temp_rndnum
			if D_st[j,i] < 0
				D_st[j,i] = 0
			else
				temp_rndnum = 5 * randn(rng,Float64)
				α_d_st[j,i] = α_d[j] + temp_rndnum
			end
		end
		if S[j] > 0
			temp_rndnum = 5 * randn(rng,Float64)
			S_st[j,i] = S[j] + temp_rndnum
			if S_st[j,i] < 0
				S_st[j,i] = 0
			else
				temp_rndnum = 5 * randn(rng,Float64)
				α_s_st[j,i] = α_s[j] + temp_rndnum
			end
		end
	end
end

D_st = transpose(D_st)
S_st = transpose(S_st)
α_s_st = transpose(α_s_st)
α_d_st = transpose(α_d_st)

# only assign spatial shift capacity if each of the nodes have a positive demand level
for i in 1:N, j in (i+1):N
	if D[i] > 0 && D[j] > 0
		temp = rand(rng)
		if temp < 0.7
			Del[i,j] = 10
			Del[j,i] = 10
			#println("Positive shift capacity between $(i) and $(j)")
		else
			#println("No spatial shift between $(i) and $(j)")
		end
	end
end

α_δ = 0.0001 * copy(Del .> 0)
