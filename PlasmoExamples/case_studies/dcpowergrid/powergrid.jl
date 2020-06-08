using Plasmo
using Ipopt

include("load_data.jl")

gamma = 1e5
oms = [1, 4, 7]
max_iter = 1000
scl = 1

g = args[:g]
N_buses = nv(g)
N_lines = ne(g)
y = args[:y]
del = args[:del]

powergrid = OptiGraph()
@optinode(powergrid,buses[1:N_buses])  #create node buses
@optinode(powergrid,lines[1:N_lines])  #create transmission lines

node_map_in = Dict((bus,ModelNode[]) for bus in buses)
node_map_out = Dict((bus,ModelNode[]) for bus in buses)

line_map = Dict()
edge_map = Dict()
B = Dict()
angle_rate = Dict()
ngens = Dict()
load_map = Dict()

for (i,edge) in enumerate(edges(g))
    line = lines[i]
    v_from = edge.src
    v_to = edge.dst

    edge_map[(v_from,v_to)] = line

    B[line] = y[v_from,v_to]
    angle_rate[line] = del[v_from,v_to]

    bus_from = buses[v_from]
    bus_to = buses[v_to]
    bus_vec = [bus_from,bus_to]
    line_map[line] = bus_vec

    push!(node_map_in[bus_to],line)
    push!(node_map_out[bus_from],line)
end

for i = 1:N_buses
    neighs = neighbors(g,i)
    bus = buses[i]

    ngens[bus] = args[:ng][i]
    load_map[bus] = args[:sd][i]
    bus.ext[:c1] = args[:c1][i]
    bus.ext[:c2] = args[:c2][i]

    bus.ext[:va_lower] = args[:val][i]
    bus.ext[:va_upper] = args[:vau][i]

    bus.ext[:gen_lower] = args[:sl][i]
    bus.ext[:gen_upper] = args[:su][i]
end

for line in lines
    bus_from = line_map[line][1]
    bus_to = line_map[line][2]
    @variable(line,va_i,start = 0)
    @variable(line,va_j ,start = 0)
    @variable(line,flow,start = 0)
    @constraint(line,flow == B[line]*(va_i - va_j))
    delta = angle_rate[line]
    @constraint(line,delta <= va_i - va_j <= -delta)
    @objective(line,Min,1/4*gamma*(va_i - va_j)^2)
end

dual_links = LinkConstraintRef[]
#primal_links = LinkConstraintRef[]
for (i,bus) in enumerate(buses)
    va_lower = bus.ext[:va_lower]
    va_upper = bus.ext[:va_upper]

    gen_lower = bus.ext[:gen_lower]
    gen_upper = bus.ext[:gen_upper]

    @variable(bus,va_lower <= va <= va_upper,start = 0)           #voltage angle
    @variable(bus,P[j=1:ngens[bus]], start = 0)
    for j = 1:length(P)
        set_lower_bound(P[j],gen_lower[j])
        set_upper_bound(P[j],gen_upper[j])  #generators at this bus
    end
    lines_in = node_map_in[bus]
    lines_out = node_map_out[bus]

    @variable(bus,power_in[1:length(lines_in)])
    @variable(bus,power_out[1:length(lines_out)])

    for (j,line) in enumerate(lines_in)
        link = @linkconstraint(powergrid,bus[:power_in][j] == line[:flow])
        push!(dual_links,link)
    end
    for (j,line) in enumerate(lines_out)
        link = @linkconstraint(powergrid,bus[:power_out][j] == line[:flow])
        push!(dual_links,link)
    end

    @constraint(bus,power_balance, sum(bus[:P][j] for j=1:ngens[bus]) + sum(power_in) - sum(power_out) - load_map[bus] == 0)
    @objective(bus,Min,sum(bus.ext[:c1][j]*bus[:P][j] + bus.ext[:c2][j]*bus[:P][j]^2 for j = 1:ngens[bus]))
end
@linkconstraint(powergrid,line_coupling_i[line = lines],line[:va_i] == line_map[line][1][:va])
@linkconstraint(powergrid,line_coupling_j[line = lines],line[:va_j] == line_map[line][2][:va])
primal_links = [line_coupling_i.data;line_coupling_j.data]
