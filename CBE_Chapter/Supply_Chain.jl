# Data sets
N = ["n1", "n2"]; P = ["p1", "p2", "p3"];
S = ["s1(p1)", "s2(p1)"]; D = ["n1(p2)", "n2(p2)", "n2(p3)"];
L = ["n1->n2(p2)", "n1->n2(p3)"]; T = ["p1-p2", "p1-p3"];

# Map nodes to supplies
sloc = Dict("n1" => [S[1], S[2]], "n2" => [])

# Map supplies to products they supply
sprod = ["p1","p1"]; sprod = Dict(zip(S,sprod));

# Define upper bounds on supplies
sub = [1000,500]; sub = Dict(zip(S,sub)); # unit

# Define cost of supplying one unit of a supply's product
sbid = [3,2.5]; sbid = Dict(zip(S, sbid)); # $/unit

# Map nodes to demands
dloc = Dict("n1" => [D[1]], "n2" => [D[2], D[3]])

# Map demands to the products they require
dprod = ["p2","p2","p3"]; dprod = Dict(zip(D,dprod));

# Define upper bounds on demands
dub = [100,200,500]; dub = Dict(zip(D,dub)); # unit

# Define price for one unit of the demand's product
dbid = [300,300,15]; dbid = Dict(zip(D,dbid)) # $/unit

# Map nodes to technologies
ξloc = Dict("n1" => [T[1], T[2]], "n2" => [])

# Define upper bounds on technology conversions
ξub = [1000,1000]; ξub = Dict(zip(T,ξub)); # unit

# Define cost of converting one unit of product
ξbid = [1,0.5]; ξbid = Dict(zip(T,ξbid)); # $/unit

# Product conversion matrices; rows are technologies, columns are products
γ = [-1.0 0.9 0.0;
     -1.0 0.0 0.3]
γ = Dict((T[i], P[j]) => γ[i, j] for j in eachindex(P), i in eachindex(T));

# Map transport data (flocs = supplying node; flocr = receiving node)
flocs = ["n1","n1"]; flocs = Dict(zip(L,flocs));
flocr = ["n2","n2"]; flocr = Dict(zip(L,flocr));

# Define upper bounds on transport
fub = [1000,1000]; fub = Dict(zip(L,fub)); # unit

# Map transport flows to the product they support
fprod = [["p2"], ["p3"]]; fprod = Dict(zip(L,fprod));

# Define cost of transporting a unit of product
fbid = [0.1,0.1]; fbid = Dict(zip(L,fbid)); # $/unit

using Plasmo, Ipopt
graph = OptiGraph()

@optinode(graph, nodes[N])
@optinode(graph, transport[L])

# Define variables/constraints/objectives for node locations
for n in N
    # Define variables based on above mappings
    @variable(nodes[n], s[sloc[n]]>=0)
    @variable(nodes[n], d[dloc[n]]>=0)
    @variable(nodes[n], ξ[ξloc[n]]>=0)

    # Define upper bounds on variables
    @constraint(nodes[n], [i in sloc[n]], s[i] <= sub[i])
    @constraint(nodes[n], [j in dloc[n]], d[j] <= dub[j])
    @constraint(nodes[n], [t in ξloc[n]], ξ[t] <= ξub[t])

    # Define expressions for summing supplies, demands, techs on a node
    @expression(nodes[n], sum_supplies[p in P],
        sum(s[i] for i in sloc[n] if sprod[i] == p)
    )
    @expression(nodes[n], sum_demands[p in P],
        sum(d[j] for j in dloc[n] if dprod[j] == p)
    )
    @expression(nodes[n], sum_techs[p in P], sum(ξ[t] * γ[t, p] for t in ξloc[n]))

    # Define cost expressions for a node
    scost = @expression(nodes[n], sum(sbid[i] * s[i] for i in sloc[n]))
    dcost = @expression(nodes[n], sum(dbid[j] * d[j] for j in dloc[n]))
    ξcost = @expression(nodes[n], sum(ξbid[t] * ξ[t] for t in ξloc[n]))

    # Set objective on each node
    @objective(nodes[n], Max, dcost - scost - ξcost)
end

# Define variables/constraints/objectives for transport nodes
for l in L
    # Define flow variable
    @variable(transport[l], f[fprod[l]] >= 0)

    # Set upper bound on flow variable
    @constraint(transport[l], [j in fprod[l]], f[j] <= fub[l])

    # Define flow in and out of nodes
    @expression(transport[l], flow_in[p in P, n in N],
        sum(f[i] for i in fprod[l] if flocr[l] == n && i == p)
    )
    @expression(transport[l], flow_out[p in P, n in N],
        sum(f[i] for i in fprod[l] if flocs[l] == n && i == p)
    )

    # Set objective (penalizing transport costs)
    @objective(transport[l], Max, - sum(fbid[l] * f[j] for j in fprod[l]))
end

# For each node, do a product balance for each product
for n in N
    for p in P
        # Node balance must be equal to zero
        node_balance = (nodes[n][:sum_supplies][p] - nodes[n][:sum_demands][p] +
            nodes[n][:sum_techs][p] + sum(transport[l][:flow_in][p, n] for l in L) -
            sum(transport[l][:flow_out][p, n] for l in L)
        )
        if node_balance != 0
            @linkconstraint(graph, node_balance == 0)
        end
    end
end

# Set objective
set_to_node_objectives(graph)

set_optimizer(graph, Ipopt.Optimizer)
set_to_node_objectives(graph)
optimize!(graph)

println(objective_value(graph))
println("Node 1 demand solutions = ", value.(graph[:nodes]["n1"][:d]))
println("Node 2 demand solutions = ", value.(graph[:nodes]["n2"][:d]))
println("Technology conversion = ", value.(graph[:nodes]["n1"][:ξ]))
