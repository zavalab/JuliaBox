# stochastic model for natural gas network
# victor m. zavala / 2013

linkData = []
nodeData = []
supData = []
demData = []
# read link data
linkinfo=readdlm("DAT_PAPER/linkinfo.tab")
(nr, ~) = size(linkinfo)
nlink = nr - 2
for i=1:nlink
    push!(linkData, LinkData(string(linkinfo[i+2,1]),string(linkinfo[i+2,2]),string(linkinfo[i+2,3]), float(linkinfo[i+2,4]),
    float(linkinfo[i+2,5]), string(linkinfo[i+2,6]),0,0,0,0, 0,0))
end

# read node data
nodeinfo=readdlm("DAT_PAPER/nodeinfo.tab")
(nr, ~) = size(nodeinfo)
nnode = nr -2
for i=1:nnode
    push!(nodeData, NodeData(string(nodeinfo[i+2,1]),float(nodeinfo[i+2,2]),float(nodeinfo[i+2,3])))
end

# read supply data
supinfo=readdlm("DAT_PAPER/supinfo.tab")
(nr, ~) = size(supinfo)
nsup = nr -2
for i=1:nsup
    push!(supData, SupplyData(string(supinfo[i+2,1]),string(supinfo[i+2,2]),float(supinfo[i+2,3]),float(supinfo[i+2,4])))
end

# read demand data
deminfo=readdlm("DAT_PAPER/demandinfo.tab")
(nr, ~) = size(deminfo)
ndem = nr -2
for i=1:ndem
    push!(demData, DemandData(string(deminfo[i+2,1]),string(deminfo[i+2,2]),float(deminfo[i+2,3]),Array{Float64}( S, Nt)))
end


# scaling factors
ffac=(1e+6*rhon)/(24*3600)                     # from scmx10-6/day to kg/s
ffac2=(3600)/(1e+4*rhon)                       # from kg/s to scmx10-4/hr
pfac=1e+5                                      # from bar to Pa
pfac2=1e-5                                     # from Pa to bar
dfac=1e-3                                      # from mm to m
lfac=1e+3                                      # from km to m
dtG=TF/(Nt)

# convert units for input data
for i in 1:length(linkData)
    linkData[i].diam *= dfac                     # from  mm to m
    linkData[i].length *= lfac                   # from  km to m
end

for i in 1:length(supData)
    supData[i].min *= ffac*ffac2	          # from scmx106/day to kg/s and then to scmx10-4/hr
    supData[i].max *= ffac*ffac2   		  # from scmx106/day to kg/s and then to scmx10-4/hr
end

for i in 1: length(demData)
    demData[i].d *= ffac*ffac2			  # from scmx106/day to kg/s and then to scmx10-4/hr
end
for i in 1: length(nodeData)
    nodeData[i].pmin *= pfac*pfac2   # from bar to Pascals and then to bar
    nodeData[i].pmax *= pfac*pfac2   # from bar to Pascals and then to bar
end

# compute friction coefficient
lam = Array{Float64}(length(linkData))
for i in 1:length(linkData)
    linkData[i].lam  = lam[i]
    lam[i] = (2*log10(3.7*linkData[i].diam/(eps*dfac)))^(-2)
end

# compute constants
nu2 = gam*z*R*Tgas/M
A = Array{Float64}(length(linkData))
for i in 1:length(linkData)
    A[i] = (1/4)*pi*linkData[i].diam*linkData[i].diam
    linkData[i].A  = A[i]
    linkData[i].c1 = (pfac2/ffac2)*(nu2/A[i])
    linkData[i].c2 = A[i]*(ffac2/pfac2)
    linkData[i].c3 = A[i]*(pfac2/ffac2)*(8*lam[i]*nu2)/(pi*pi*(linkData[i].diam^5))
    linkData[i].dx = linkData[i].length/(Nx-1)

end
c4 = (1/ffac2)*(Cp*Tgas)


stochdeminfo=readdlm("DAT_PAPER/stoch_scen_paper.dat")
for k in SCENG
    for t in TIMEG
    	demData[1].stochd[k,t] = stochdeminfo[k,t]
    	for j in 2:ndem
	    demData[j].stochd[k, t] = demData[j].d
	    nnz += 1
	end
    end
end

for i=1:nlink
    linkDict[linkData[i].name] = linkData[i]
end

for i=1:nnode
    nodeDict[nodeData[i].name] = nodeData[i]
end

for i=1:nsup
    supDict[supData[i].name] = supData[i]
end

for i=1:ndem
    demDict[demData[i].name] = demData[i]
end
