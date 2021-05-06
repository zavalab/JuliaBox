using DelimitedFiles, PowerModels, LightGraphs, Plots, Metis

function add_branch!(g,G,B,i,j,gg,bb)
    add_edge!(g,i,j)
    G[i,j]=G[j,i]=gg
    B[i,j]=B[j,i]=bb
end

function get_network(str,K)
    dat = build_ref(PowerModels.parse_file(str))[:it][:pm][:nw][0]

    g = Graph(length(dat[:bus]))
    G = Dict{Tuple{Int,Int},Float64}()
    B = Dict{Tuple{Int,Int},Float64}()

    cnt = 0
    buslist = Int[]
    buslistinv = Dict{Int,Int}()
    mapper = Dict(1=>1,2=>2,3=>3,9=>4)
    
    for (i,bus) in dat[:bus]
        cnt += 1
        push!(buslist,i)
        buslistinv[i] = cnt
    end
    for (i,br) in dat[:branch]
        gg,bb = calc_branch_y(br)
        add_branch!(g,G,B,buslistinv[br["f_bus"]],buslistinv[br["t_bus"]],gg,bb)
    end
    
    hasgen = Vector{Bool}(undef,length(dat[:bus]))
    for (i,gen) in dat[:gen]
        hasgen[gen["gen_bus"]] =true
    end
    
    return g,G,B,Metis.partition(g,K),hasgen
end




g,G,B,part,hagen = get_network("pglib_opf_case500_tamu.m",4)
pos = readdlm("pos.csv",',')
vals= 0.4*ones(500)


plt = plot_graph(g,pos,vals)
savefig(plt,"hello.pdf")

