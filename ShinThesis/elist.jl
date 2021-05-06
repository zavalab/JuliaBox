
using DelimitedFiles

name = "case30000_goc"
writedlm("elist/$name.csv", [["Source" "Target"]; vcat([[e.src e.dst] for e in edges(g)]...)],',')
