#******************************************************************************#
#********** FILE OPTIONS, FOLDER PATHWAYS, AND HARD-CODED PARAMETERS **********#
#******************************************************************************#
################################################################################
### NOTES ###
#=
> Records solution variables to individual .csv files

> This file is intended to be run following a TimedMarketsModel.jl run;
  so no packages are imported here.

> Includes:
s_out[n,p,z]
d_out[n,p,z]
v_out[n,p,z]
sl_out[i]
dl_out[j]
vl_out[k]
f_out[n,m,p,z]
x_out[t,n,q,z]
g_out[t,n,p,z]
cp_out[n,p,z]
cp_storage[n,p,z]
cp_transport[n,m,p,z]
cp_tech[t,n,z]
PhiSupply[i]
PhiDemand[j]
PhiStorage[k]
PhiTransport[n,m,p]
PhiTech[t,n]
=#
################################################################################
### DIRECTORY ###
# Output folder name
RecordFolder = "/records"

################################################################################
### WRITE VARIABLES TO FILES ###

# sl[i]
filename = open(DataFolder*RecordFolder*"/supply.csv","w")
print(filename, "# 1.Supply ID| 2.Supplier ID| 3.Node| 4.Product| 5.Period| 6.Variable Value")
for i in S
    print(filename, "\n"*i*"|"*sup_sup[i]*"|"*sup_node[i]*"|"*sup_prod[i]*"|"*sup_time[i]*"|"*string(sl_out[i]))
end
close(filename)

# dl[j]
filename = open(DataFolder*RecordFolder*"/demand.csv","w")
print(filename, "# 1.Demand ID| 2.Node| 3.Product| 4.Period| 5.Variable Value")
for j in D
    print(filename, "\n"*j*"|"*dem_node[j]*"|"*dem_prod[j]*"|"*dem_time[j]*"|"*string(dl_out[j]))
end
close(filename)

# vl[k,z]
filename = open(DataFolder*RecordFolder*"/storage.csv","w")
print(filename, "# 1.Storage ID| 2.Node| 3.Product| 4.Period| 5.Variable Value")
for k in K
    for z in Z
        print(filename, "\n"*k*"|"*sto_node[k]*"|"*sto_prod[k]*"|"*z*"|"*string(vl_out[k,z]))
    end
end
close(filename)

# s[n,p,z]
filename = open(DataFolder*RecordFolder*"/nodal_supply.csv","w")
print(filename, "# 1.Node| 2.Product| 3.Period| 4.Variable Value")
for n in N
    for p in P
        for z in Z
            print(filename, "\n"*n*"|"*p*"|"*z*"|"*string(s_out[n,p,z]))
        end
    end
end
close(filename)

# d[n,p,z]
filename = open(DataFolder*RecordFolder*"/nodal_demand.csv","w")
print(filename, "# 1.Node| 2.Product| 3.Period| 4.Variable Value")
for n in N
    for p in P
        for z in Z
            print(filename, "\n"*n*"|"*p*"|"*z*"|"*string(d_out[n,p,z]))
        end
    end
end
close(filename)

# v[n,p,z]
filename = open(DataFolder*RecordFolder*"/nodal_storage.csv","w")
print(filename, "# 1.Node| 2.Product| 3.Period| 4.Variable Value")
for n in N
    for p in P
        for z in Z
            print(filename, "\n"*n*"|"*p*"|"*z*"|"*string(v_out[n,p,z]))
        end
    end
end
close(filename)

# f[n,m,p,z]
filename = open(DataFolder*RecordFolder*"/transport.csv","w")
print(filename, "# 1.Arc| 2.Product| 3.Period| 4.Variable Value")
for a in A
    for p in P
        for z in Z
            print(filename, "\n"*a*"|"*p*"|"*z*"|"*string(f_out[a,p,z]))
        end
    end
end
close(filename)

# x[t,n,q,z]
filename = open(DataFolder*RecordFolder*"/consumption.csv","w")
print(filename, "# 1.Technology| 2.Node| 3.Product| 4.Period| 5.Variable Value")
for t in T
    for n in N
        for q in Q
            for z in Z
                print(filename, "\n"*t*"|"*n*"|"*q*"|"*z*"|"*string(x_out[t,n,q,z]))
            end
        end
    end
end
close(filename)

# g[t,n,p,z]
filename = open(DataFolder*RecordFolder*"/generation.csv","w")
print(filename, "# 1.Technology| 2.Node| 3.Product| 4.Period| 5.Variable Value")
for t in T
    for n in N
        for p in P
            for z in Z
                print(filename, "\n"*t*"|"*n*"|"*p*"|"*z*"|"*string(g_out[t,n,p,z]))
            end
        end
    end
end
close(filename)

# cp[n,p,z]
filename = open(DataFolder*RecordFolder*"/nodal_price.csv","w")
print(filename, "# 1.Node| 2.Product| 3.Period| 4.Variable Value")
for n in N
    for p in P
        for z in Z
            print(filename, "\n"*n*"|"*p*"|"*z*"|"*string(cp_out[n,p,z]))
        end
    end
end
close(filename)

# cp_storage[n,p,z]
filename = open(DataFolder*RecordFolder*"/storage_price.csv","w")
print(filename, "# 1.Node| 2.Product| 3.Period| 4.Variable Value")
for n in N
    for p in P
        for z in Zz
            print(filename, "\n"*n*"|"*p*"|"*z*"|"*string(cp_storage[n,p,z]))
        end
    end
end
close(filename)

# cp_transport[n,m,p,z]
filename = open(DataFolder*RecordFolder*"/transport_price.csv","w")
print(filename, "# 1.Node 1| 2.Node 2| 3.Product| 4.Period| 5.Variable Value")
for a in A
    for p in P
        for z in Z
            print(filename, "\n"*a*"|"*p*"|"*z*"|"*string(cp_transport[a,p,z]))
        end
    end
end
close(filename)

# cp_tech[t,n,z]
filename = open(DataFolder*RecordFolder*"/technology_price.csv","w")
print(filename, "# 1.Technology| 2.Node| 3.Product| 4.Period| 5.Variable Value")
for t in T
    for n in N
        for z in Z
            print(filename, "\n"*t*"|"*n*"|"*z*"|"*string(cp_tech[t,n,z]))
        end
    end
end
close(filename)

# PhiSupply[i]
filename = open(DataFolder*RecordFolder*"/supply_profit.csv","w")
print(filename, "# 1.Supply ID| 2.Supplier ID| 3.Node| 4.Product| 5.Period| 6.Variable Value")
for i in S
    print(filename, "\n"*i*"|"*sup_sup[i]*"|"*sup_node[i]*"|"*sup_prod[i]*"|"*sup_time[i]*"|"*string(PhiSupply[i]))
end
close(filename)

# PhiDemand[i]
filename = open(DataFolder*RecordFolder*"/demand_profit.csv","w")
print(filename, "# 1.Demand ID| 2.Node| 3.Product| 4.Period| 5.Variable Value")
for j in D
    print(filename, "\n"*j*"|"*dem_node[j]*"|"*dem_prod[j]*"|"*dem_time[j]*"|"*string(PhiDemand[j]))
end
close(filename)

# PhiStorage[i]
filename = open(DataFolder*RecordFolder*"/storage_profit.csv","w")
print(filename, "# 1.Storage ID| 2.Node| 3.Product| 4.Variable Value")
for k in K
    print(filename, "\n"*k*"|"*sto_node[k]*"|"*sto_prod[k]*"|"*string(PhiStorage[k]))
end
close(filename)

# PhiTransport[n,m,p]
filename = open(DataFolder*RecordFolder*"/transport_profit.csv","w")
print(filename, "# 1.Arc| 2.Product| 3.Variable Value")
for a in A
    for p in P
        print(filename, "\n"*a*"|"*p*"|"*string(PhiTransport[a,p]))
    end
end
close(filename)

# PhiTech[t,n]
filename = open(DataFolder*RecordFolder*"/technology_profit.csv","w")
print(filename, "# 1.Technology| 2.Node| 3.Variable Value")
for t in T
    for n in N
        print(filename, "\n"*t*"|"*n*"|"*string(PhiTech[t,n]))
    end
end
close(filename)

################################################################################
### END OF CODE ###
println(PrintSpacer,"\nRecords Saved!\n",PrintSpacer)