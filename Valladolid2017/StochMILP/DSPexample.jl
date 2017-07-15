
using JuMP, Dsp, MPI

# Comment out this line if you want to run in serial
MPI.Init()

xi = [[7,7] [11,11] [13,13]]

# create JuMP.Model with number of blocks
m = Model(3)

@variable(m, 0 <= x[i=1:2] <= 5, Int)
@objective(m, Min, -1.5 * x[1] - 4 * x[2])
for s = 1:3
    # create a JuMP.Model block linked to m with id s and probability 1/3
    blk = Model(m, s, 1/3)
    @variable(blk, y[j=1:4], Bin)
    @objective(blk, Min, -16 * y[1] + 19 * y[2] + 23 * y[3] + 28 * y[4])
    @constraint(blk, 2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= xi[1,s] - x[1])
    @constraint(blk, 6 * y[1] + y[2] + 3 * y[3] + 2 * y[4] <= xi[2,s] - x[2])
end

solve_types = [:Dual, :Benders, :Extensive]
status = solve(m, solve_type = solve_types[2])

getobjectivevalue(m)

# Comment out this line if you want to run in serial
MPI.Finalize()


