
# Loading packages:
using JuMP 
using Distributions 
using Ipopt
using PyPlot
set_default_plot_size(20cm, 15cm)

# Generate random data: 
N = 1000
srand(0)
μ = 0; σ = 2; 
d = Normal(μ,σ)
R = rand(d,N);

# Plotting:
plt[:hist](R, bins = 30)
grid("on")
xlabel(L"\xi")
ylabel(L"p(\xi)")
savefig("ex1data.pdf")

# ex1gauss.mod 
function ex1GaussModel(xip)
    m = Model(solver=IpoptSolver(print_level=0))
    @variable(m, x)
    @objective(m, Min, (x-xip)^2 - x*xip)
    solve(m)
    obj = getobjectivevalue(m)
    x = getvalue(x)
    return obj,x
end
        
# solve problems with data points
solex1gauss = zeros(N)
solxex1gauss = zeros(N)

for i=1:N 
    (solex1gauss[i], solxex1gauss[i]) = ex1GaussModel(R[i])
end 

# Plotting: 
sol = solex1gauss;
solx = solxex1gauss;
plt[:hist](sol, bins = 30)
grid("on")
xlabel(L"f(x^*(\xi),\xi)")
ylabel(L"p(f(\cdot))")
savefig("ex1solgauss.pdf")

plt[:hist](solx, bins = 30)
grid("on")
xlabel(L"x^*(\xi)")
ylabel(L"p(x^*(\xi))")
savefig("ex1solxgauss.pdf")

# ex1gausscons.mod

function ex1GaussConsModel(xip)
    m = Model(solver=IpoptSolver(print_level=0))
    @variable(m, -1 <= x <= 1)
    @objective(m, Min, (x-xip)^2 - x*xip)
    
    solve(m)
    obj = getobjectivevalue(m)
    x = getvalue(x)
    
    return obj,x
end
        
# solve problems with data points
solex1gausscons = zeros(N)
solxex1gausscons = zeros(N)

for i=1:N 
    (solex1gausscons[i], solxex1gausscons[i]) = ex1GaussConsModel(R[i])
end 

# Plotting: 
sol = solex1gausscons
solx = solxex1gausscons
plt[:hist](sol, bins = 30)
grid("on")
xlabel(L"f(x^*(\xi),\xi)")
ylabel(L"p(f(x^*(\xi),\xi))")
savefig("ex1solgausscons.pdf")

plt[:hist](solx, bins = 30)
grid("on")
xlabel(L"x^*(\xi)")
ylabel(L"p(x^*(\xi))")
savefig("ex1solxgausscons.pdf") 
