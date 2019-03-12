using DifferentialEquations
using PyPlot
using Interpolations
using Distributions
using JuMP
using Ipopt
using JLD
include("setup.jl")


C = 10
FR_band0 = P_nominal*C                          # Initial FR Band Per area
MPC(u0, "Heuristic")



#=
function f_CC(out,du,u,p,t)
    ratetest = p[2]
    f_common(out,du,u,p,t)
    if Sei
        it = u[Ncp+Ncn+5]
    else
        iint = u[Ncp+Ncn+1]
        it = iint
    end
    out[end] = it - TC*ratetest
end

function CC(u0, ratetest,endvolt1,endvolt2)
    function stop_cond(u,t,integrator)
        phi_p = u[Ncp+Ncn+2]
        phi_n = u[Ncp+Ncn+3]

            cf = u[Ncp+Ncn+Nsei+8]
            capacity_remain = 1 - cf/Qmax
            csn_avg = u[Ncp+1]
            soc = csn_avg/csnmax

	if t <= 100
	    return false
	end
	if phi_p - phi_n >= (endvolt1 - 0.005) || phi_p - phi_n <= (endvolt2+0.005)
	    return true
        end    
	if soc >= capacity_remain - 0.005
	    return true
	end    
	return false
    end
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(stop_cond,affect!)

    du0 = zeros(Ncp+Ncn+4+Nsei+Ncum)
    tspan = (0.0, 10000.0)
    p = [ratetest]
    prob = DAEProblem(f_CC,du0,u0,tspan,p,differential_vars=differential_vars)
    sol = DifferentialEquations.solve(prob, IDA(), callback=cb, reltol=1e-8, abstol=1e-8)
    return sol
end

function f_CV(out,du,u,p,t)
    endvolt1 = p[2]
    phi_p = u[Ncp+Ncn+2]
    phi_n = u[Ncp+Ncn+3]
    f_common(out,du,u,p,t)
    out[end] = phi_p - phi_n - endvolt1
end

function CV(u0, endvolt1)
    function stop_cond(u,t,integrator)
        if Sei
            it = u[Ncp+Ncn+5]
        else
            iint = u[Ncp+Ncn+1]
            it = iint
        end
        return (it - 0.1)
    end
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(stop_cond,affect!)

    du0 = zeros(Ncp+Ncn+4+Nsei+Ncum)
    tspan = (0.0, 10000.0)
    p = [endvolt1]
    prob = DAEProblem(f_CV,du0,u0,tspan,p,differential_vars=differential_vars)
    sol = DifferentialEquations.solve(prob, IDA(), callback=cb)
    return sol
end


function CCCV(crate)
    ICs = copy(u0)
    ICs[Ncp+Ncn+10] = TC*crate
    sol = CC(ICs, 1, crate, 4.2, 2.5)
    t = sol.t
    index = 1
    Q = sol[Ncp+Ncn+6+Nsei, index:end]
    pot = sol[Ncp+Ncn+4, index:end]
    ICs = sol[end]
    sol_CV = CV(ICs, 1, 4.2)
    t_CV = sol_CV.t
    Q_CV = sol_CV[Ncp+Ncn+6+Nsei, index:end]
    pot_CV = sol_CV[Ncp+Ncn+4, index:end]
    return Q, pot, Q_CV, pot_CV
end


Q05C,pot05C, Q05V,pot05V = CCCV(0.5)
Q10C,pot10C, Q10V,pot10V = CCCV(1.0)
Q15C,pot15C, Q15V,pot15V = CCCV(1.5)
Q20C,pot20C, Q20V,pot20V = CCCV(2.0)

fig = figure("CCCV",figsize=(10,10))
plot(Q05C,pot05C, color="green", label="rate = 0.5, CC")
plot(Q10C,pot10C, color="green", label="rate = 1.0, CC")
plot(Q15C,pot15C, color="green", label="rate = 1.5, CC")
plot(Q20C,pot20C, color="green", label="rate = 2.0, CC")

plot(Q05V,pot05V, color="blue", label="rate = 0.5, CV")
plot(Q10V,pot10V, color="blue", label="rate = 1.0, CV")
plot(Q15V,pot15V, color="blue", label="rate = 1.5, CV")
plot(Q20V,pot20V, color="blue", label="rate = 2.0, CV")

fsl = 18
fsa = 18
xlabel("Q",fontsize=fsl)
ylabel("pot",fontsize=fsl)
xlim(0, 35)
ax=gca()
setp(ax[:get_yticklabels](),fontsize=fsa)
setp(ax[:get_xticklabels](),fontsize=fsa)
legend()
ax[:legend](fontsize=fsa)
grid("on")
PyPlot.tight_layout()
savefig("CCCV.pdf")
=#


