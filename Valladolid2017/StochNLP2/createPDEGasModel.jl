function createPDEGasModel(s)
     m = Model(solver=IpoptSolver())
         @variable(m, nodeDict[j].pmin<=p[j in NODE, TIMEG]<=nodeDict[j].pmax, start= 50)         # node pressure - [bar]
         @variable(m, 0<=dp[j = LINK,  TIMEG; linkDict[j].ltype == "a"]<=100, start= 10)          # compressor boost - [bar]
         @variable(m, 1<=fin[LINK, TIMEG]<=500, start= 100)                                       # flow in pipe - [scmx10-4/hr]
         @variable(m, 1<=fout[LINK, TIMEG]<=500, start= 100)                                      # flow out pipe - [scmx10-4/hr]
         @variable(m, 0.01<=sG[j in SUP, TIMEG]<=supDict[j].max, start = 10)                      # supply flow - [scmx10-4/hr]
         @variable(m, dem[DEM, TIMEG],    start=100)                                              # demand flow - [scmx10-4/hr]
         @variable(m, 0<=pow[j = LINK, TIMEG; linkDict[j].ltype == "a"]<=3000, start= 1000)       # compressor power - [kW]
         @variable(m, slack[LINK, TIMEG, DIS]>=0, start= 10)                                      # auxiliary variable
         # define spatio-temporal variables
         @variable(m, 10<=px[LINK, TIMEG, DIS]<=100, start= 50)                             # link pressure profile - [bar]
         @variable(m, 1<=fx[LINK, TIMEG, DIS]<=100, start= 100)                             # link flow profile - [scmx10-4/hr]

    # compressor equations
         @NLconstraint(m, powereq[j = LINK, t = TIMEG; linkDict[j].ltype == "a"], pow[j,t] == c4*fin[j,t]*(((p[linkDict[j].startloc,t] + dp[j,t])/(p[linkDict[j].startloc,t]))^om-1))

    # node balance [mass]
         @constraint(m,nodemeq[i in NODE, t = TIMEG], sum(       fout[j,t] for j in LINK if linkDict[j].endloc==i)
                                                                 - sum(       fin[j,t] for  j in LINK if linkDict[j].startloc==i)
                                                                 + sum(        sG[j,t] for  j in SUP if supDict[j].loc == i )
                                                                 - sum(       dem[j,t] for  j in DEM if demDict[j].loc == i )
                                                                 ==0)

    # flow equations for passive and active links
         @constraint(m, flow[j = LINK, t = TIMEGm, k = 1:(Nx-1)], (px[j,t+1,k]-px[j,t,k])/dtG + linkDict[j].c1*(fx[j,t+1,k+1]-fx[j,t+1,k])/(linkDict[j].dx)==0)

     # boundary conditions flow
         @constraint(m, flow_start[j = LINK, t = TIMEG], fx[j,t,1]==fin[j,t])
         @constraint(m, flow_end[j = LINK, t = TIMEG], fx[j,t,Nx]==fout[j,t])

     # pressure equations for passive and active links
     @constraint(m, press[j = LINK, t = TIMEGm,k = 1:(Nx-1)], (fx[j,t+1,k]-fx[j,t,k])/dtG ==
        - linkDict[j].c2*(px[j,t+1,k+1]-px[j,t+1,k])/linkDict[j].dx - slack[j,t+1,k])
     @NLconstraint(m, slackeq[j = LINK, t = TIMEG, k = 1:Nx],  slack[j,t,k]*px[j,t,k] - linkDict[j].c3*fx[j,t,k]*fx[j,t,k] == 0);

     # boundary conditions pressure, passive links
         @constraint(m, presspas_start[j = LINK, t = TIMEG; linkDict[j].ltype == "p"], px[j,t,1] ==  p[linkDict[j].startloc,t])
         @constraint(m,   presspas_end[j = LINK, t = TIMEG; linkDict[j].ltype == "p"], px[j,t,Nx] == p[linkDict[j].endloc,t])

     # boundary conditions, active links
      @constraint(m, pressact_start[j = LINK, t = TIMEG; linkDict[j].ltype == "a"], px[j,t,1] ==  p[linkDict[j].startloc,t] + dp[j,t])
      @constraint(m,   pressact_end[j = LINK, t = TIMEG; linkDict[j].ltype == "a"], px[j,t,Nx] == p[linkDict[j].endloc,t])

     # fix pressure at supply nodes
         @constraint(m, suppres[j in SUP, t = TIMEG], p[supDict[j].loc,t] == nodeDict[supDict[j].loc].pmin)

     # discharge pressure for compressors
     @constraint(m, dispress[j in LINK,t in TIMEG; linkDict[j].ltype=="a"],  p[linkDict[j].startloc,t]+dp[j,t] <= nodeDict[linkDict[j].startloc].pmax)

     # line pack constraints
         @constraint(m, line_packT,  sum(sum(fx[j,Nt,k] for k in DIS)*linkDict[j].dx for j in LINK) >= sum(sum(fx[j,1,k] for k in DIS)*linkDict[j].dx for j = LINK))

     # ss constraints
         @constraint(m, flow_ss[j = LINK, t =0, k = 1:(Nx-1)], (fx[j,t+1,k+1]-fx[j,t+1,k])==0)
         @constraint(m, pres_ss[j = LINK, t =0, k = 1:(Nx-1)],  - linkDict[j].c2*(px[j, t+1,k+1]-px[j,t+1,k])/linkDict[j].dx - slack[j,t+1,k] == 0)

         @objective(m, Min, 1e-3*(1.0/S)*sum(ce*pow[j,t]*(dtG/3600) for j = LINK , t = TIMEG if linkDict[j].ltype == "a")
                           +1e-3*(1.0/S)*sum(cd*(dem[j,t]-demDict[j].stochd[s,t])^2 for j in DEM, t in TIMEG))

    return m
end
