function get_model_linear(N,q,b)
    n = 2
    p = 1
    nd= 2
    A = [1 1;0 1]
    B = [0;b]
    Q = [q 0;0 0]
    R = 1
    Qf= [1 0;0 1]
    x0= [0;0]

    m = Model(()->MadNLP.Optimizer(print_level=MadNLP.ERROR))
    @variable(m,x[1:N+1,1:n])
    @variable(m,u[1:N,1:p])
    @NLparameter(m,d[i=1:N,j=1:nd]== (j==1 ? sin(2pi/N*i) : 0))
    @constraint(m,[i=1:n],x[1,i]==x0[i])
    @NLconstraint(m,l[i=1:N,j=1:n],x[i+1,j] ==
                  sum(A[j,k]*x[i,k] for k=1:n) + sum(B[j,k]*u[i,k] for k=1:p) + d[i,j])
    @objective(m,Min,.5*sum(x[i,j]*x[i,k]*Q[j,k] for i=1:N for j=1:n for k=1:n)
               + .5*sum(u[i,j]*u[i,k]*R[j,k] for i=1:N for j=1:p for k=1:p)
               + .5*sum(x[N+1,j]*x[N+1,k]*Qf[j,k] for j=1:n for k=1:n))

    m[:n] = n; m[:p] = p; m[:nd] = nd; m[:d] = d
    
    return m
end

function get_model_nonlinear(N,q,b; dt = 1)
    n = 9
    p = 4
    nd= 9
    
    x0= zeros(9)
    Q = [1,1,1,q,q,q,1,1,1]
    Qf= ones(9)
    R = ones(4)
    
    m = Model(()->MadNLP.Optimizer(print_level=MadNLP.ERROR))
    @variable(m,x[1:N+1,1:n],start = 0)
    @variable(m,u[1:N,1:p],start = 0)
    @NLparameter(m,d[i=1:N+1,j=1:nd]== (j==1 ? sin(2pi/N*i) : 0))
    @constraint(m,[i=1:n],x[1,i]==x0[i])
    l = Array{ConstraintRef,2}(undef,N,n)
    for i=1:N
        l[i,1]=@NLconstraint(
            m,x[i+1,1] == x[i,1] + (x[i,2])*dt)
        l[i,2]=@NLconstraint(
            m, x[i+1,2] == x[i,2] + 
            (u[i,1]*cos(x[i,7])*sin(x[i,8])*cos(x[i,9])+u[i,1]*sin(x[i,7])*sin(x[i,9]))*dt)
        l[i,3]=@NLconstraint(
            m, x[i+1,3] == x[i,3] + 
            (x[i,4])*dt)
        l[i,4]=@NLconstraint(
            m, x[i+1,4] == x[i,4] + 
            (u[i,1]*cos(x[i,7])*sin(x[i,8])*sin(x[i,9])-u[i,1]*sin(x[i,7])*cos(x[i,9]))*dt)
        l[i,5]=@NLconstraint(
            m, x[i+1,5] == x[i,5] + 
            (x[i,6])*dt)
        l[i,6]=@NLconstraint(
            m, x[i+1,6] == x[i,6] + 
            (u[i,1]*cos(x[i,7])*cos(x[i,8])-9.8)*dt)
        l[i,7]=@NLconstraint(
            m, x[i+1,7] == x[i,7] + 
            (b*u[i,2]*cos(x[i,7])/cos(x[i,8])+u[i,3]*sin(x[i,7])/cos(x[i,8]))*dt)
        l[i,8]=@NLconstraint(
            m, x[i+1,8] == x[i,8] + 
            (-b*u[i,2]*sin(x[i,7])+u[i,3]*cos(x[i,7]))*dt)
        l[i,9]=@NLconstraint(
            m, x[i+1,9] == x[i,9] + 
            (b*u[i,2]*cos(x[i,7])*tan(x[i,8])+u[i,3]*sin(x[i,7])*tan(x[i,8])+u[i,4])*dt)
    end    
    @NLobjective(m,Min, .5*sum(Q[j]*(x[i,j]-d[i,j])^2 for i=1:N for j=1:n)
                 + .5*sum(R[j]*(u[i,j]^2) for i=1:N for j=1:p)
                 + .5*sum(Q[j]*(x[N+1,j]-d[N+1,j])^2 for j=1:n))

    m[:n] = n; m[:p] = p; m[:nd] = nd; m[:d] = d; m[:l] = l
    
    return m
end
