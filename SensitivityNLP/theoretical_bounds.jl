using SparseArrays

function upsilon_rho(m,ps)
    
    nlp = m.moi_backend.optimizer.model.nlp
    ips = m.moi_backend.optimizer.model.ips

    tmp = copy(nlp.x)
    nlp.x .= 1:nlp.n
    pind = vcat([Int.(value.(p)[:]) for p in ps]...)
    xind = setdiff(1:nlp.n,pind)
    nlp.x .= tmp
    
    jac_sparsity_I = Vector{Int}(undef,nlp.nnz_jac)
    jac_sparsity_J = Vector{Int}(undef,nlp.nnz_jac)
    jac_nz = Vector{Float64}(undef,nlp.nnz_jac)
    
    nlp.jac_sparsity!(jac_sparsity_I,jac_sparsity_J)
    nlp.con_jac!(jac_nz,nlp.x)

    jac = sparse(jac_sparsity_I,jac_sparsity_J,jac_nz,nlp.m,nlp.n)
    # println("Jacobian")
    # println(jac)
    
    hess_sparsity_I = Vector{Int}(undef,nlp.nnz_hess)
    hess_sparsity_J = Vector{Int}(undef,nlp.nnz_hess)
    hess_nz = Vector{Float64}(undef,nlp.nnz_hess)
    
    nlp.hess_sparsity!(hess_sparsity_I,hess_sparsity_J)
    nlp.lag_hess!(hess_nz,nlp.x,nlp.l,1.)

    hess = sparse(hess_sparsity_I,hess_sparsity_J,hess_nz,nlp.n,nlp.n)
    # println("Hessian")
    # println(hess)
    
    kkt = [Symmetric(hess) jac'; jac spzeros(nlp.m,nlp.m)]

    H = kkt[[xind;ips.n+1:end],[xind;ips.n+1:end]]
    R = kkt[[xind;ips.n+1:end],pind]

    # println("H")
    # println(H)
    # println("R")
    # println(R)
    
    factH = svd(Array(H))
    factR = svd(Array(R))
    osigH = factH.S[1]
    usigH = factH.S[end]
    osigR = factR.S[1]

    return osigH*osigR/usigH^2, (osigH^2-usigH^2)/(osigH^2+usigH^2)
end
