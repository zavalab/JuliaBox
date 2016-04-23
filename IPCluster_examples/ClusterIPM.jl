# Julia interface for ClusterIPM
# Nai-Yuan Chiang && Kibaek Kim - 2015

# Check if JuMP and StochJuMP are installed.
using JuMP

using StochJuMP

macro checked_lib(libname, path)
    (Libdl.dlopen_e(path) == C_NULL) && error("Unable to load \n\n$libname ($path)\n\n")
    quote const $(esc(libname)) = $path end
end

# Load dependencies
@checked_lib libcoinhsl "/opt/IPCluster/src/ThirdParty/coinhsl/lib/lib/libcoinhsl.so"
@checked_lib libclusteripm "/opt/IPCluster/lib/libClusterIPM.so.1.0"

macro CluIPM_ccall(func, args...)
	@unix_only return quote
		ccall(($func, "/opt/IPCluster/lib/libClusterIPM.so.1.0"), $(args...))
	end
end

full=0;
structure=1;
shur=2;

decrease=0; 
dynamic=1;

vcint(a::Vector{Int}) = Base.convert(Vector{Cint},a)


type CluIPM
	p::Ptr{Void}
	function CluIPM()
		p = @CluIPM_ccall("createClusterIPM", Ptr{Void}, ())
		ipm_prob = new(p)
		finalizer(ipm_prob, freeCluIPM)
		return ipm_prob
	end
end

function freeCluIPM(ipm_prob::CluIPM)
	if ipm_prob.p == C_NULL
		return
	end
	@CluIPM_ccall("freeClusterIPM", Void, (Ptr{Void},), ipm_prob.p)
	ipm_prob.p = C_NULL
	return
end

function checkConstraintTypes(m::JuMP.Model)
    numRows = length(m.linconstr)
    eq_idx   = Int[]
    sizehint!(eq_idx, numRows)
    ineq_idx = Int[]
    sizehint!(ineq_idx, numRows)
    for it in 1:numRows
        if m.linconstr[it].lb == m.linconstr[it].ub
            push!(eq_idx, it)
        else
            push!(ineq_idx, it)
        end
    end
	@assert length(eq_idx)   == numRows
    @assert length(ineq_idx) == 0
end

function get_sparse_Hes(m::JuMP.Model)
    n = m.numCols
    vars1 = Int[x.col for x in m.obj.qvars1]
    vars2 = Int[x.col for x in m.obj.qvars2]
    coeff_copy = copy(m.obj.qcoeffs)
    for i in 1:length(vars1)
        if vars1[i] == vars2[i] # "terms" form
            coeff_copy[i] *= 2
        end
        if vars1[i] > vars2[i]
            vars1[i], vars2[i] = vars2[i], vars1[i]
        end
    end
    Q = sparse((vars1), (vars2), coeff_copy, n, n)
    istriu(Q) && (Q == Q')
    #@assert istril(Q)
    irow_hes, jcol_hes, val_hes = findnz(Q)
    return (irow_hes), (jcol_hes), val_hes
end


function prepConstrMatrix(m::JuMP.Model)
    if !haskey(m.ext, :Stochastic)
        return JuMP.prepConstrMatrix(m)
    end

    stoch = StochJuMP.getStochastic(m)
    if stoch.parent == nothing
       return JuMP.prepConstrMatrix(m)
    else
        rind = Int[]
        cind = Int[]
        value = Float64[]
        linconstr = deepcopy(m.linconstr)
        for (nrow,con) in enumerate(linconstr)
            aff = con.terms
            for (var,id) in zip(reverse(aff.vars), length(aff.vars):-1:1)
                push!(rind, nrow)
                if m.linconstr[nrow].terms.vars[id].m == stoch.parent
                    push!(cind, var.col)
                elseif m.linconstr[nrow].terms.vars[id].m == m
                    push!(cind, stoch.parent.numCols + var.col)
                end
                push!(value, aff.coeffs[id])
                splice!(aff.vars, id)
                splice!(aff.coeffs, id)
            end
        end
    end
    return sparse(rind, cind, value, length(m.linconstr), stoch.parent.numCols + m.numCols)
end


function getDataFormat(model::Model)

	# objective coefficients
	obj, rlb, rub = JuMP.prepProblemBounds(model)

	# do test to check if rlb = rub
	#checkConstraintTypes(model)

	# Get a column-wise sparse matrix for constraint coeff matrix (Jacobian)
	jac = prepConstrMatrix(model)
	
	# sparse description of constraint coeff matrix (Jacobian)
	irow_jac, jcol_jac, val_jac = findnz(jac)
	irow_jac_CIndex = convert(Vector{Cint}, irow_jac )
	jcol_jac_CIndex = convert(Vector{Cint}, jcol_jac )
	nnz_jac = convert(Cint, length(irow_jac))
		
	# get sparse description of hessian
	irow_hes, jcol_hes, val_hes = get_sparse_Hes(model)
	irow_hes_CIndex = convert(Vector{Cint}, irow_hes )
	jcol_hes_CIndex = convert(Vector{Cint}, jcol_hes )
	nnz_hes = convert(Cint, length(irow_hes))

	return obj, model.colLower, model.colUpper, rlb, nnz_jac, (irow_jac_CIndex), (jcol_jac_CIndex), val_jac, nnz_hes, (irow_hes_CIndex), (jcol_hes_CIndex), val_hes
end

function get_subJac_data(owner::JuMP.Model, master::JuMP.Model)

	# objective coefficients
	#obj, rlb, rub = JuMP.prepProblemBounds(owner)
	# do test to check if rlb = rub
	#checkConstraintTypes(owner)

	# Get sparse description of constraint coeff matrix (Jacobian)
	#jac = prepConstrMatrix(owner)
	#irow_jac, jcol_jac, val_jac = findnz(jac)
	#nnz_jac = convert(Cint, length(irow_jac))

	  
	# split jacobian into TEz and TEw
    	#numRows = convert(Int, length(rlb))
	numRows =  length(owner.linconstr)
    	rowptr  = Array(Int, numRows+1)

    # get a vague idea of how large submatrices will be
    #nnz = nnz_jac
    nnz = 0
    for it in 1:numRows
        nnz += length(owner.linconstr[it].terms.coeffs)
    end

    rowval = Int[]
    colval   = Int[]
    sizehint!(colval, nnz)
    rownzval = Float64[]
    sizehint!(rownzval, nnz)
    
    nnz = 0
    tmprow   = JuMP.IndexedVector(Float64, master.numCols)
    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    for it in 1:numRows
        rowptr[it] = nnz + 1
        coeffs = owner.linconstr[it].terms.coeffs
        vars = owner.linconstr[it].terms.vars
        for (ct,ind) in enumerate(coeffs)
            if vars[ct].m == master
                JuMP.addelt!(tmprow, vars[ct].col, ind)
            end
        end
        for i in 1:tmprow.nnz
            nnz += 1
            idx = tmpnzidx[i]
	    push!(rowval, it)
            push!(colval, idx)
            push!(rownzval, tmpelts[idx])
        end
        JuMP.empty!(tmprow)
    end
    rowptr[numRows+1] = nnz + 1

     irow_mat_CIndex = convert(Vector{Cint}, rowval )
     jcol_mat_CIndex = convert(Vector{Cint}, colval )
     nnz_mat = convert(Cint, length(rowval))
    return nnz_mat, (irow_mat_CIndex), (jcol_mat_CIndex), rownzval

   #=
    nnz = 0
    tmprow   = JuMP.IndexedVector(Float64, master.numCols)
    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    for it in 1:numRows
        rowptr[it] = nnz + 1
        coeffs = owner.linconstr[it].terms.coeffs
        vars = owner.linconstr[it].terms.vars
        for (it,ind) in enumerate(coeffs)
            if vars[it].m == master
                JuMP.addelt!(tmprow, vars[it].col, ind)
            end
        end
        for i in 1:tmprow.nnz
            nnz += 1
            idx = tmpnzidx[i]
            push!(colval, idx)
            push!(rownzval, tmpelts[idx])
        end
        JuMP.empty!(tmprow)
    end
    rowptr[numRows+1] = nnz + 1

    mat = SparseMatrixCSC(master.numCols, numRows, rowptr, colval, rownzval)
	# sparse description of this matrix (sub-mat of Jacobian)
	irow_mat, jcol_mat, val_mat = findnz(mat)
	irow_mat_CIndex = convert(Vector{Cint}, irow_mat )
	jcol_mat_CIndex = convert(Vector{Cint}, jcol_mat )
	nnz_mat = convert(Cint, length(irow_mat))

    return nnz_mat, (irow_mat_CIndex), (jcol_mat_CIndex), val_mat
    =#
end

function Set_NumScenario(ipm_prob::CluIPM, numScens::Integer)
	@CluIPM_ccall("Set_NumScenario", Void, (Ptr{Void}, Cint), ipm_prob.p, convert(Cint, numScens))
end

function CluIPMloadProblem(ipm_prob::CluIPM, master::Model)

	# get 1st stage problem
	stoch  = getStochastic(master)
	
	nscen  = convert(Cint, stoch.num_scen)
	ncols1 = convert(Cint, master.numCols)
	nrows1 = convert(Cint, length(master.linconstr))

	Set_NumScenario(ipm_prob, nscen)
	
	# get 1st stage problem data
        obj, xlb, xub, cb, nnz_jac,irow_jac, jcol_jac, val_jac, nnz_hes, irow_hes, jcol_hes, val_hes = getDataFormat(master)

	@CluIPM_ccall("Set_FirstStageData", Void, 
	  (  Ptr{Void}, Cint, Cint, 
		 Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
		 Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
		 Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}		),
		 ipm_prob.p, ncols1, nrows1, 
		 cb, obj, xlb, xub, 
		 nnz_jac, irow_jac, jcol_jac, val_jac,
		 nnz_hes, irow_hes, jcol_hes, val_hes	)


	# get 2nd stage fixed data from first scenario
	ncols2 = 0
	nrows2 = 0
	for s in 1:nscen
	    child=stoch.children[s]
	    ncols2 = convert(Cint, child.numCols)
	    nrows2 = convert(Cint, length(child.linconstr))
	    xlb2=child.colLower
	    xub2=child.colUpper
            obj2, rlb2, rub2 = JuMP.prepProblemBounds(child)
	    nnz_jac_c1r2, irow_jac_c1r2, jcol_jac_c1r2, val_jac_c1r2 = get_subJac_data(child, master)
	    nnz_jac_c2r2, irow_jac_c2r2, jcol_jac_c2r2, val_jac_c2r2 = get_subJac_data(child, child)
	    irow_hes2, jcol_hes2, val_hes2 = get_sparse_Hes(child)
	    irow_hes2_CIndex = convert(Vector{Cint}, irow_hes2 )
	    jcol_hes2_CIndex = convert(Vector{Cint}, jcol_hes2 )
	    nnz_hes2 = convert(Cint, length(irow_hes2))


		
	    @CluIPM_ccall("Set_SecondStageData_noRhs", Void,
	     (  Ptr{Void}, Cint, Cint, 
		 Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
		 Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
		 Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
		 Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint),
		 ipm_prob.p, ncols2, nrows2, 
		 obj2, xlb2, xub2, 
		 nnz_jac_c1r2, irow_jac_c1r2, jcol_jac_c1r2, val_jac_c1r2,
		 nnz_jac_c2r2, irow_jac_c2r2, jcol_jac_c2r2, val_jac_c2r2,
		 nnz_hes2, vcint(irow_hes2), vcint(jcol_hes2), val_hes2,convert(Cint, s-1) 	)
#        end

	# get 2nd stage data
#	for s in 1:nscen
		# get model
#		child=stoch.children[s]
#		@assert ncols2   == convert(Cint, child.numCols)
#		@assert nrows2   == convert(Cint, length(child.linconstr))

		# get model data
		xlb2=child.colLower
		xub2=child.colUpper
#		obj, rlb, rub = JuMP.prepProblemBounds(child)
#		checkConstraintTypes(child)

		@CluIPM_ccall("Set_SecondStageData_rhs", Void,
		  ( Ptr{Void}, Cint, Cint, 
			Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint),
			ipm_prob.p, ncols2, nrows2, 
		 	rlb2, xlb2, xub2, convert(Cint, s-1) )
	end

end

function SetUpCluIPM(ipm_prob::CluIPM, comm)
	@CluIPM_ccall("SetUp_ClusterIPM", Void, (Ptr{Void},MPI.Comm), ipm_prob.p, comm)
end

function solveCluIPM(ipm_prob::CluIPM, comm)
	@CluIPM_ccall("ClusterIPM_Solve", Void, (Ptr{Void}, MPI.Comm), ipm_prob.p, comm)
end

function Set_Algorithm(ipm_prob::CluIPM, alg::Integer)
	@CluIPM_ccall("Set_Algorithm", Void, (Ptr{Void}, Cint), ipm_prob.p, convert(Cint, alg))
end

function Set_Par_mu_update(ipm_prob::CluIPM, muAlg::Integer)
	@CluIPM_ccall("Set_Par_mu_update", Void, (Ptr{Void}, Cint), ipm_prob.p, convert(Cint, muAlg))
end

#function getObjValue(ipm_prob::CluIPM)
#	return @CluIPM_ccall("getObjValue", Cdouble, (Ptr{Void},), ipm_prob.p)
#end



function CluIPM_solve(master::JuMP.Model)
         MPI.Init()
         comm = MPI.COMM_WORLD
         println("[$(MPI.Comm_rank(comm))/$(MPI.Comm_size(comm))] create problem ")
	 prob = CluIPM()
	 CluIPMloadProblem(prob, m)
         SetUpCluIPM(prob, comm)
	 solveCluIPM(prob, comm)
         MPI.Finalize()
         return Int32(1)
end










