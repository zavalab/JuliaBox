# Created by Sungho Shin (sungho.shin@wisc.edu)

# To reproduce the results in "Graph-Based Modeling and Decomposition of Energy Infrastructures",
# export JULIA_NUM_THREADS=20 before starting the julia session

@__DIR__() in LOAD_PATH || push!(LOAD_PATH,@__DIR__)
using JuMP, Plasmo, MadNLP, AdchemCaseStudy, DelimitedFiles
const ACS = AdchemCaseStudy

# Solvers
with_ma57()=MadNLP.Optimizer(
    blas_num_threads=20,
    disable_garbage_collector=true,
    fixed_variable_treatment="relax_bounds",
    linear_solver="ma57")
with_pardiso_mkl()=MadNLP.Optimizer(
    blas_num_threads=20,
    disable_garbage_collector=true,
    fixed_variable_treatment="relax_bounds",
    linear_solver="pardisomkl",
    pardisomkl_num_threads=20)
with_schwarz()=MadNLP.Optimizer(
    blas_num_threads=20,
    disable_garbage_collector=true,
    fixed_variable_treatment="relax_bounds",
    linear_solver="schwarz",
    schwarz_num_parts=20)
with_schwarz_custom()=MadNLP.Optimizer(
    blas_num_threads=20,
    disable_garbage_collector=true,
    fixed_variable_treatment="relax_bounds",
    linear_solver="schwarz",
    schwarz_custom_partition=true,
    schwarz_num_parts_upper=20)

function benchmark(m::JuMP.Model,optimizer::Function)
    JuMP.set_optimizer(m,optimizer)
    @time optimize!(m)
    empty!(m.moi_backend.model_cache.optattr)
    m.nlp_data != nothing && (m.nlp_data.evaluator = nothing)
    ips = m.moi_backend.optimizer.model.ips
    return ips.cnt.total_time, ips.cnt.linear_solver_time, ips.cnt.eval_function_time
end
function benchmark(m::Plasmo.OptiGraph,optimizer::Function)
    @time MadNLP.optimize!(m;option_dict = optimizer().option_dict)
    ips = m.optimizer
    return ips.cnt.total_time, ips.cnt.linear_solver_time, ips.cnt.eval_function_time
end

# Benchmark
record = []
numvar = []

for Nt in [1,3,7,14,30,60,180]
    # natural gas transient
    for i=1:4
        if i==1
            m = ACS.ng_transient_6a(Nt=Nt,plasmo=false)
        elseif i==2
            m = ACS.ng_transient_6a(Nt=Nt,plasmo=true)
        elseif i==3
            m = ACS.mp_strg_opf_14_ac(Nt=Nt,plasmo=false)
        else
            m = ACS.mp_strg_opf_14_ac(Nt=Nt,plasmo=true)
        end
        
        push!(record,benchmark(m,with_ma57)...); GC.gc()
        push!(record,benchmark(m,with_pardiso_mkl)...); GC.gc()
        push!(record,benchmark(m,m isa OptiGraph ? with_schwarz_custom : with_schwarz)...); GC.gc()
        push!(numvar,num_variables(m)); GC.gc()
    end
end

writedlm("output/adchem_record.csv",record,',')
writedlm("output/adchem_numvar.csv",numvar,',')
