using JuMP, MadNLP, Plots, LaTeXStrings; pgfplotsx()

include("cstudy_include.jl")

ns= 30
sig = .1

for (name,N,j,get_model) = [
    ("linear",100,50,get_model_linear),
    ("nonlinear",100,50,get_model_nonlinear)
]
    println("Solving $name problem")
    
    for (case,q,b,color) in [(1,1,1,:blue),(2,0,0,:red)]
        println("... for q=$q, b=$b")
        
        m = get_model(N,q,b)
        println("...... reference problem")
        optimize!(m)
        
        xsol = value.(m[:x])
        usol = value.(m[:u])
        lsol = dual.(m[:l])

        xsam = []
        usam = []
        lsam = []
        dref = value.(m[:d][j,:])

        println("...... perturbed problems")
        for i=1:ns
            println("......... sample #$i")
            set_value.(m[:d][j,:],dref+sig*randn(m[:nd]))
            optimize!(m)
            push!(xsam, value.(m[:x]))
            push!(usam, value.(m[:u]))
            push!(lsam, dual.(m[:l]))
        end

        mkpath("fig")
        println("...... plotting")
        for (var,sym,sol,sam,l) in  [("x","x",xsol,xsam,m[:n]),
                                     ("u","u",usol,usam,m[:p]),
                                     ("l","\\lambda",lsol,lsam,m[:n])]
            for i=1:l
                println("......... var $var[$i]")
                plt = plot(xlabel=L"$i$",ylabel="\$$(sym)_i[$i]\$",xlim=(0,N),
                           size=(450,150),legend=:none,framestyle=:box);
                for j=1:ns
                    plot!(plt,0:N-1,sam[j][1:N,i],color=:lightgray);
                end
                plot!(plt,0:N-1,sol[1:N,i],color=color,linestyle=:dash);
                plot!(plt,[j-1],linetype=:vline,linestyle=:dot,color=:black);
                savefig(plt,"fig/cstudy-$name-case-$case-$var-$i.pdf")
            end
        end
        println()
    end
end
