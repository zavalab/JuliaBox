using Revise
using PlasmoData, JLD, MLUtils, LIBSVM, Random, Statistics

### Construct matrices as graphs and get EC curves ###

# Load in the SO2 data; so2_data is size (288, 134, 134, 3)
so2_grayscale_data = load((@__DIR__)*"/../LC_data/so2_grayscale_data.jld")["data"]
so2_classes = load((@__DIR__)*"/../LC_data/so2_classes.jld")["classes"]

diagonal = true
scale = false

# Define threshold range for each EC curve
thresh = 0:.001:1

# Define a matrix for the EC curves; EC curves will be concatenated
so2_ECs = Array{Float64, 2}(undef, length(thresh), 288)

for i in 1:288
    # Build a graph from a 3-D array (134 x 134 x 3)
    mat_graph = matrix_to_graph(so2_grayscale_data[i, :, :], diagonal = diagonal)
    so2_ECs[:, i] = run_EC_on_nodes(mat_graph, thresh, scale = scale)
    #println("Done with iteration $i")
end

### Perform 5 fold CV with SVMs on EC data ###

Random.seed!(10)
# shuffle data
Xs, ys = shuffleobs((so2_ECs, so2_classes))

# define a function for calculating accuracy
function get_accuracy(yhat, ytest)
    num_errors = 0
    for i in 1:length(yhat)
        if yhat[i] != ytest[i]
            num_errors += 1
        end
    end
    return 1 - num_errors/length(yhat)

end

# Perform 5-fold CV
accuracy_values = []
for (train_data, val_data) in kfolds((Xs, ys); k = 5)
    model = svmtrain(train_data[1], train_data[2], kernel = Kernel.Linear)
    yhat, decision_values = svmpredict(model, val_data[1])

    accuracy = get_accuracy(yhat, val_data[2])
    push!(accuracy_values, accuracy)
end

t_processing = @elapsed begin

    # Define threshold range for each EC curve
    thresh = 0:.001:1

    # Define a matrix for the EC curves; EC curves will be concatenated
    so2_ECs = Array{Float64, 2}(undef, length(thresh), 288)

    for i in 1:288
        # Build a graph from a 3-D array (134 x 134 x 3)
        mat_graph = matrix_to_graph(so2_grayscale_data[i, :, :], diagonal = diagonal)
        so2_ECs[:, i] = run_EC_on_nodes(mat_graph, thresh, scale = scale)
        #println("Done with iteration $i")
    end

    ### Perform 5 fold CV with SVMs on EC data ###

    Random.seed!(10)
    # shuffle data
    Xs, ys = shuffleobs((so2_ECs, so2_classes))

    # define a function for calculating accuracy
    function get_accuracy(yhat, ytest)
        num_errors = 0
        for i in 1:length(yhat)
            if yhat[i] != ytest[i]
                num_errors += 1
            end
        end
        return 1 - num_errors/length(yhat)

    end

    # Perform 5-fold CV
    accuracy_values = []
    for (train_data, val_data) in kfolds((Xs, ys); k = 5)
        model = svmtrain(train_data[1], train_data[2], kernel = Kernel.Linear)
        yhat, decision_values = svmpredict(model, val_data[1])

        accuracy = get_accuracy(yhat, val_data[2])
        push!(accuracy_values, accuracy)
    end

    println(accuracy_values)
    println(mean(accuracy_values))
    println(std(accuracy_values))
end

println("Time was $t_processing")
println("diagonal = $diagonal and scale = $scale")
