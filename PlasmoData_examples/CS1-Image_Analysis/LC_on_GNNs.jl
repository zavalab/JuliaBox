using Revise
using PlasmoData, Graphs, Statistics, LinearAlgebra, SparseArrays, JLD, Random
using GraphNeuralNetworks, Flux, MLUtils
using Flux:onecold, onehotbatch
using Flux.Losses: crossentropy
using Flux.Losses: logitcrossentropy
using Flux.Data: DataLoader
using CUDA

# I have adapted the code from GraphNeuralNetworks.jl example which can be found at:
# https://github.com/CarloLucibello/GraphNeuralNetworks.jl/blob/master/examples/graph_classification_tudataset.jl
# If you do not have a GPU available, you may need to change the "usecuda" argument in line 70 to false

# The data used in this example is from the paper Bao et al. 2022
# https://doi.org/10.1021/jacs.2c03424

# Define function to give testing and training loss and accuracy
function eval_loss_accuracy(model, data_loader, device)
    loss = 0.
    acc = 0.
    ntot = 0
    for (g, y) in data_loader
        g, y = (g, y) |> device
        n = length(y[1,:])
        ŷ = model(g, g.ndata.x)# |> vec

        y_test = zeros(n)
        for i in 1:n
            argmax_ind     = argmax(y[:, i])[1]
            argmax_ind_hat = argmax(ŷ[:, i])[1]

            if argmax_ind == argmax_ind_hat
                y_test[i] = 1
            end
        end

        loss += logitcrossentropy(ŷ, y; ) * n
        acc  += mean(y_test) * n
        ntot += n
    end
    return (loss = round(loss/ntot, digits=4), acc = round(acc*100/ntot, digits=2))
end

# Define function to load in the appropriate data and build the GNN graphs
function getdataset()
    #data    = load((@__DIR__)*"/SO2_grayscale_data.jld")["data"]
    data = load((@__DIR__)*"/LC_data/so2_data.jld")["data"]
    y_label = load((@__DIR__)*"/LC_data/so2_classes.jld")["classes"]

    graphs = []

    for i in 1:288
        datagraph = matrix_to_graph(data[i, :, :, :])
        gnn_graph = GNNGraph(datagraph.g,
            ndata = get_node_data(datagraph)'
        )
        push!(graphs, gnn_graph)
    end

    y = Flux.onehotbatch(y_label, 1:4)

    return graphs, y
end

# arguments for the `train` function
Base.@kwdef mutable struct Args
    η = 1f-3            # learning rate
    batchsize = 32      # batch size (number of graphs in each batch)
    epochs = 250        # number of epochs
    seed = 17           # set seed > 0 for reproducibility
    usecuda = true      # if true use cuda (if available)
    nhidden = 32        # dimension of hidden features
    infotime = 1        # report every `infotime` epochs
end

# Define function to train the data over each epoch
function train(train_data_overall, args)

    train_data, test_data = splitobs(train_data_overall, at=195, shuffle=true)


    #train_data, test_data = splitobs(dataset, at=NUM_TRAIN, shuffle=true)

    train_loader = DataLoader(train_data; args.batchsize, shuffle=true, collate=true)
    test_loader = DataLoader(test_data; args.batchsize, shuffle=false, collate=true)

    # DEFINE MODEL

    nin = size(train_data[1][1].ndata.x, 1)
    nhidden = args.nhidden

    model = GNNChain(GraphConv(nin => nhidden, relu),
                     GraphConv(nhidden => nhidden, relu),
                     GlobalPool(mean),
                     Dense(nhidden, 4))  |> device

    ps = Flux.params(model)
    opt = Adam(args.η)

    optimal_params = []
    for i in 1:length(ps)
        len = length(ps[i])
        push!(optimal_params, CuArray{Float32, 2}(undef, (len, 1)))# |> device
    end

    # LOGGING FUNCTION

    function report(epoch)
        train = eval_loss_accuracy(model, train_loader, device)
        test = eval_loss_accuracy(model, test_loader, device)
        println("Epoch: $epoch   Train: $(train)   Test: $(test)")
        return test.loss
    end

    # TRAIN
    #report(0)
    test_val_initial = report(0)
    for epoch in 1:args.epochs
        for (g, y) in train_loader
            g, y = (g, y) |> device
            gs = Flux.gradient(ps) do
                ŷ = model(g, g.ndata.x) #|> vec
                logitcrossentropy(ŷ, y)
            end
            Flux.Optimise.update!(opt, ps, gs)
        end

        test_val_next = report(epoch)

        if test_val_next < test_val_initial
            for i in 1:length(ps)
                optimal_params[i][:, 1] .= ps[i][:]
            end
        end

    end

    return optimal_params
end


args = Args()
args.seed > 0 && Random.seed!(args.seed)

if args.usecuda && CUDA.functional()
    device = gpu
    args.seed > 0 && CUDA.seed!(args.seed)
    @info "Training on GPU"
else
    device = cpu
    @info "Training on CPU"
end
dataset = getdataset()
Xs, ys = shuffleobs(dataset)

time_data_processing = @elapsed begin
dataset = getdataset()

Random.seed!(10)
Xs, ys = shuffleobs(dataset)
end
train((Xs, ys), args)

accuracy_values = []

time_CV = @elapsed begin
    for (train_data, test_data) in kfolds((Xs, ys); k = 5)

        println()
        println("STARTING NEXT ITERATION")
        op1 = train(train_data, args)
        trained_model = GNNChain(GraphConv(3 => args.nhidden, relu),
            GraphConv(args.nhidden => args.nhidden, relu),
            GlobalPool(mean),
            Dense(args.nhidden, 4)
        )  |> device

        ps_test = Flux.params(trained_model)

        for i in 1:length(ps_test)
            ps_test[i][:] .= op1[i][:, 1]
        end

        loader = DataLoader(test_data; args.batchsize, shuffle=false, collate=true)

        test_loss_acc = eval_loss_accuracy(trained_model, loader, device)
        push!(accuracy_values, test_loss_acc.acc)
    end
end

println(accuracy_values)
println("The average accuracy is ", mean(accuracy_values), " %")
println("The time to perform 5 fold CV is ", time_CV)
println("The time for data processing is ", time_data_processing)
