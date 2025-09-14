num_data = 100
num_hidden = 10

# Exemple random
num_features = 10
function random_data(num_data, num_features)
    X = rand(num_data, num_features)
    y = rand(num_data)
    return X, y
end
X, y = random_data(num_data, num_features)

# Moon
import MLJBase, Tables
function random_moon(num_data; noise = 0.1)
    X_table, y_cat = MLJBase.make_moons(num_data, noise = noise)
    X = Tables.matrix(X_table)
    y = 2(float.(y_cat.refs) .- 1.5)
    return X, y
end
X, y = random_moon(num_data)

using Colors
function plot_moon(model, W, X, y)
	col = [Colors.JULIA_LOGO_COLORS.red, Colors.JULIA_LOGO_COLORS.blue]
	scatter(X[:, 1], X[:, 2], markerstrokewidth=0, color = col[round.(Int, (3 .+ y) / 2)], label = "")
    x1 = range(minimum(X[:, 1]), stop = maximum(X[:, 1]), length = 30)
    x2 = range(minimum(X[:, 2]), stop = maximum(X[:, 2]), length = 30)
    contour!(x1, x2, (x1, x2) -> model(W, [x1, x2]')[1], label = "", colorbar_ticks=([1], [0.0]))
end

W1 = rand(size(X, 2), num_hidden)
W2 = rand(num_hidden)
W = (W1, W2)

include(joinpath(@__DIR__, "forward.jl"))

mse(y_est, y) = sum((y_est - y).^2) / length(y)

identity_activation(W, X) = X * W[1] * W[2]
plot_moon(identity_activation, W, X, y)

forward_diff(mse, identity_activation, W, X, y)

include(joinpath(@__DIR__, "train.jl"))

losses = train!(forward_diff, mse, identity_activation, W, X, y, 0.01, 10)

using Plots
plot(eachindex(losses), losses, label = "")

# Train some more

train!(forward_diff, mse, identity_activation, W, X, y, 0.01, 10, losses)
plot(eachindex(losses), losses, label = "")

########## tanh ############

function tanh_activation(W, X)
    W1, W2 = W
    hidden_layer = tanh.(X * W1)
    return hidden_layer * W2
end

plot_moon(tanh_activation, W, X, y)

# On peut réexécuter les lignes qui initialisent `W` ici
# FIXME no method matching `tanh(::Dual)`
forward_diff(mse, tanh_activation, W, X, y)

losses = train!(forward_diff, mse, tanh_activation, W, X, y, 0.01, 10)

# Train more!

train!(forward_diff, mse, tanh_activation, W, X, y, 0.01, 100, losses)

plot(eachindex(losses), losses, label = "")

# ReLU
function relu(x)
    return max(0, x)
end

function relu_activation(W, X)
    W1, W2 = W
    hidden_layer = relu.(X * W1)
    return hidden_layer * W2
end

plot_moon(relu_activation, W, X, y)

# On peut réexécuter les lignes qui initialisent `W` ici
# FIXME no method matching `isless(::Dual, ::Int)`
# 3 choix ici:
# 1) implémenter `isless(::Dual, ::Real)`
# 2) implémenter `max(::Real, ::Dual)`
# 3) implémenter `relu(::Dual)`
forward_diff(mse, relu_activation, W, X, y)

losses = train!(forward_diff, mse, relu_activation, W, X, y, 0.01, 10)

# Train more !

train!(forward_diff, mse, relu_activation, W, X, y, 0.01, 100, losses)

plot(eachindex(losses), losses, label = "")

# Cross entropy
function softmax(x)
    exps = exp.(x .- maximum(x, dims = 2))
    # Le -max(x) est utilisé dans la fonction softmax pour éviter les problèmes 
    # d'instabilité numérique lorsque les entrées sont grandes. Cela améliore la 
    # stabilité des calculs sans affecter les résultats finaux.
    return exps ./ sum(exps, dims = 2)
end

function one_hot_encode(labels::Vector)
    classes = unique!(sort(labels))
    return classes' .== labels
end

Y_encoded = one_hot_encode(y)

function relu_softmax(W, X)
    W1, W2 = W
    hidden_layer = relu.(X * W1)
    logits = hidden_layer * W2
    return softmax(logits)
end

function cross_entropy(Y_est, Y)
    @assert size(Y_est) == size(Y_encoded)
    return -sum(log.(Y_est) .* Y) / size(Y, 1)
end

W = (rand(size(X, 2), num_hidden), rand(num_hidden, size(Y_encoded, 2)))

plot_moon(relu_softmax, W, X, y)

# On peut réexécuter les lignes qui initialisent `W` ici
# FIXME no method matching `isless(::Dual, ::Int)`
forward_diff(cross_entropy, relu_softmax, W, X, Y_encoded)

losses = train!(forward_diff, cross_entropy, relu_softmax, W, X, Y_encoded, 0.01, 10)

# Train more !

train!(forward_diff, cross_entropy, relu_softmax, W, X, Y_encoded, 0.01, 100, losses)
