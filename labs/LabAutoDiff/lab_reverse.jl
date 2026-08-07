include("lab.jl")
using Test, LinearAlgebra

X, y = random_moon(num_data)
W1 = rand(size(X, 2), num_hidden)
W2 = rand(num_hidden)
W = (W1, W2)

include(joinpath(@__DIR__, "reverse.jl"))

∇f = @time forward_diff(mse, identity_activation, W, X, y)
∇r = @time reverse_diff(mse, identity_activation, W, X, y)

# We should get a difference at the order of `1e-15` unless we got it wrong:
norm.(∇f .- ∇r)
@test all(∇f .≈ ∇r)


∇f = @time forward_diff(mse, tanh_activation, W, X, y)
∇r = @time reverse_diff(mse, tanh_activation, W, X, y)

# We should get a difference at the order of `1e-15` unless we got it wrong:
norm.(∇f .- ∇r)
@test all(∇f .≈ ∇r)

∇f = @time forward_diff(mse, relu_activation, W, X, y)
∇r = @time reverse_diff(mse, relu_activation, W, X, y)

norm.(∇f .- ∇r)
@test all(∇f .≈ ∇r)

Y_encoded = one_hot_encode(y)
W = (rand(size(X, 2), num_hidden), rand(num_hidden, size(Y_encoded, 2)))

∇f = @time forward_diff(cross_entropy, relu_softmax, W, X, Y_encoded)
∇r = @time reverse_diff(cross_entropy, relu_softmax, W, X, Y_encoded)

norm.(∇f .- ∇r)
@test all(∇f .≈ ∇r)
