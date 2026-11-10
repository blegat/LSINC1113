function train!(diff, loss, model, W, X, y, η, num_iters, losses = [loss(model(W, X), y)])
	for _ in 1:num_iters
		∇ = diff(loss, model, W, X, y)
        for i in eachindex(W)
            W[i] .= W[i] .- η .* ∇[i]
        end
		push!(losses, loss(model(W, X), y))
	end
	return losses
end
