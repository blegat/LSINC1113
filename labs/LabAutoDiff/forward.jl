struct Dual
    value::Float64
    derivative::Float64
end
Base.broadcastable(d::Dual) = Ref(d)
Base.zero(::Dual) = Dual(0, 0)
Base.zero(::Type{Dual}) = Dual(0, 0)
Base.:*(α::Number, x::Dual) = Dual(α * x.value, α * x.derivative)
Base.:*(x::Dual, α::Number) = Dual(x.value * α, x.derivative * α)
Base.:+(x::Dual, y::Dual) = Dual(x.value + y.value, x.derivative + y.derivative)
Base.:-(x::Dual, y::Number) = Dual(x.value - y, x.derivative)
Base.:/(x::Dual, α::Number) = Dual(x.value / α, x.derivative / α)
Base.:*(x::Dual, y::Dual) = Dual(x.value * y.value, x.value * y.derivative + x.derivative * y.value)
Base.:^(x::Dual, n::Integer) = Base.power_by_squaring(x, n)

using OneHotArrays

function forward_diff(loss, model, W::Tuple, X, y, dk, di, dj)
    dW = [
        Dual.(
            W[k],
            if k == dk
                onehot(di, axes(W[k], 1)) * onehot(dj, axes(W[k], 2))'
            else
                0
            end
        )
        for k in eachindex(W)
    ]
    return loss(model(dW, X), y).derivative
end

function forward_diff(loss, model, W::Tuple, X, y, dk)
	return [forward_diff(loss, model, W, X, y, dk, i, j) for i in axes(W[dk], 1), j in axes(W[dk], 2)]
end

function forward_diff(loss, model, W::Tuple, X, y)
	return [forward_diff(loss, model, W, X, y, k) for k in eachindex(W)]
end
