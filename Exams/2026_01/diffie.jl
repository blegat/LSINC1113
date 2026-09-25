function modinv(a, n)
	g, x, y = gcdx(a, n)
	return mod(x, n)
end

function smallest_primitive_root(p)
    for g in 2:(p-1)
        if length(Set(g.^(0:(p-2)) .% p)) == p-1
            return g
        end
    end
    error("$p has no primitive root")
end

n = 4
p = n^2 + 1
g = smallest_primitive_root(p)
inv_g = modinv(g, p)

gi = g.^(0:(n-1)) .% p
gj = g.^(0:n:(n^2-n)) .% p
inv_gj = inv_g.^(0:n:(n^2-n)) .% p


sort(vec(gi .* gj' .% p))

a = 3
b = 5

A = g^a % p
B = g^b % p

K = A^b % p
K2 = B^a % p

println(K)
println(K2)
