using Primes

p = big(Primes.prime(10^6))
g = primitive_root(p)
g_a, g_b, key = generate_keys(g, p, rand(big.(1:(p-2))), rand(big.(1:(p-2))))
