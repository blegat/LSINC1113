invmod(11, 8)
11 * invmod(11, 8)
invmod(8, 11)
8 * invmod(8, 11)
1 * 11 * invmod(11, 8) + 0 * 8 * invmod(8, 11)

function chinese_rem(a, n)
    # a: list of remainders
    # n: list of moduli (assumed pairwise coprime)
    N = prod(n)
    x = 0
    for (ai, ni) in zip(a, n)
        Ni = div(N, ni)
        # Find the modular inverse of Ni modulo ni
        # That is, find mi such that (Ni * mi) % ni == 1
        mi = invmod(Ni, ni)
        x += ai * Ni * mi
        println("$ai \\cdot $Ni \\cdot $mi + ")
    end
    println("\\equiv $x \\pmod{$N}")
    return mod(x, N)
end

chinese_rem([1, 0], [8, 11])

function prime_factors(n)
    factors = Int[]
    i = 2
    while i * i <= n
        while n % i == 0
            push!(factors, i)
            n = div(n, i)
        end
        i += 1
    end
    if n > 1
        push!(factors, n)
    end
    return factors
end

chinese_rem([0, 1, 2], prime_factors(10759))

prime_factors(225)
prime_factors(365)
prime_factors(687)
prime_factors(4333)
prime_factors(10759)
prime_factors(30685)
prime_factors(60190)
