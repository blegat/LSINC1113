function encode(message::String)
    num = zero(BigInt)
    for byte in reverse(Vector{UInt8}(message))
        num <<= 8
        num |= byte
    end
    return num
end

function decode(num::BigInt)
    message = UInt8[]
    while !iszero(num)
        push!(message, convert(UInt8, num & typemax(UInt8)))
        num >>= 8
    end
    return String(message)
end

@assert decode(encode("Hello world")) == "Hello world"

decrypt(encrypted::BigInt, key::BigInt) = decode(encrypted - key)

encrypt(message::String, key::BigInt) = encode(message) + key

key = encode("key")

@assert decrypt(encrypt("Hello world", key), key) == "Hello world"
