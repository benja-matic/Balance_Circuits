###Local Circuit With Random Connections

function random_weights(N, sigma, p)
    W = zeros(N,N)
    for i = 1:N
        W[i,:] = randn(N)*sigma
    end

    px = rand(N,N) .< p
    W = W .* px
    W /= p

    return W
