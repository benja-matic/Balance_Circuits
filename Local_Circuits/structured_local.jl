###Local Circuit With Structured Connections
###Respects Dale's Law

function von_mises_dist(x, k, mu, N)
  a = exp(k*cos(x-mu))/(N*besseli(0, k))
end

function random_weights(N, p, thresh)
    k = p*N #each cell sees on average k inputs
    ks = sqrt(k)
    W = zeros(N,N)
    for i = 1:N
        W[i,:] = randn(N) * thresh / ks
    end

    px = rand(N,N) .< p
    W = W .* px
    # W /= (p*N/thresh) #synapses shrink with number of connections, grow with threshold

    return W
end

function sparse_rep(W,N)
    flat = []
    for i = 1:N
        mi = find(W[:,i])
        push!(flat, mi)
    end
    return flat
end

function LIF_delta_solve_CSC(W, CSR, S, runtime, h)

    vth = 20
    ntotal = round(Int, runtime/h)
    S = fill(S, N)
    v = rand(N)*vth
    tau_m = 20.
    #dV/dt = S - v/tau_m; Euler: V = V + dV/dt * h; Euler: V += S - V/tau_m * h
    #compute h/tau_m once in advance
    m_leak = h/tau_m
    t = []
    r = []
    input = zeros(N, ntotal)
    volta = zeros(N, ntotal)

    for iter = 1:ntotal

        volta[:, iter] = v
        v += S - v*m_leak
        ves = v .> 20.
        vsm = sum(ves)

        if vsm > 0
            spe = find(ves)
            for j = 1:vsm
                js = spe[j]
                v[CSR[js]] += W[CSR[js], js]
                input[CSR[js], iter] += W[CSR[js], js]
                push!(r, js)
                push!(t, iter)
            end
        end

        v -= 20.*ves
    end

    return t, r, input, volta
end
function LIF_delta_solve(W, S, runtime, h)

    vth = 20
    ntotal = round(Int, runtime/h)
    S = fill(S, N)
    v = rand(N)*vth
    tau_m = 20.
    #dV/dt = S - v/tau_m; Euler: V = V + dV/dt * h; Euler: V += S - V/tau_m * h
    #compute h/tau_m once in advance
    m_leak = h/tau_m
    t = []
    r = []
    input = zeros(N, ntotal)
    volta = zeros(N, ntotal)

    for iter = 1:ntotal

        volta[:, iter] = v
        v += S - v*m_leak
        ves = v .> 20.
        vsm = sum(ves)

        if vsm > 0
            spe = find(ves)
            for j = 1:vsm
                v[CSR[j]] += W[CSR[j], j]
                input[CSR[j], iter] += W[CSR[j],j]
                push!(r, spe[j])
                push!(t, iter)
            end
        end

        v -= 20.*ves
    end

    return t, r, input, volta
end
runtime = 100 #ms
h = 0.1
N = 2000
rt = runtime/1000. #convert to seconds
p = 0.2
S = 2.1
vth = 20. #sticking with 20 but be sure same in W and LIF_delta_solve
W = random_weights(N, p, vth)
CSC = sparse_rep(W,N)
t, r, input, volta = LIF_delta_solve_CSC(W, CSC, S, runtime, h)
