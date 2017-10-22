###Local Circuit With Random Connections
###Does not respect Dale's law

function random_weights(N, p, thresh)
    k = p*N #each cell sees on average k inputs
    ks = sqrt(k)
    W = zeros(N,N)

    for i = 1:N
        W[i,:] = -1.5 + 2.5*randn(N) * thresh / ks
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

function interpolate_spike(v2, v1, vth)
    x = (v1 - v2)
    t = (vth - v2)/x
    return t
end

function LIF_delta_solve_CSC(W, CSR, S, runtime, h)

    vth = 20
    ntotal = round(Int, runtime/h)
    S = fill(S, N) .* h
    synapse = zeros(N)
    v = rand(N)*vth
    v_buff = v
    tau_m = 20.
    tau_s = 2.
    #dV/dt = S - v/tau_m; Euler: V = V + dV/dt * h; Euler: V += S - V/tau_m * h
    #compute h/tau_m once in advance
    m_leak = h/tau_m
    s_leak = h/tau_s
    t = []
    r = []
    input = zeros(N, ntotal)
    volta = zeros(N, ntotal)

    for iter = 1:ntotal

        volta[:, iter] = v
        v += S + synapse - v*m_leak
        input[:, iter] = synapse
        ves = v .> 20.
        vsm = sum(ves)

        if vsm > 0
            spe = find(ves)
            for j = 1:vsm
                #delta_h = interpolate_spike(v[spe[j]], v_buff[spe[j]], vth)
                #ls = exp(delta_h/tau_s)
                js = spe[j]
                synapse[CSR[js]] += W[CSR[js], js] #.* ls
                push!(t, iter)
                push!(r, js)
            end
        end

        v -= 20.*ves
        v_buff = v
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
                v += W[:, spe[j]]
                input[:, iter] += W[:,spe[j]]
                push!(r, spe[j])
                push!(t, iter)
            end
        end

        v -= 20.*ves
    end

    return t, r, input, volta
end

runtime = 10000 #ms
fbinsize = 100 #ms
cbinsize = 50 #ms
h = 0.05
N = 2000
rt = runtime/1000. #convert to seconds
p = 0.0
S = 2.0
vth = 20. #sticking with 20 but be sure same in W and LIF_delta_solve

W = random_weights(N, p, vth)
CSC = sparse_rep(W,N)
@time t, r, input, volta = LIF_delta_solve_CSC(W, CSC, S, runtime, h)

include("analyze_results.jl")

tic()
Neurons = Neuron_finder(r, 10, 400)
Rates = [length(find(r .== i))/rt for i=1:N]
countF = count_train_intron(fbinsize, t, r, Neurons, length(Neurons), false)
spike_correlations = rand_pair_cor(cbinsize, t, r, Neurons, 1000)
FANO_mean, FANO_median, FANO_std = fano_train(countF, -5)
CVS = CV_ISI_ALLTIME(Neurons, t, r)
Input_M = [mean(input[i,:]) for i=1:N]
Input_V = [var(input[i,:]) for i=1:N]

println("FANO MEAN:",FANO_mean)
println("CV MEAN:", mean(CVS))
println("COR MEAN:", mean(spike_correlations))
println("RATE MEAN: ", mean(Rates))
println("INPUT MEAN: ", mean(Input_M))
println("INPUT VAR: ", mean(Input_V))

toc()

figure(1)
plot(t .* h, r, "g.", ms = 1.)
title("Raster Plot")
xlabel("Time (ms)")
ylabel("Neuron #")

figure(2)
plt[:hist](Rates, 100)
title("Histogram of Firing Rates")
xlabel("Firing Rate (Hz)")

figure(3)
plt[:hist](CVS, 100)
title("Histogram of CV(isi)")
xlabel("CV(isi)")

figure(4)
plt[:hist](spike_correlations, 100)
title("Histogram of Pairwise Spike Count Correlations")
xlabel("Correlation")

figure(5)
plot(W[500,:][:])
title("Sample Row of Weights Matrix")
xlabel("Synaptic Strength Received by This Neuron (mV)")

figure(6)
plt[:hist](input[500,:][:], 100)
title("Histogram of Synaptic Input to Sample Neuron")
xlabel("Synaptic Input/dt (mV/$(h) ms)")
