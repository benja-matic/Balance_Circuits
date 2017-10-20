###Local Circuit With Random Connections
###Does not respect Dale's law

function random_weights_dale(Ne, Ni, p, thresh, Jee, Jei, Jie, Jii)

    k = p*(Ne+Ni) #each cell sees on average k inputs
    ks = sqrt(k)
    scaling = thresh/ks

    Jee /= Ne*p
    Jei /= Ni*p
    Jie /= Ne*p
    Jii /= Ni*p

    W = zeros(N,N)
    wee = zeros(Ne,Ne)
    wei = zeros(Ne,Ni)
    wie = zeros(Ni,Ne)
    wii = zeros(Ni,Ni)

    for i = 1:Ne
        wee[i,:] = Jee
        wei[i,:] = Jei
    end

    for i = 1:Ni
        wie[i,:] = Jie
        wii[i,:] = Jii
    end

    W[1:Ne, 1:Ne] = wee
    W[1:Ne, Ne+1:end] = wei
    W[Ne+1:end, 1:Ne] = wie
    W[Ne+1:end, Ne+1:end] = wii

    px = rand(N,N) .< p
    W = W .* px

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

function LIF_delta_solve_CSC_dale(W, CSR, S, runtime, h)

    vth = 20
    ntotal = round(Int, runtime/h)
    S = fill(S, N) .* h
    S[Ne+1:end] = 0. #only excitatory neurons receive drive
    v = rand(N)*vth
    tau_m = 20.
    #dV/dt = S - v/tau_m; Euler: V = V + dV/dt * h; Euler: V += S - V/tau_m * h
    #compute h/tau_m once in advance
    m_leak = h/tau_m
    te = []
    re = []
    ti = []
    ri = []
    input = zeros(N, ntotal)
    volta = zeros(N, ntotal)

    for iter = 1:ntotal

        volta[:, iter] = v
        v += S .- v*m_leak
        ves = v .> 20.
        vsm = sum(ves)

        if vsm > 0
            spe = find(ves)
            for j = 1:vsm
                js = spe[j]
                v[CSR[js]] += W[CSR[js], js]
                input[CSR[js], iter] += W[CSR[js], js]
                if js > Ne
                    push!(ti, iter)
                    push!(ri, js)
                else
                    push!(te, iter)
                    push!(re, js)
                end
            end
        end

        v -= 20.*ves
    end

    return te, re, ti, ri, input, volta
end

runtime = 2000 #ms
fbinsize = 100 #ms
cbinsize = 50 #ms
h = 0.05
Ne = 1600
Ni = 400
N = Ne+Ni
rt = runtime/1000. #convert to seconds
p = 0.12
S = 2.1
vth = 20. #sticking with 20 but be sure same in W and LIF_delta_solve

Jee = 80.
Jie = 400.
Jei = -200.
Jii = -100.

W = random_weights_dale(Ne, Ni, p, vth, Jee, Jei, Jie, Jii)
CSC = sparse_rep(W,N)
@time te, re, ti, ri, input, volta = LIF_delta_solve_CSC_dale(W, CSC, S, runtime, h)

figure(1)
plot(te, re, "g.", ms = 1.)
figure(2)
plot(ti, ri, "b.", ms = 1.)

include("analyze_results.jl")

E_Neurons = Neuron_finder(re, 10, 400)
I_Neurons = Neuron_finder(ri, 10, 250)
if ((E_Neurons == -5) | (I_Neurons == -5))
  println("Too Little Spiking")
  return "garbage"
else
E_rate = length(re)/Ne/rt
I_rate = length(ri)/Ni/rt
E_countF = count_train_intron(fbinsize, te, re, E_Neurons, length(E_Neurons), false)
I_countF = count_train_intron(fbinsize, ti, ri, I_Neurons, length(I_Neurons), false)
E_spike_correlation = rand_pair_cor(cbinsize, te, re, E_Neurons, 1000)
I_spike_correlation = rand_pair_cor(cbinsize, ti, ri, I_Neurons, 250)
E_FANO_mean, E_FANO_median, E_FANO_std = fano_train(E_countF, -5)
I_FANO_mean, I_FANO_median, I_FANO_std = fano_train(I_countF, -5)
E_CV_mean, E_CV_median, E_CV_STD = CV_ISI_ALLTIME(E_Neurons, te, re)
I_CV_mean, I_CV_median, I_CV_STD = CV_ISI_ALLTIME(I_Neurons, te, re)

println("FANO E:",E_FANO_mean, ", FANO I: ", I_FANO_mean)
println("CV E:", E_CV_mean, ", CV I: ", I_CV_mean)
println("COR E:", E_spike_correlation, ", COR I: ", I_spike_correlation)
println("rate e: ", E_rate, ", rate i: ", I_rate)
end
