###Local Circuit With Structured Connections
###Respects Dale's Law

function von_mises_dist(x, k, mu, N)
  a = exp(k*cos(x-mu))/(N*besseli(0, k))
end

function weights(Ne,Ni,kee,kei,kie_L,kii,Aee,Aei,Aie_L,Aii,pee,pei,pie,pii)
# Construct weight functions
W = zeros(Ne+Ni, Ne+Ni)
wee = zeros(Ne,Ne)
wei = zeros(Ne,Ni)
wie = zeros(Ni,Ne)
wii = zeros(Ni,Ni)

we = 2*collect(0:Ne-1)*pi/Ne
wi = 2*collect(0:Ni-1)*pi/Ni

for i = 1:Ne
    wee[i,:]=Aee*circshift(von_mises_dist(we, kee, 0, Ne), i-1)
    wei[i,:]=Aei*circshift(von_mises_dist(wi, kei, 0, Ni),div(Ni*i,Ne)-1) #
end

for i = 1:Ni
    wie[i,:]=Aie_L*circshift(von_mises_dist(we, kie_L, 0, Ne),div(Ne*i,Ni)-1)# + Aie_D*circshift(von_mises_dist(we, kie_D, pi, Ne),div(Ne*i,Ni)-1)
    wii[i,:]=Aii*circshift(von_mises_dist(wi, kii, 0, Ni),i-1)
end

# keep connections with probability pab
wee = wee.*(rand(Ne,Ne).<pee)
wei = wei.*(rand(Ne,Ni).<pei)
wie = wie.*(rand(Ni,Ne).<pie)
wii = wii.*(rand(Ni,Ni).<pii)
#println("weights complete")

W[1:Ne, 1:Ne] = wee
W[1:Ne, Ne+1:end] = -wei
W[Ne+1:end, 1:Ne] = wie
W[Ne+1:end, Ne+1:end] = -wii

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
h = .1
tau = 20.
Ne = 1600
Ni = 400
runtime = 10000 #ms
rt = runtime/1000.
ntotal = round(runtime/h) #time points
fbinsize = 4000
cbinsize = 1000

kee = 0.
kei = 0.
kie = 0.
kii = 0.

Aee = 200.
Aei = 1000.
Aie = 1000.
Aii = 500.

p = 0.2
s = 2.1

Aee /= p
Aei /= p
Aie /= p
Aii /= p

vth = 20. #sticking with 20 but be sure same in W and LIF_delta_solve

W = weights(Ne,Ni,kee, kei, kie, kii, Aee, Aei, Aie, Aii, p,p,p,p);
CSC = sparse_rep(W,N)
@time te, re, ti, ri, input, volta = LIF_delta_solve_CSC_dale(W, CSC, S, runtime, h)

include("analyze_results.jl")

include("analyze_results.jl")

E_Neurons = Neuron_finder(re, 10, 400)
I_Neurons = Neuron_finder(ri, 10, 250)
if ((E_Neurons == -5) | (I_Neurons == -5))
  println("Too Little Spiking")
  return "garbage"
else

E_Rates = [length(find(re .== i))/rt for i=1:Ne]
I_Rates = [length(find(ri .== i))/rt for i=1:Ni]
E_rate = length(re)/Ne/rt
I_rate = length(ri)/Ni/rt
E_countF = count_train_intron(fbinsize, te, re, E_Neurons, length(E_Neurons), false)
I_countF = count_train_intron(fbinsize, ti, ri, I_Neurons, length(I_Neurons), false)
E_spike_correlations = rand_pair_cor(cbinsize, te, re, E_Neurons, 1000)
I_spike_correlations = rand_pair_cor(cbinsize, ti, ri, I_Neurons, 250)
E_FANO_mean, E_FANO_median, E_FANO_std = fano_train(E_countF, -5)
I_FANO_mean, I_FANO_median, I_FANO_std = fano_train(I_countF, -5)
ECVS = CV_ISI_ALLTIME(E_Neurons, te, re)
ICVS = CV_ISI_ALLTIME(I_Neurons, ti, ri)

println("FANO E:",E_FANO_mean, ", FANO I: ", I_FANO_mean)
println("CV E:", mean(ECVS), ", CV I: ", mean(ICVS))
println("COR E:", mean(E_spike_correlations), ", COR I: ", mean(I_spike_correlations))
println("rate e: ", E_rate, ", rate i: ", I_rate)
end

figure(1)
plot(te .* h, re, "g.", ms = 1.)
title("Excitatory Raster")
xlabel("Time (ms)")
ylabel("Neuron #")

figure(2)
plot(ti, ri, "b.", ms = 1.)
title("Inhibitory Raster")
xlabel("Time (ms)")
ylabel("Neuron #")

figure(3)
plt[:hist](E_Rates, 100)
title("Histogram of Excitatory Firing Rates")
xlabel("Firing Rate (Log(Hz))")

figure(4)
plt[:hist](ECVS, 100)
title("Histogram of CV(isi) in Excitatory Cells")
xlabel("CV(isi)")

figure(5)
plt[:hist](E_spike_correlations, 100)
title("Histogram of Pairwise Spike Count Correlations")
xlabel("Correlation")

figure(6)
plot(W[500,:][:])
title("Sample Row of Weights Matrix")
xlabel("Synaptic Strength Received by This Neuron (mV)")

figure(7)
plt[:hist](input[500,:][:], 100)
title("Histogram of Synaptic Input to Sample Neuron")
xlabel("Synaptic Input/dt (mV/$(h) ms)")
