include("molda_fast.jl")
include("Analyze.jl")

srand(4321)

Aie_NLS = linspace(0,80,10)

Aee = 10.
Aei = 30.
Aie = 40.
Aie_NL = 50.
Aii = 50.

s_strength = 3.08
fe = s_strength
fe2 = fe
fi = 0.
fi2 = fi

# p = .2
k = 800

N = 4000
IFRAC = 2.
N_local = Int64(round(N/2)) #divide the network into two EI circuits
Ni_local = Int64(round(N_local/IFRAC))
Ne_local = Int64(N_local-Ni_local)
Ne2 = Ne_local*2
N2 = Int64(round(N/2)) #divide the network into two EI circuits
NiL = Int64(round(N2/IFRAC))
NeL = Int64(N2-NiL)
Ne2 = NeL*2
Ni2 = NiL*2

vth = 20
tau_m = 20.
tau_s = 2.
tau_a = 1375.
g_a = .4

min_e_neurons = 20
min_i_neurons = 50
runtime = 10000 #ms
h = .1
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
netd_binsize = 50/h
end_trans = 0.
rt = ((ntotal - end_trans)/1000.)*h

W = homogenous_4x4_weights(N, IFRAC, k, Aee, Aei, Aie, Aie_NL, Aii);
CSR = sparse_rep(W, N);

@time t, r = molda_euler_lif_CSR(h, runtime, N, IFRAC, W, CSR, fe, fi, fe2, fi2, vth, tau_m, tau_s, tau_a, g_a)

e_m = find(r .<= Ne2);
i_m = find(r .> Ne2);
te = t[e_m];
re = r[e_m];
ti = t[i_m];
ri = r[i_m];

E_R_top = [length(find(re .== i))/rt for i=NeL+1:Ne2];
E_R_bot = [length(find(re .== i))/rt for i=1:NeL];
I_R_top = [length(find(ri .== i))/rt for i=Ne2+NiL+1:N];
I_R_bot = [length(find(ri .== i))/rt for i=Ne2+1:Ne2+NiL];

MER1 = mean(E_R_bot);#*(1/1000.)
MER2 = mean(E_R_top);#*(1/1000.)
MIR1 = mean(I_R_bot);#*(1/1000.)
MIR2 = mean(I_R_top);#*(1/1000.)

top_e_neurons, bot_e_neurons = Neurons_tb_ns(re, NeL, 10, min_e_neurons);
top_i_neurons, bot_i_neurons = Neurons_tb_ns(ri, NiL, 10, min_i_neurons);

CV_ETOP = CV_ISI_ALLTIME(top_e_neurons, te, re);
CV_EBOT = CV_ISI_ALLTIME(bot_e_neurons, te, re);
CV_ITOP = CV_ISI_ALLTIME(top_i_neurons, ti, ri);
CV_IBOT = CV_ISI_ALLTIME(bot_i_neurons, ti, ri);
