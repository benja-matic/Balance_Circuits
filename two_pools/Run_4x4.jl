include("4x4.jl")
include("Analyze.jl")

srand(4321)

Aee = 30.
Aei = 400.
Aie = 400.
Aie_NL = 600.
Aii = 200.

s_strength = 3.08
p = .2

N = 4000
N2 = div(N, 2)
NeL = div(4*N2, 5)
Ne2 = NeL*2

vth = 20
tau_m = 20.
tau_s = 2.
tau_a = 575.
g_a = 0


min_e_neurons = 20
min_i_neurons = 50
runtime = 5000 #ms
h = .1
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
netd_binsize = 50/h
end_trans = 2000
rt = ((ntotal - end_trans)/1000.)*h

W = homogenous_4x4_weights(N, p, Aee, Aei, Aie, Aie_NL, Aii);
CSR = sparse_rep(W, N);

@time t, r = euler_lif_CSR_4x4(h, runtime, N, W, CSR, s_strength, vth, tau_m, tau_s, tau_a, g_a)


ex = find(r .<= Ne2)
te_pt = t[ex]
re_pt = r[ex]

ix = find(r .> Ne2)
ti_pt = t[ix]
ri_pt = r[ix]

wta_ness, bias = score_analysis(re_pt, Ne2)
top_e_neurons, bot_e_neurons = Neurons_tb_ns(re_pt, NeL, 10, min_e_neurons)
I_Neurons = Neuron_finder(ri_pt, 10, min_i_neurons)

CV_TOP = CV_ISI_ALLTIME(top_e_neurons, te_pt, re_pt)
CV_BOT = CV_ISI_ALLTIME(bot_e_neurons, te_pt, re_pt)
CV_INH = CV_ISI_ALLTIME(I_Neurons, ti_pt, ri_pt)

# E_R_top = [length(find(re_pt .== i))/rt for i in top_e_neurons]
# E_R_bot = [length(find(re_pt .== i))/rt for i in bot_e_neurons]
E_R_top = [length(find(re_pt .== i))/rt for i=NeL+1:Ne2]
E_R_bot = [length(find(re_pt .== i))/rt for i=1:NeL]
I_R = [length(find(ri_pt .== i))/rt for i=Ne2+1:N]
