include("Euler_W.jl")
include("Analyze.jl")

srand(4321)
function sd_2_k(sd)
  return 1./(2*pi*((sd)^2))
end

kee = .26 #WTA
kee = .96 #Normalization
kei = .93
kie = .97
kii = .5

Aee = 84
Aei = 314
Aie = 1319
Aii= 689

s_strength = 3.08
p = .34

N = 4000
Ne = div(N*4, 5)
Ni = div(N, 5)

vth = 20
tau_m = 20.
tau_s = 2.
tau_a = 500.
g_a = .002

kee = sd_2_k(kee)
kei = sd_2_k(kei)
kie = sd_2_k(kie)
kii = sd_2_k(kii)

s1 = s_strength
s2 = s_strength

min_e_neurons = 20
min_i_neurons = 50
runtime = 14000 #ms
h = .1
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
netd_binsize = 150/h
end_trans = 4000
rt = ((ntotal - end_trans)/1000.)*h
half_e = div(Ne, 2)
half_i = div(Ni, 2)
FPe = div(Ne,5)

Aee /= p
Aei /= p
Aie /= p
Aii /= p
#tic
c=0

W = Weights(Ne,Ni,kee,kei,kie,kii,Aee,Aei,Aie,Aii,p)
CSR = sparse_rep(W, N)

@time t, r, Input = euler_lif_CSR(h, runtime, Ne, W, CSR, s1, s2, vth, tau_m, tau_s, tau_a, g_a)

#Analyze Simulations With WTA or Normalization
wta_ness, bias, E_Rates, I_Rates, CV_TOP, CV_BOT, CV_INH, RSC_TOP, RSC_BOT, RSC_INH, ETI, EBI, III, F_TOP, F_BOT, F_INH, te_pt, re_pt, e_top_pt, e_bot_pt, i_pop_pt = WTAN_Analysis(t, r, Input, end_trans, Ne, Ni, half_e, min_e_neurons, min_i_neurons, fbinsize, cbinsize, FPe)
#Analyze Simulations with Rivalry
# te_pt, re_pt, TN, BN, d, cvd, flags, times, cvdlp, t2, f2, cvd2, top, tdom, bot, bdom, nmz, tnmz, FF_TOP, FF_BOT, cwTu, cwTd, cwBu, cwBd, CV_TU, CV_BU, CV_TD, CV_BD, etwm, ebwm, etlm, eblm = Rivalry_Analysis(t, r, Input, end_trans, Ne, Ni, half_e, min_e_neurons, min_i_neurons, fbinsize, cbinsize, netd_binsize, FPe)
