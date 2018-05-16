include("molda_fast.jl")
include("Analyze.jl")

srand(4321)

Aie_NLS = linspace(0,80,10)

Aee = 10.5
Aei = 20.
Aie = 30.
Aie_NL = 30.
Aii = 45.

s_strength = 5.08
fe = s_strength
fe2 = fe
fi = 0.#fe - .15
fi2 = fi

# p = .2
k = 200

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
tau_a = 350.
g_a = 0.43

min_e_neurons = 20
min_i_neurons = 50
runtime = 50000 #ms
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

# E_R_top = [length(find(re .== i))/rt for i=NeL+1:Ne2];
# E_R_bot = [length(find(re .== i))/rt for i=1:NeL];
# I_R_top = [length(find(ri .== i))/rt for i=Ne2+NiL+1:N];
# I_R_bot = [length(find(ri .== i))/rt for i=Ne2+1:Ne2+NiL];
#
# MER1 = mean(E_R_bot);#*(1/1000.)
# MER2 = mean(E_R_top);#*(1/1000.)
# MIR1 = mean(I_R_bot);#*(1/1000.)
# MIR2 = mean(I_R_top);#*(1/1000.)
#
# top_e_neurons, bot_e_neurons = Neurons_tb_ns(re, NeL, 10, min_e_neurons);
# top_i_neurons, bot_i_neurons = Neurons_tb_ns(ri, NiL, 10, min_i_neurons);
#
# CV_ETOP = CV_ISI_ALLTIME(top_e_neurons, te, re);
# CV_EBOT = CV_ISI_ALLTIME(bot_e_neurons, te, re);
# CV_ITOP = CV_ISI_ALLTIME(top_i_neurons, ti, ri);
# CV_IBOT = CV_ISI_ALLTIME(bot_i_neurons, ti, ri);


TN, BN = Neurons_tb_ns(re, NeL, 10, 100) #neurons in either pool who fired at least 10 spkes in simulation
ntd, nts = nt_diff_H(te, re, ntotal, NeL, netd_binsize)
s = ntd ./ nts #signal for dominances
flags, times = WLD_01(s, -.333, .333)
top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times, netd_binsize) #find win, lose, and draw times
tbf, rbf = ligase(bot, bdom, te, re, BN) #bottom pool up states
ttf, rtf = ligase(top, tdom, te, re, TN) #top pool up states
tbdf, rbdf = ligase(top, tdom, te, re, BN) #bottom pool down states
ttdf, rtdf = ligase(bot, bdom, te, re, TN) #top pool down states
countFT = count_train_intron(fbinsize, ttf, rtf, TN, length(TN), false)
countFB = count_train_intron(fbinsize, tbf, rbf, BN, length(BN), false)
countFBD = count_train_intron(fbinsize, tbdf, rbdf, BN, length(BN), false)
countFTD = count_train_intron(fbinsize, ttdf, rtdf, TN, length(TN), false)
FF_TOP = fano_train(countFT, -5)
FF_BOT = fano_train(countFB, -5)
FF_TOPD = fano_train(countFTD, -5)
FF_BOTD = fano_train(countFBD, -5)
#correlations
cwTu = rand_pair_cor(cbinsize, ttf, rtf, TN, 1000)
cwBu = rand_pair_cor(cbinsize, tbf, rbf, BN, 1000)
cwBd = rand_pair_cor(cbinsize, ttdf, rtdf, TN, 1000)
cwTd = rand_pair_cor(cbinsize, tbdf, rbdf, BN, 1000)

CV_TU = CV_ISI(top, TN, te, re)
CV_BU = CV_ISI(bot, BN, te, re)
CV_BD = CV_ISI(top, BN, tbdf, rbdf)
CV_TD = CV_ISI(bot, TN, ttdf, rtdf)


d = convert(Array{Float64}, diff(netd_binsize/(1000./h) .* times))
cvd = cv(d)

LP = .3

dx = []
for i in d
    if i > LP
        push!(dx, i)
    end
end
dx = convert(Array{Float64}, dx)
cvdlp = cv(dx)
