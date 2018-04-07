include("Euler_W.jl")
include("Analyze.jl")

srand(4321)
function sd_2_k(sd)
  return 1./(2*pi*((sd)^2))
end

kee = .26 #WTA
# kee = .96 #Normalization
kei = .93
kie = .97
kii = .5

# Aee = 210.
# Aee = 120.
# Aie = 3500.
# Aei = 564.
# Aie = 5000.
# Aii= 689

Aee = 84
Aei = 314
Aie = 1319
Aii = 739.
# kee = .48#parse(Float64, ARGS[1])
# kei = .93#parse(Float64, ARGS[2])
# kie = .97#parse(Float64, ARGS[3])
# kii = .9#parse(Float64, ARGS[4])
#
# Aee = 54#piarse(Float64, ARGS[5])
# Aei = 104#parse(Float64, ARGS[6])
# Aie = 180#parse(Float64, ARGS[7])
# Aii= 150#parse(Float64, ARGS[8])

kee = sd_2_k(kee)
kei = sd_2_k(kei)
kie = sd_2_k(kie)
kii = sd_2_k(kii)

s_strength = 3.08
p = .34

N = 4000
Ne = div(N*4, 5)
Ni = div(N, 5)

vth = 20
tau_m = 20.
tau_s = 2.
tau_a = 625.
# tau_a = 425.
# g_a = .001
g_a = .002
# g_a = .0075
# kee = sd_2_k(kee)
# kei = sd_2_k(kei)
# kie = sd_2_k(kie)
# kii = sd_2_k(kii)

s1 = s_strength
s2 = s_strength

min_e_neurons = 20
min_i_neurons = 50
runtime = 100000 #ms
h = .1
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
netd_binsize = 50/h
end_trans = 2000
rt = ((ntotal - end_trans)/1000.)*h
half_e = div(Ne, 2)
half_i = div(Ni, 2)
FPe = div(Ne,5)
# Aie = [500, 1000, 1500, 2000, 5000]

# Angle = 500.
ang = 0

Aee /= p
Aei /= p
Aie /= p
Aii /= p
#tic
c=0

# alts = []
# MDTs = []
# MDBs = []

# for i in Aie

w2 = 5
input_width = .15 #units of k

W = Weights(Ne,Ni,kee,kei,kie,kii,Aee,Aei,Aie,Aii,p)
CSR = sparse_rep(W, N)
@time t, r = euler_lif_CSR(h, runtime, Ne, W, CSR, s1, s2, w2, input_width, vth, tau_m, tau_s, tau_a, g_a, ang)

em = find(r .<= Ne)
re = r[em]
te = t[em]

ntd, nts = nt_diff_H(te, re, ntotal, half_e, netd_binsize)
s = ntd ./ nts #signal for dominances
flags, times = WLD_01(s, -.333, .333)

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

t2, f2 = splice_reversions(flags, times) ###Another way to get rid of rapid switches that aren't really there
fw = find(f2 .== "win")
fl = find(f2 .== "lose")
tx = t2 .* netd_binsize
d = convert(Array{Float64}, diff(tx))
cvd2 = cv(d) ### Another estimate of CVD

TN, BN = Neurons_tb_ns(re, half_e, 10, 100) #neurons in either pool who fired at least 10 spkes in simulation
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
