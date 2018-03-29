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

Aee = 84
Aei = 314
Aie = 1319
Aii= 689

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
tau_a = 575.
# tau_a = 425.
# g_a = 0
# g_a = .001
g_a = .002

# kee = sd_2_k(kee)
# kei = sd_2_k(kei)
# kie = sd_2_k(kie)
# kii = sd_2_k(kii)

s1 = s_strength
s2 = s_strength

min_e_neurons = 20
min_i_neurons = 50
runtime = 20000 #ms
h = .1
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
netd_binsize = 50/h
end_trans = 4000
rt = ((ntotal - end_trans)/1000.)*h
half_e = div(Ne, 2)
half_i = div(Ni, 2)
FPe = div(Ne,5)
Angle = 200

Aee /= p
Aei /= p
Aie /= p
Aii /= p
#tic
c=0

W = Weights(Ne,Ni,kee,kei,kie,kii,Aee,Aei,Aie,Aii,p)
CSR = sparse_rep(W, N)

# @time t, r, Input, adapt, Exc, Inh, Vlt = euler_lif_CSR(h, runtime, Ne, W, CSR, s1, s2, vth, tau_m, tau_s, tau_a, g_a, angle)
@time t, r = euler_lif_CSR(h, runtime, Ne, W, CSR, s1, s2, vth, tau_m, tau_s, tau_a, g_a, Angle)

ex = find(r .<= Ne)
te = t[ex]
re = r[ex]


FPe = div(Ne,5)
P2s = round(Int,3*Ne/4-FPe) - Angle
P2e = round(Int,3*Ne/4+FPe) - Angle

ntd, nts = nt_diff_angle(te, re, ntotal, P2s, P2e, netd_binsize)
s = ntd ./nts
flags, times = WLD_01(s, -.333, .333)
top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times, netd_binsize) #find win, lose, and draw times
MDT = tdom/length(top)
MDB = bdom/length(bot)
MDN = tnmz/length(nmz)
MDD = ntotal/length(times)
println("$(MDT), $(MDB), $(length(times))")

plot(t, r, "g.", ms = 1.)
plot(times .* netd_binsize, fill(P2s, length(times)), "r.")
#
# ntd, nts = nt_diff_H(te, re, ntotal, half_e, netd_binsize)
# s = ntd ./ nts #signal for dominances
# flags, times = WLD_01(s, -.333, .333)
#
d = convert(Array{Float64}, diff(netd_binsize/(1000./h) .* times))
cvd = cv(d) ###Raw estimate of CVD, likely to include very rapid switches which should really be smoothed out

LP = .3

dx = []
for i in d
    if i > LP
        push!(dx, i)
    end
end
dx = convert(Array{Float64}, dx)
cvdlp = cv(dx) ###Low-pass filter measure of CVD

###Use this code for spiking statistics, which includes only dominance and suppression times; mixed percept time is absorbed into either of those
# t2, f2 = splice_reversions(flags, times) ###Another way to get rid of rapid switches that aren't really there
# fw = find(f2 .== "win")
# fl = find(f2 .== "lose")
# tx = t2 .* netd_binsize
# d2 = convert(Array{Float64}, diff(tx))
# cvd2 = cv(d2) ### Another estimate of CVD


#Analyze Simulations With WTA or Normalization
# wta_ness, bias, E_Rates, I_Rates, CV_TOP, CV_BOT, CV_INH, RSC_TOP, RSC_BOT, RSC_INH, ETI, EBI, III, FF_TOP, FF_BOT, FF_INH, te_pt, re_pt, e_top_pt, e_bot_pt, i_pop_pt = WTAN_Analysis(t, r, Input, end_trans, Ne, Ni, half_e, min_e_neurons, min_i_neurons, fbinsize, cbinsize, FPe)

# write_array("N_TOP_IM.txt", ETI)
# write_array("N_BOT_IM.txt", EBI)
# write_array("N_TOP_RSC.txt", RSC_TOP)
# write_array("N_BOT_RSC.txt", RSC_BOT)
# write_array("N_TOP_FF.txt", FF_TOP)
# write_array("N_BOT_FF.txt", FF_BOT)
# write_array("N_TOP_CV.txt", CV_TOP)
# write_array("N_BOT_CV.txt", CV_BOT)

# write_array("W_TOP_IM.txt", ETI)
# write_array("W_BOT_IM.txt", EBI)
# write_array("W_TOP_RSC.txt", RSC_TOP)
# write_array("W_BOT_RSC.txt", RSC_BOT)
# write_array("W_TOP_FF.txt", FF_TOP)
# write_array("W_BOT_FF.txt", FF_BOT)
# write_array("W_TOP_CV.txt", CV_TOP)
# write_array("W_BOT_CV.txt", CV_BOT)



#Analyze Simulations with Rivalry
# te_pt, re_pt, TN, BN, d, cvd, flags, times, cvdlp, t2, f2, cvd2, top, tdom, bot, bdom, nmz, tnmz, FF_TOP, FF_BOT, cwTu, cwTd, cwBu, cwBd, CV_TU, CV_BU, CV_TD, CV_BD, etwm, ebwm, etlm, eblm = Rivalry_Analysis(t, r, Input,adapt, end_trans, Ne, Ni, half_e, min_e_neurons, min_i_neurons, fbinsize, cbinsize, netd_binsize, FPe)

# write_array("R_TOP_W_IM.txt", ETI)
# write_array("R_BOT_W_IM.txt", EBI)
# write_array("R_TOP_W_RSC.txt", cwTu)
# write_array("R_BOT_W_RSC.txt", cwBu)
# write_array("R_TOP_W_FF.txt", FF_TOP)
# write_array("R_BOT_W_FF.txt", FF_BOT)
# write_array("R_TOP_W_CV.txt", CV_TU)
# write_array("R_BOT_W_CV.txt", CV_BU)

# write_array("R_TOP_L_IM.txt", ETI)
# write_array("R_BOT_L_W_IM.txt", EBI)
# write_array("R_TOP_L_RSC.txt", cwTd)
# write_array("R_BOT_L_RSC.txt", cwBd)
# write_array("R_TOP_L_FF.txt", FF_TOP)
# write_array("R_BOT_L_FF.txt", FF_BOT)
# write_array("R_TOP_L_CV.txt", CV_TD)
# write_array("R_BOT_L_CV.txt", CV_BD)
