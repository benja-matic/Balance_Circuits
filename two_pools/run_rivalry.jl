#println("GOT THIS FAR")
#include("measure_balance.jl")
include("euler_b_2_rivalry.jl")
srand(4321)
function sd_2_k(sd)
  return 1./(2*pi*((sd)^2))
end
######Rivalry with regular spiking
kee = .48#parse(Float64, ARGS[1])
kei = .93#parse(Float64, ARGS[2])
kie_L = .97#parse(Float64, ARGS[3])
kii = .9#parse(Float64, ARGS[4])

Aee = 54#parse(Float64, ARGS[5])
Aei = 104#parse(Float64, ARGS[6])
Aie_L = 180#parse(Float64, ARGS[7])
Aii= 150#parse(Float64, ARGS[8])

#####Rivalry with irregular spiking
# kee = .26#parse(Float64, ARGS[1])
# kei = .93#parse(Float64, ARGS[2])
# kie_L = .97#parse(Float64, ARGS[3])
# kii = .5#parse(Float64, ARGS[4])
#
# Aee = 84#parse(Float64, ARGS[5])
# Aei = 314#parse(Float64, ARGS[6])
# Aie_L = 1319#parse(Float64, ARGS[7])
# Aii= 689#parse(Float64, ARGS[8])

# kee = .18#parse(Float64, ARGS[1])
# kei = .93#parse(Float64, ARGS[2])
# kie_L = .97#parse(Float64, ARGS[3])
# kii = .5#parse(Float64, ARGS[4])
#
# Aee = 184#parse(Float64, ARGS[5])
# Aei = 514#parse(Float64, ARGS[6])
# Aie_L = 2319#parse(Float64, ARGS[7])
# Aii= 689#parse(Float64, ARGS[8])

kee = sd_2_k(kee)
kei = sd_2_k(kei)
kie_L = sd_2_k(kie_L)
kii = sd_2_k(kii)

s_strength = 2.2#parse(Float64, ARGS[9])
p = .34#parse(Float64, ARGS[10])
Ne = 3200#parse(Int64, ARGS[11])
# Ne = 1600
Ni = div(Ne, 4)
h = .1
vth = 20
tau_m = 20.
tau_ee = 2
tau_ei = 2
tau_ie = 2
tau_ii = 2
tau_ae = 700.
tau_ai = 1000.
# g_e = 0.
g_e = .001
g_i = .0
pee = p
pei = p
pie = p
pii = p

s1 = s_strength
s2 = s_strength

# Ne = 2000
# Ni = 500
# runtime = 500000
runtime = 20000 #ms
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 50/h
vbinszie = div(fbinsize, 2)
half_e = div(Ne, 2)

Aee /= p
Aei /= p
Aie_L /= p
Aii /= p
tic()
c=0
wee,wei,wie,wii = weights(Ne,Ni,kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, p,p,p,p);
# for g_e = .0001:0.005:.05
te,re,ti,ri,kill_flag, e_top, e_bot, adapt_top, adapt_bot = euler_lif(h,runtime,Ne,Ni,wee,wei,wie,wii,s1,s2, vth, tau_m, tau_ee, tau_ei, tau_ie, tau_ii, tau_ae, tau_ai, g_e, g_i);

ntd, nts = nt_diff_H(te, re, ntotal, half, 150/h)
s = ntd./nts #signal for dominances
flags, times = WLD_01(s, -.333, .333)

d = convert(Array{Float64}, diff(1500/10000. .* times))
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

println("CVD is: $(cvd)")
println("CVD_LP is: $(cvdlp)")

# t2, f2 = splice_reversions(flags, times)
# fw = find(f2 .== "win")
# fl = find(f2 .== "lose")
# tx = t2*250.
# d = convert(Array{Float64}, diff(250*t2))
# cvd = cv(d)
# a = []
# dt = diff(times*250.)
# dta = convert(Array{Float64}, dt)



# #
# #
# TN, BN = Neurons_tb_ns(re, half, 10, 100) #neurons in either pool who fired at least 10 spkes in 50 seconds
# top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times) #find win, lose, and draw times
# tbf, rbf = ligase(bot, bdom, te, re, BN) #bottom pool up states
# ttf, rtf = ligase(top, tdom, te, re, TN) #top pool up states
# tbdf, rbdf = ligase(top, tdom, te, re, BN) #bottom pool down states
# ttdf, rtdf = ligase(bot, bdom, te, re, TN) #top pool down states
# # countFT = count_train_intron(fbinsize, ttf, rtf, TN, length(TN), false)
# # countFB = count_train_intron(fbinsize, tbf, rbf, BN, length(BN), false)
# cbinsize = 1000/h
# #correlations
# cwTu = rand_pair_cor(cbinsize, ttf, rtf, TN, 1000)
# cwBu = rand_pair_cor(cbinsize, tbf, rbf, BN, 1000)
# cwBd = rand_pair_cor(cbinsize, ttdf, rtdf, TN, 1000)
# cwTd = rand_pair_cor(cbinsize, tbdf, rbdf, BN, 1000)

# CVSTu = CV_ISI_D(top, TN, te, re)
# CVSBu = CV_ISI_D(bot, BN, te, re)
# CVSTd = CV_ISI_D(bot, TN, te, re)
# CVSBd = CV_ISI_D(top, BN, te, re)

# subplot(221)
# plt[:hist](cwTu, 50)
# axvline(mean(cwTu), linestyle = "dashed", color = "g")
# gca()
# xlim(-1, 1)
# title("Up-State")
# ylabel("Pool 1")
# xticks([])
# subplot(222)
# plt[:hist](cwTd, 50)
# axvline(mean(cwTd), linestyle = "dashed", color = "g")
# gca()
# title("Down-State")
# xlim(-1, 1)
# xticks([])
# subplot(223)
# plt[:hist](cwBu, 50)
# axvline(mean(cwBu), linestyle = "dashed", color = "g")
# gca()
# ylabel("Pool 2")
# xlabel("Pairwise rsc")
# xlim(-1, 1)
# subplot(224)
# plt[:hist](cwBd, 50)
# axvline(mean(cwBd), linestyle = "dashed", color = "g")
# gca()
# xlabel("Pairwise rsc")
# xlim(-1, 1)


# plt[:hist](d, 50)

# subplot(221)
# title("Up-State")
# ax1 = gca()
# plt[:hist](CVSTu, 50)
# xlabel("CV(isi)")
# subplot(222)
# ax2 = gca()
# title("Down-State")
# plt[:hist](CVSTd, 50)
# xlabel("CV(isi)")
# subplot(223)
# ax3 = gca()
# plt[:hist](cwTu, 50)
# xlabel("Spike Count Correlation")
# subplot(224)
# ax4 = gca()
# plt[:hist](cwTd, 50)
# xlabel("Spike Count Correlation")




# TN, BN = Neurons_tb_ns(re, half, 10., 100)
#
# CVST = CV_ISI_ALLTIME(TN, te, re)
# CVSB = CV_ISI_ALLTIME(BN, te, re)
# cbinsize = 1000/h
# rpcT = rand_pair_cor(cbinsize, te, re, TN, 1000)
# rpcB = rand_pair_cor(cbinsize, te, re, BN, 1000)
#
# subplot(121)
# plt[:hist](rpcB, 50)
# axvline(mean(rpcB), linestyle = "dashed", color = "g")
# gca()
# title("Winner")
# subplot(122)
# plt[:hist](rpcT)
# axvline(mean(rpcT), linestyle = "dashed", color = "g")
# gca()
# title("Loser")
#


# AN = collect(union(Set(TN), Set(BN)))
# CVSA = CV_ISI_ALLTIME(AN, te, re)
# cbinsize = 50/h
# rpcA = rand_pair_cor(cbinsize, te, re, AN, 1000)
# plt[:hist](rpcA, 50)
# title("Pairwise Correlations in Normalization Regime")
# xlabel("Pariwise rsc")
# axvline(mean(rpcA), linestyle = "dashed", color = "g")

#
# subplot(121)
# ax2 = gca()
# plt[:hist](CVSA, 50)
# xlabel("CV(isi)")
# subplot(122)
# ax3 = gca()
# plt[:hist](rpcA, 50)
# xlabel("Spike Count Correlation")

# subplot(221)
# title("Winner")
# ax1 = gca()
# plt[:hist](CVSB, 50)
# xlabel("CV(isi)")
# subplot(222)
# ax2 = gca()
# title("Loser")
# plt[:hist](CVST, 50)
# xlabel("CV(isi)")
# subplot(223)
# ax3 = gca()
# plt[:hist](rpcB, 50)
# xlabel("Spike Count Correlation")
# subplot(224)
# ax4 = gca()
# plt[:hist](rpcT, 50)
# xlabel("Spike Count Correlation")




# FD = [length(find(re.==i)) for i = 1:Ne]
# if kill_flag == true
#   println("##RESULT -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
#   cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10
# else
# #ntd, nts = ntf(te, re, ntotal, half, 250)
# ntd, nts = nt_diff(te, re, ntotal, half, 25/h)
# s = ntd./nts #signal for dominances
# #flags, times = win_lose_draw(s)
# flags, times = WLD_01(s)
# top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times)
# rhz = round(.1*runtime/1000.)
# if (tdom > ntotal*.1) & (bdom > ntotal*.1)
#   #rivalry
#   println("RIVALING")
#   TN, BN = Neurons_tb_ns(re, half, rhz, 100)
#   if ((TN == -5) | (BN == -5))
#     println("##RESULT -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
#     cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10
#   else
#   t2, f2 = splice_reversions(flags, times)
#   fw = find(f2 .== "win")
#   fl = find(f2 .== "lose")
#   tx = t2*250.
#   d = convert(Array{Float64}, diff(250*t2))
#   cvd = cv(d)
#   MDT, MDB = (tdom*1./length(bot)), (bdom*1./length(top))
#   tbf, rbf = ligase(bot, bdom, te, re, BN)
#   ttf, rtf = ligase(top, tdom, te, re, TN)
#   tbdf, rbdf = ligase(top, tdom, te, re, BN)
#   ttdf, rtdf = ligase(bot, bdom, te, re, TN)
#   #countCT = count_train_intron(cbinsize, ttf, rtf, TN, 50, true)
#   #countCB = count_train_intron(cbinsize, tbf, rbf, BN, 50, true)
#   countFT = count_train_intron(fbinsize, ttf, rtf, TN, length(TN), false)
#   countFB = count_train_intron(fbinsize, tbf, rbf, BN, length(BN), false)
#   #cwT, cwB = correlate_within(countCT, -5), correlate_within(countCB, -5)
#   cwT = rand_pair_cor(cbinsize, ttf, rtf, TN, 1000)
#   cwB= rand_pair_cor(cbinsize, tbf, rbf, BN, 1000)
#   MFT, medFT, stdFT = fano_train(countFT, -5)
#   MFB, medFB, stdFB = fano_train(countFB, -5)
#   mcvt, medcvt, stdcvt = CV_ISI(top, TN, te, re)
#   mcvb, medcvb, stdcvb = CV_ISI(bot, BN, te, re)
#   mcvtd, medcvtd, stdcvtd = CV_ISI(top, BN, tbdf, rbdf)
#   mcvbd, medcvbd, stdcvbd = CV_ISI(bot, TN, ttdf, rtdf)
#
#   mtu, vtu, mbd, vbd, baltu, balbd, e_top_up, e_bot_dn = input_analysis(top, e_top, e_bot)
#   mtd, vtd, mbu, vbu, baltd, balbu, e_top_dn, e_bot_up = input_analysis(bot, e_top, e_bot)
#   e_tu, e_bd = ligate_inputs(top, e_top, e_bot)
#   e_bu, e_td = ligate_inputs(bot, e_top, e_bot)
#   mtu, stdtu, sktu, kutu = moments(e_tu)
#   mbd, stdbd, skbd, kubd = moments(e_bd)
#   mtd, stdtd, sktd, kutd = moments(e_td)
#   mbu, stdbu, skbu, kubu = moments(e_bu)
#   println("##RESULT 1, $(MDT), $(MDB), $(cvd), $(cwT), $(cwB), $(mtu), $(stdtu), $(sktu), $(kutu), $(mbd), $(stdbd), $(skbd), $(kubd), $(mtd), $(stdtd), $(sktd), $(kutd), $(mbu), $(stdbu), $(skbu), $(kubu), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
# end
#
# elseif tnmz > 160000
#   println("NORMALIZATION")
#   #normalization
#   TN, BN = Neurons_tb_ns(re, half, rhz, 100)
#   if ((TN == -5) | (BN == -5))
#     println("##RESULT -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
#     cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8
#   else
#   tbf, rbf = ligase(nmz, tnmz, te, re, BN)
#   ttf, rtf = ligase(nmz, tnmz, te, re, TN)
#   #countCT = count_train_intron(cbinsize, ttf, rtf, TN, 50, true)
#   #countCB = count_train_intron(cbinsize, tbf, rbf, BN, 50, true)
#   countFT = count_train_intron(fbinsize, ttf, rtf, TN, length(TN), false)
#   countFB = count_train_intron(fbinsize, tbf, rbf, BN, length(BN), false)
#   cwT = rand_pair_cor(cbinsize, ttf, rtf, TN, 1000)
#   cwB = rand_pair_cor(cbinsize, tbf, rbf, BN, 1000)
#   #cwT, cwB = correlate_within(countCT, -5), correlate_within(countCB, -5)
#   MFT, medFT, stdFT = fano_train(countFT, -5)
#   MFB, medFB, stdFB = fano_train(countFB, -5)
#   mcvt, medcvt, stdcvt = CV_ISI(nmz, TN, te, re)
#   mcvb, medcvb, stdcvb = CV_ISI(nmz, BN, te, re)
#   ts, bs = ligate_inputs(nmz, e_top, e_bot)
#   mtu, stdtu, sktu, kutu = moments(ts)
#   mbd, stdbd, skbd, kubd = moments(bs)
#   #mtu, vtu, mbd, vbd, baltu, balbd = input_analysis(nmz, e_top, e_bot)
#   #mtd, stdtd, mbu, vbu, baltd, balbu = -7, -7, -7, -7, -7, -7
#   mtd, stdtd, sktd, kutd, mbu, stdbu, skbu, kubu = -7, -7, -7, 7, -7, -7, -7, -7
#   println("##RESULT 0, -7, -7, $(tnmz), $(cwT), $(cwB), $(mtu), $(stdtu), $(sktu), $(kutu), $(mbd), $(stdbd), $(skbd), $(kubd), $(mtd), $(stdtd), $(sktd), $(kutd), $(mbu), $(stdbu), $(skbu), $(kubu), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
# end
#
#   #mtd, vtd, mbu, vbu = input_analysis(bot, e_top, e_bot)
# elseif (tdom > .8*ntotal) | (bdom > .8*ntotal)
#   println("WINNER TAKE ALL")
#   TN, BN = Neurons_tb_ns(re, half, rhz, 100)
#   n1, n2 = length(TN), length(BN)
#   println("here's $(n1), and $(n2)")
#   if ((TN == -5) & (BN == -5))
#     println("##RESULT -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
#     cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8
#   else
#   if n1 > n2
#     Neurons = TN
#     s1 = e_top
#     s2 = e_bot
#     exon = top
#     dom = tdom
#   else
#     Neurons = BN
#     s1 = e_bot
#     s2 = e_top
#     exon = bot
#     dom = bdom
#   end
#   tbf, rbf = ligase(exon, dom, te, re, Neurons)
#   #countCB = count_train_intron(cbinsize, tbf, rbf, Neurons, 10, true)
#   countFB = count_train_intron(fbinsize, tbf, rbf, Neurons, length(Neurons), false)
#   cwB = rand_pair_cor(cbinsize, tbf, rbf, Neurons, 1000)
#   #cwB = correlate_within(countCB, -5)
#   MFB, medFB, stdFB = fano_train(countFB, -5)
#   mcvb, medcvb, stdcvb = CV_ISI(exon, Neurons, te, re)
#   #ts, bs = ligate_inputs(exon, s1, s2)
#   mtd, stdtd, sktd, kutd = moments(s2)
#   mbu, stdbu, skbu, kubu = moments(s1)
#   mcvt, medcvt, stdcvt, MFT, medFT, stdFT = -7, -7, -7, -7, -7, -7
#   mtu, stdtu, sktu, kutu, mbd, stdbd, skbd, kubd = -7, -7, -7, -7, -7, -7, -7, 7
#   #mtu, vtu, mbd, vbd, baltu, balbd = input_analysis(exon, e_top, e_bot)
#   #mtd, vtd, mbu, vbu, baltd, balbu = -7, -7, -7, -7, -7, -7
#   cwT = -7
#   println("##RESULT 2, -7, -7, $(tnmz), -7, $(cwB), $(mtu), $(stdtu), $(sktu), $(kutu), $(mbd), $(stdbd), $(skbd), $(kubd), $(mtd), $(stdtd), $(sktd), $(kutd), $(mbu), $(stdbu), $(skbu), $(kubu), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
# end
# else
#   #some kind of error
#   println("##RESULT -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9,-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
# # end
# end
# end

# toc()
# half = div(Ne, 2)
# quarter = div(half, 2)
# rt = half + quarter
# rb = quarter
#
# figure(1)
# plot(te, re, "g.", ms = 1.)
# title("cwT=$(cwT), cwB=$(cwB), MFT=$(MFT), MFB=$(MFB), cvt=$(mcvt), cvb=$(mcvb)", fontsize = 6)
# xlabel("kee=$(kee), kei=$(kei), kie=$(kie_L), kii=$(kii), Aee=$(Aee*p), Aei=$(Aei*p), Aie=$(Aie_L*p), Aii=$(Aii*p), S=$(s_strength), p=$(p), g=$(g_e)", fontsize = 6)
# # figure(2)
# plot(e_top, "r")
# x=te[find(re.==rt)]
# plot(x, fill(3, length(x)), "g.", ms=20)
# plot(t2*250., fill(5, length(t2)), "b.", ms = 50)
# title("Asynchronous-Irregular Rivalry Inputs: Single Representative", fontsize = 48)
# xticks(fontsize = 32)
# yticks(fontsize = 32)
# xlabel("Time (dt = 0.1ms); mean=$(mean(e_top[1000,:])), std=$(std(e_top[1000,:]))S", fontsize = 42)
# ylabel("Input (Drive, Adaptation, & Recurent Synapses)", fontsize = 42)
# # plot(e_bot, "b")
# # title("vtu=$(vtu), mtu=$(mtu); vtd=$(vtd), mtd=$(mtd)", fontsize = 24)
# # xlabel("vbu=$(vbu), mbu=$(mbu); vbd=$(vbd), mbd=$(mbd)", fontsize = 24)
# # ylabel("kee=$(kee), kei=$(kei), kie=$(kie_L), kii=$(kii), Aee=$(Aee*p), Aei=$(Aei*p), Aie=$(Aie_L*p), Aii=$(Aii*p), S=$(s_strength), p=$(p), g=$(g_e)", fontsize = 14)
# figure(3)
# plt[:hist](diff(t2)*250, 15)
# println(cv(diff(t2)), "is the CVD")
#
# # xinp = []
# # c = 1
# # while c < length(e_top_up)-20
# #   x = sum(e_top_up[c:c+19])
# #   push!(xinp, x)
# #   c+=20
# # end
#
# xin = convert(Array{Float64}, xinp)
# #xin = zeros(div(length(e_top_up),2))
# figure(4)
# title("Histogram of Input Net Amplitude In 2ms bins, Rivalry Network", fontsize = 48)
# #plt[:hist](xin, 1000)
# plt[:hist](e_top_up, 1000)
# xticks(fontsize = 32)
# yticks(fontsize = 32)
# xlabel("Summed Input in 2ms bin; mu = $(mean(xin)), sigma=$(std(xin))", fontsize = 42)
# ylabel("Frequency", fontsize = 42)
# etu = zeros(length(P1), length(e_tu));
# ebd = zeros(length(P1), length(e_tu));
# ebu = zeros(length(P1), length(e_bu));
# etd = zeros(length(P1), length(e_bu));
# for i in eachindex(P1)
#   etu[i,1:end], ebd[i,1:end] = ligate_inputs(top, input[P2[i],:][:], input[P1[i],:][:])
#   ebu[i,1:end], etd[i,1:end] = ligate_inputs(bot, input[P1[i],:][:], input[P2[i],:][:])
# end

# etwm = zeros(length(P2))
# ebwm = zeros(length(P2))
# etlm = zeros(length(P2))
# eblm = zeros(length(P2))
# etws = zeros(length(P2))
# ebws = zeros(length(P2))
# etls = zeros(length(P2))
# ebls = zeros(length(P2))
# for i in eachindex(P2)
#   etwm[i] = mean(etu[i,:])
#   etlm[i] = mean(etd[i,:])
#   ebwm[i] = mean(ebu[i,:])
#   eblm[i] = mean(ebd[i,:])
#   etws[i] = std(etu[i,:])
#   etls[i] = std(etd[i,:])
#   ebws[i] = std(ebu[i,:])
#   ebls[i] = std(ebd[i,:])
# end
# write_raster(etwm, "riv_top_win_means.txt")
# write_raster(ebwm, "riv_bot_win_means.txt")
# write_raster(etlm, "riv_top_los_means.txt")
# write_raster(eblm, "riv_bot_los_means.txt")
# write_raster(etws, "riv_top_win_stds.txt")
# write_raster(ebws, "riv_bot_win_stds.txt")
# write_raster(etls, "riv_top_los_stds.txt")
# write_raster(ebls, "riv_bot_los_stds.txt")
#2624


function nt_diff_H(t, r, ntotal, half, netd_binsize)

  netd_bins = collect(1:netd_binsize:ntotal)
  ntd = zeros(length(netd_bins)-1)
  nts = zeros(length(netd_bins)-1)

  for j = 2:length(netd_bins)
    tf = find(netd_bins[j-1] .<= t .< netd_bins[j])
    T = sum(r[tf] .> half)
    B = sum(r[tf] .<= half)
    NTD = T - B #################Differences in Spikes
    ntd[j-1] = NTD
    nts[j-1] = length(tf)
    #println(T+B==length(tf))
  end

  return ntd, nts
end

function WLD_01(s, tl, th)
  if maximum(s) < tl
    return ["lose", "end"], [1, length(s)]
  elseif minimum(s) >= th
    return ["win", "end"], [1, length(s)]
  elseif (tl <= maximum(s) <= th) & (tl <= minimum(s) <= th)
    return ["draw", "end"], [1, length(s)]
  end
  times = []
  flags = []
  if s[1] >= th
    flag = "win"
  elseif tl <= s[1] <= th
    flag = "draw"
  elseif s[1] < tl
    flag = "lose"
  end

  push!(times, 1)
  push!(flags, flag)

  s2 = s[2:end]
  for i in eachindex(s2)
    f1 = comp_01(s2[i], tl, th)
    if f1 != flag
      flag = f1
      push!(times, i+1)
      push!(flags, flag)
    end
  end
  push!(times, length(s))
  push!(flags, "end")
  return flags, times
end

function comp_01(x, tl, th)
  if x > th
    return "win"
  elseif tl <= x <= th
    return "draw"
  elseif x < tl
    return "lose"
  else
    return "weird"
  end
end
