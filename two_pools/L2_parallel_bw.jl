#println("GOT THIS FAR")
#include("measure_balance.jl")
include("euler_b_2_rivalry.jl")
srand(1234)
function sd_2_k(sd)
  return 1./(2*pi*((sd)^2))
end

kee = .26#parse(Float64, ARGS[1])
kei = .93#parse(Float64, ARGS[2])
kie_L = .97#parse(Float64, ARGS[3])
kii = .5#parse(Float64, ARGS[4])

Aee = 84#parse(Float64, ARGS[5])
Aei = 314#parse(Float64, ARGS[6])
Aie_L = 1319#parse(Float64, ARGS[7])
Aii= 689#parse(Float64, ARGS[8])

kee = sd_2_k(kee)
kei = sd_2_k(kei)
kie_L = sd_2_k(kie_L)
kii = sd_2_k(kii)

s_strength = 3.08*1.02
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
tau_ae = 750.
tau_ai = 1000.
g_e = .0015
g_i = .0
pee = p
pei = p
pie = p
pii = p

#0.8338799042222288, 0.17736295260301269, 0.2330863885042438, 0.23766272735583688, 132.632740955, 106.477147728, 851.027483616, 376.44001464, 4.66530540123, 0.44446060668, 3944.464966, 0.0095363078743, 0.0

# kee = sd_2_k(0.37634223780717)
# #kee = sd_2_k(.8)
# kei = sd_2_k(0.74471527543)
# kie_L = sd_2_k(0.996261417223)
# kii = sd_2_k(0.5336451638526)
#s_strength = 1.8
# kee = sd_2_k(.37)
# kei = sd_2_k(1)
# kie_L = sd_2_k(1)
# kii = sd_2_k(.5)

# kee = sd_2_k(.403678)
# kei = sd_2_k(.94728)
# kie_L = sd_2_k(.82632)
# kii = sd_2_k(.81833)

# kee = sd_2_k(.2341)
# kei = sd_2_k(.4973)
# kie_L = sd_2_k(.555)
# kii = sd_2_k(.2189)
# Aee = 100.#94.#80.19552235
# Aei = 211#314.75109066
# Aie_L = 808.52219365
# Aii = 648.#130.916487007

# Aee = 132.6#94.#80.19552235
# Aei = 106.4#314.75109066
# Aie_L = 851.052219365
# Aii = 576.4#130.916487007
#
# # p = .2
# s_strength = 2.
# #p = 0.165403843441
# p = .12
# s_strength = 2.144452373606
s1 = 3.08
s2 = s_strength

# Ne = 2000
# Ni = 500
runtime = 100000 #ms
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
vbinszie = div(fbinsize, 2)
half = div(Ne, 2)

Aee /= p
Aei /= p
Aie_L /= p
Aii /= p
tic()
c=0
wee,wei,wie,wii = weights(Ne,Ni,kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, p,p,p,p);
# for g_e = .0001:0.005:.05
te,re,ti,ri,kill_flag, e_top, e_bot = euler_lif(h,runtime,Ne,Ni,wee,wei,wie,wii,s1,s2, vth, tau_m, tau_ee, tau_ei, tau_ie, tau_ii, tau_ae, tau_ai, g_e, g_i);
FD = [length(find(re.==i)) for i = 1:Ne]
if kill_flag == true
  println("##RESULT -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
  cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10
else
#ntd, nts = ntf(te, re, ntotal, half, 250)
ntd, nts = nt_diff(te, re, ntotal, half, 25/h)
s = ntd./nts #signal for dominances
#flags, times = win_lose_draw(s)
flags, times = WLD_01(s)
top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times)
rhz = round(.1*runtime/1000.)
if (tdom > ntotal*.1) & (bdom > ntotal*.1)
  #rivalry
  println("RIVALING")
  TN, BN = Neurons_tb_ns(re, half, rhz, 100)
  if ((TN == -5) | (BN == -5))
    println("##RESULT -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
    cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10
  else
  t2, f2 = splice_reversions(flags, times)
  d = convert(Array{Float64}, diff(250*t2))
  cvd = cv(d)
  MDT, MDB = (tdom*1./length(bot)), (bdom*1./length(top))
  tbf, rbf = ligase(bot, bdom, te, re, BN)
  ttf, rtf = ligase(top, tdom, te, re, TN)
  tbdf, rbdf = ligase(top, tdom, te, re, BN)
  ttdf, rtdf = ligase(bot, bdom, te, re, TN)
  #countCT = count_train_intron(cbinsize, ttf, rtf, TN, 50, true)
  #countCB = count_train_intron(cbinsize, tbf, rbf, BN, 50, true)
  countFT = count_train_intron(fbinsize, ttf, rtf, TN, length(TN), false)
  countFB = count_train_intron(fbinsize, tbf, rbf, BN, length(BN), false)
  #cwT, cwB = correlate_within(countCT, -5), correlate_within(countCB, -5)
  cwT = rand_pair_cor(cbinsize, ttf, rtf, TN, 1000)
  cwB= rand_pair_cor(cbinsize, tbf, rbf, BN, 1000)
  MFT, medFT, stdFT = fano_train(countFT, -5)
  MFB, medFB, stdFB = fano_train(countFB, -5)
  mcvt, medcvt, stdcvt = CV_ISI(top, TN, te, re)
  mcvb, medcvb, stdcvb = CV_ISI(bot, BN, te, re)
  mtu, vtu, mbd, vbd, baltu, balbd, e_top_up, e_bot_dn = input_analysis(top, e_top, e_bot)
  mtd, vtd, mbu, vbu, baltd, balbu, e_top_dn, e_bot_up = input_analysis(bot, e_top, e_bot)
  e_tu, e_bd = ligate_inputs(top, e_top, e_bot)
  e_bu, e_td = ligate_inputs(bot, e_top, e_bot)
  mtu, stdtu, sktu, kutu = moments(e_tu)
  mbd, stdbd, skbd, kubd = moments(e_bd)
  mtd, stdtd, sktd, kutd = moments(e_td)
  mbu, stdbu, skbu, kubu = moments(e_bu)
  println("##RESULT 1, $(MDT), $(MDB), $(cvd), $(cwT), $(cwB), $(mtu), $(stdtu), $(sktu), $(kutu), $(mbd), $(stdbd), $(skbd), $(kubd), $(mtd), $(stdtd), $(sktd), $(kutd), $(mbu), $(stdbu), $(skbu), $(kubu), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
end

elseif tnmz > 160000
  println("NORMALIZATION")
  #normalization
  TN, BN = Neurons_tb_ns(re, half, rhz, 100)
  if ((TN == -5) | (BN == -5))
    println("##RESULT -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
    cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8
  else
  tbf, rbf = ligase(nmz, tnmz, te, re, BN)
  ttf, rtf = ligase(nmz, tnmz, te, re, TN)
  #countCT = count_train_intron(cbinsize, ttf, rtf, TN, 50, true)
  #countCB = count_train_intron(cbinsize, tbf, rbf, BN, 50, true)
  countFT = count_train_intron(fbinsize, ttf, rtf, TN, length(TN), false)
  countFB = count_train_intron(fbinsize, tbf, rbf, BN, length(BN), false)
  cwT = rand_pair_cor(cbinsize, ttf, rtf, TN, 1000)
  cwB = rand_pair_cor(cbinsize, tbf, rbf, BN, 1000)
  #cwT, cwB = correlate_within(countCT, -5), correlate_within(countCB, -5)
  MFT, medFT, stdFT = fano_train(countFT, -5)
  MFB, medFB, stdFB = fano_train(countFB, -5)
  mcvt, medcvt, stdcvt = CV_ISI(nmz, TN, te, re)
  mcvb, medcvb, stdcvb = CV_ISI(nmz, BN, te, re)
  ts, bs = ligate_inputs(nmz, e_top, e_bot)
  mtu, stdtu, sktu, kutu = moments(ts)
  mbd, stdbd, skbd, kubd = moments(bs)
  #mtu, vtu, mbd, vbd, baltu, balbd = input_analysis(nmz, e_top, e_bot)
  #mtd, stdtd, mbu, vbu, baltd, balbu = -7, -7, -7, -7, -7, -7
  mtd, stdtd, sktd, kutd, mbu, stdbu, skbu, kubu = -7, -7, -7, 7, -7, -7, -7, -7
  println("##RESULT 0, -7, -7, $(tnmz), $(cwT), $(cwB), $(mtu), $(stdtu), $(sktu), $(kutu), $(mbd), $(stdbd), $(skbd), $(kubd), $(mtd), $(stdtd), $(sktd), $(kutd), $(mbu), $(stdbu), $(skbu), $(kubu), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
end

  #mtd, vtd, mbu, vbu = input_analysis(bot, e_top, e_bot)
elseif (tdom > .8*ntotal) | (bdom > .8*ntotal)
  println("WINNER TAKE ALL")
  TN, BN = Neurons_tb_ns(re, half, rhz, 100)
  n1, n2 = length(TN), length(BN)
  println("here's $(n1), and $(n2)")
  if ((TN == -5) & (BN == -5))
    println("##RESULT -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
    cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8
  else
  if n1 > n2
    Neurons = TN
    s1 = e_top
    s2 = e_bot
    exon = top
    dom = tdom
  else
    Neurons = BN
    s1 = e_bot
    s2 = e_top
    exon = bot
    dom = bdom
  end
  tbf, rbf = ligase(exon, dom, te, re, Neurons)
  #countCB = count_train_intron(cbinsize, tbf, rbf, Neurons, 10, true)
  countFB = count_train_intron(fbinsize, tbf, rbf, Neurons, length(Neurons), false)
  cwB = rand_pair_cor(cbinsize, tbf, rbf, Neurons, 1000)
  #cwB = correlate_within(countCB, -5)
  MFB, medFB, stdFB = fano_train(countFB, -5)
  mcvb, medcvb, stdcvb = CV_ISI(exon, Neurons, te, re)
  #ts, bs = ligate_inputs(exon, s1, s2)
  mtd, stdtd, sktd, kutd = moments(s2)
  mbu, stdbu, skbu, kubu = moments(s1)
  mcvt, medcvt, stdcvt, MFT, medFT, stdFT = -7, -7, -7, -7, -7, -7
  mtu, stdtu, sktu, kutu, mbd, stdbd, skbd, kubd = -7, -7, -7, -7, -7, -7, -7, 7
  #mtu, vtu, mbd, vbd, baltu, balbd = input_analysis(exon, e_top, e_bot)
  #mtd, vtd, mbu, vbu, baltd, balbu = -7, -7, -7, -7, -7, -7
  cwT = -7
  println("##RESULT 2, -7, -7, $(tnmz), -7, $(cwB), $(mtu), $(stdtu), $(sktu), $(kutu), $(mbd), $(stdbd), $(skbd), $(kubd), $(mtd), $(stdtd), $(sktd), $(kutd), $(mbu), $(stdbu), $(skbu), $(kubu), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
end
else
  #some kind of error
  println("##RESULT -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9,-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(g_e)")
# end
end
end

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
