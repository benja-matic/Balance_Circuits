include("Euler_W.jl")
include("Analyze.jl")

srand(4321)
function sd_2_k(sd)
  return 1./(2*pi*((sd)^2))
end

# kee = .26 #WTA
# kee = .2
# kee = .96 #Normalization
# kei = .93
# kie = .97
# kii = .5

kee = .95
kei = .95
kie = .95
kii = .95

Aee = 84
# Aee = 210.
# Aee = 120.
Aei = 314
Aie = 3500.
# Aei = 564.
# Aie = 5000.
# Aie = 1319
# Aii= 689
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
tau_a = 575.
# tau_a = 425.
g_a = 0
# g_a = .001
# g_a = .002
# g_a = .0075
# kee = sd_2_k(kee)
# kei = sd_2_k(kei)
# kie = sd_2_k(kie)
# kii = sd_2_k(kii)

s1 = s_strength
s2 = s_strength

min_e_neurons = 20
min_i_neurons = 50
runtime = 10000 #ms
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
# Angle = [0, 100, 200, 300, 400]
# for angle in Angle
  # @time t, r, Input, adapt, Exc, Inh, Vlt = euler_lif_CSR(h, runtime, Ne, W, CSR, s1, s2, vth, tau_m, tau_s, tau_a, g_a, angle)
@time t, r, Input = euler_lif_CSR(h, runtime, Ne, W, CSR, s1, s2, w2, input_width, vth, tau_m, tau_s, tau_a, g_a, ang)

ex = find(r .<= Ne)
te_pt = t[ex]
re_pt = r[ex]

ix = find(r .> Ne)
ti_pt = t[ix]
ri_pt = r[ix]
half_e = div(Ne, 2)

wta_ness, bias = score_analysis(re_pt, Ne)
top_e_neurons, bot_e_neurons = Neurons_tb_ns(re_pt, half_e, 10, min_e_neurons)
I_Neurons = Neuron_finder(ri_pt, 10, min_i_neurons)

CV_TOP = CV_ISI_ALLTIME(top_e_neurons, te_pt, re_pt)
CV_BOT = CV_ISI_ALLTIME(bot_e_neurons, te_pt, re_pt)
CV_INH = CV_ISI_ALLTIME(I_Neurons, ti_pt, ri_pt)

# E_R_top = [length(find(re_pt .== i))/rt for i in top_e_neurons]
# E_R_bot = [length(find(re_pt .== i))/rt for i in bot_e_neurons]
E_R_top = [length(find(re_pt .== i))/rt for i=half_e+1:Ne]
E_R_bot = [length(find(re_pt .== i))/rt for i=1:half_e]
I_R = [length(find(ri_pt .== i))/rt for i=Ne+1:N]

push!(IR, mean(I_R))
push!(ERT, mean(E_R_top))
push!(ERB, mean(E_R_bot))
push!(CVET, mean(CV_TOP))
push!(CVEB, mean(CV_BOT))
push!(CVIT, mean(CV_INH))
# IR = []
# ERT = []
# ERB = []
# CVET = []
# CVEB = []
# CVIT = []
# end
# wta_ness, bias, E_Rates, I_Rates, CV_TOP, CV_BOT, CV_INH, RSC_TOP, RSC_BOT, RSC_INH, FF_TOP, FF_BOT, FF_INH, te_pt, re_pt = WTAN_no_input(t, r, end_trans, Ne, Ni, half_e, min_e_neurons, min_i_neurons, fbinsize, cbinsize, FPe)
Aie = [500, 1000, 1500, 2000, 2500, 3000, 3500]
# E_diff = ERT .- ERB
# EDA = [abs(i) for i in E_diff]
subplot(211)
ylabel("Rate (Hz)")
plot(Aie, IR, ".", label = "Inhibitory")
plot(Aie, ERT, ".", label = "Excitatory Pool 1")
plot(Aie, ERB, ".", label = "Excitatory Pool 2")
# plot(Aie, EDA, ".", label = "Difference of Excitatory Pools")
legend()
# plot(Aie, ERT, ".", label = "Excitatory 2")
subplot(212)
plot(Aie, CVIT, ".", label = "Inhibitory")
plot(Aie, CVEB, ".", label = "Excitatory 1")
plot(Aie, CVET, ".", label = "Excitatory 2")
legend()
ylabel("CV_isi")
xlabel("Aie")



half_pool = div(Ne, w2)
pool_size = half_pool*2
incoming = 2*collect(0:pool_size-1)*pi/pool_size
fe1 = s1*pool_size .* circshift(von_mises_dist(incoming, input_width, 0, pool_size), half_pool)
fe2 = s2*pool_size .* circshift(von_mises_dist(incoming, input_width, 0, pool_size), half_pool)

drive = zeros(N)
P1s = round(Int,Ne/4-(half_pool -1)) + ang
P1e = round(Int,Ne/4+half_pool) + ang
P2s = round(Int,3*Ne/4-(half_pool-1)) - ang
P2e = round(Int,3*Ne/4+half_pool) - ang

drive[P1s:P1e] .+= fe1*h
drive[P2s:P2e] .+= fe2*h

maxdrive = maximum(fe1)
for i in eachindex(drive)
    if drive[i] > maxdrive
        drive[i] = maxdrive
    end
end

# subplot(121)
# plot(te, re, "g.", ms = 1.)
# subplot(122)
# plot(drive, collect(1:length(drive)))


# figure(2)
# plot(t, r, "g.", ms = 1.)
# #
#
# ntd, nts = nt_diff_angle(te, re, ntotal, P2s, P2e, netd_binsize)
# s = ntd ./nts
# flags, times = WLD_01(s, -.333, .333)
# top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times, netd_binsize) #find win, lose, and draw times
# MDT = tdom/length(top)
# MDB = bdom/length(bot)
# MDN = tnmz/length(nmz)
# MDD = ntotal/length(times)
# println("$(MDT), $(MDB), $(length(times))")
# push!(alts, length(times))
# push!(MDTs, MDT)
# push!(MDBs, MDB)
# end
# Angle = [0, 100, 200, 300, 400, 500, 600]
# degs = 90 - (Angle*(180/Ne))

#
# plot(t, r, "g.", ms = 1.)
# plot(times .* netd_binsize, fill(P2s, length(times)), "r.")


# MD = ntotal ./ alts
# MD ./= 10000.
# subplot(211)
# plot(degs, alts, ".", ms = 20.)
# ylabel("# Alternations")
# subplot(212)
# plot(degs, MD, ".", ms = 20.)
# xlabel("Angle (degrees)")
# ylabel("1/# Alternations")
#
# ntd, nts = nt_diff_H(te, re, ntotal, half_e, netd_binsize)
# s = ntd ./ nts #signal for dominances
# flags, times = WLD_01(s, -.333, .333)
#
# d = convert(Array{Float64}, diff(netd_binsize/(1000./h) .* times))
# cvd = cv(d) ###Raw estimate of CVD, likely to include very rapid switches which should really be smoothed out
#
# LP = .3
#
# dx = []
# for i in d
#     if i > LP
#         push!(dx, i)
#     end
# end
# dx = convert(Array{Float64}, dx)
# cvdlp = cv(dx) ###high-pass filter measure of CVD

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
