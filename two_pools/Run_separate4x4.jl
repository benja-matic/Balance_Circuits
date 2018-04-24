include("4x4_separate_synapses.jl")
include("Analyze.jl")

srand(4321)


wta = []
ERS1 = []
ERS2 = []
IRS1 = []
IRS2 = []
ERT1 = []
ERT2 = []
IRT1 = []
IRT2 = []
b1L = []
b2L = []
c1L = []
c2L = []
d1L = []
d2L = []
WEE1 = []
WEI1 = []
WIE1 = []
WIEL1 = []
WII1 = []
WEE2 = []
WEI2 = []
WIE2 = []
WIEL2 = []
WII2 = []
REX1 = []
RIX1 = []
REX2 = []
RIX2 = []

SL1 = []#mean long range input coming into pool 1
SL2 = []

E1mn = []
I1mn = []
E2mn = []
I2mn = []

# ECV1 = []
# ECV2 = []
# ICV1 = []
# ICV2 = []
# means = []

# RET = []
# RIT = []
# RES = []
# RIS = []
# ECV1 = []
# ECV2 = []
# ICV1 = []
# ICV2 = []

# Aie_NL = [200, 300, 400, 500, 600, 700, 800, 900, 1000]
# for i in Aie_NL
# Aies = [10, 20, 30, 40, 60, 65, 70, 75]
# for i in Aies

# Ss = [1.5, 2., 2.5, 3., 3.5, 4.]
# Ss = [-.5, -.25, 0., 0.25, 0.5]
# # Ss = [2.8, 2.9, 3., 3.1, 3.2]
# for i in Ss
# Aie_NLS = linspace(0,80,10)
# Aie_NLS = [0, 10, 20, 30, 40, 50, 60, 65, 70, 75, 80]
# Aie_NLS = [71.5, 71.75, 72., 72.25, 72.5, 72.75, 73, 73.25, 73.5, 73.75, 74, 74.25, 74.5, 74.75, 75]
# for i in Aie_NLS
Aee = 5.
Aei = 40.
Aie = 70.
Aie_NL = 85.
Aii = 20.

#Aie = [200, 400, 600, 800, 1000]

s_strength = 3.08
fe = s_strength
fe2 = fe + .01
fi = 0.
fi2 = fi

# p = .2
k = 800

N = 16000
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
tau_a = 575.
g_a = 0

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

@time te, re, ti, ri, SEE, SEI, SIE, SIEL, SII = euler_lif_CSR_4x4_s(h, runtime, N, IFRAC, W, CSR, fe, fi, fe2, fi2, vth, tau_m, tau_s, tau_a, g_a)

wta_ness, bias = score_analysis(re, Ne2)
top_e_neurons, bot_e_neurons = Neurons_tb_ns(re, NeL, 10, min_e_neurons)
top_i_neurons, bot_i_neurons = Neurons_tb_ns(ri, NiL, 10, min_i_neurons)

CV_ETOP = CV_ISI_ALLTIME(top_e_neurons, te, re)
CV_EBOT = CV_ISI_ALLTIME(bot_e_neurons, te, re)
CV_ITOP = CV_ISI_ALLTIME(top_i_neurons, ti, ri)
CV_IBOT = CV_ISI_ALLTIME(bot_i_neurons, ti, ri)

E_R_top = [length(find(re .== i))/rt for i=NeL+1:Ne2]
E_R_bot = [length(find(re .== i))/rt for i=1:NeL]
I_R_top = [length(find(ri .== i))/rt for i=NiL+1:Ni2]
I_R_bot = [length(find(ri .== i))/rt for i=1:NiL]

# NSET = length(find(re .> NeL))
# NSEB = length(find(re .<= NeL))
# NSIT = length(find(ri .> NiL))
# NSIB = length(find(ri .<= NiL))

#leaving it in Hz, passing it into sim_2_theory, and the sim_2_theory output into theory_rates will just leave theory_rates output in Hz
MER1 = mean(E_R_bot)#*(1/1000.)
MER2 = mean(E_R_top)#*(1/1000.)
MIR1 = mean(I_R_bot)#*(1/1000.)
MIR2 = mean(I_R_top)#*(1/1000.)

println("##RESULT $(wta_ness), $(bias), $(mean(E_R_bot)), $(mean(E_R_top)), $(mean(I_R_top)), $(mean(I_R_bot)), $(mean(CV_ETOP)), $(mean(CV_EBOT)), $(mean(CV_ITOP)), $(mean(CV_IBOT))")

WEE_1, WEE_2, WIE_1, WIE_2, WIEL_1, WIEL_2, WEI_1, WEI_2, WII_1, WII_2, FE, FI = sim_2_theory(SEE, SEI, SIE, SIEL, SII, fe, fi, 1., MER1, MER2, MIR1, MIR2, 100)

push!(WEE1, WEE_1)
push!(WEI1, WEI_1)
push!(WIE1, WIE_1)
push!(WIEL1, WIEL_1)
push!(WII1, WII_1)

push!(WEE2, WEE_2)
push!(WEI2, WEI_2)
push!(WIE2, WIE_2)
push!(WIEL2, WIEL_2)
push!(WII2, WII_2)

println("##PARAMETERS $(WEE_1), $(WEE_2), $(WIE_1), $(WIE_2), $(WIEL_1), $(WIEL_2), $(WEI_1), $(WEI_2), $(WII_1), $(WII_2), $(FE), $(FI)")

RE1_THEORY, RI1_THEORY, RE2_THEORY, RI2_THEORY = theory_rates(abs(WEE_1), abs(WEE_2), abs(WIE_1), abs(WIE_2), abs(WIEL_1), abs(WIEL_2), abs(WEI_1), abs(WEI_2), abs(WII_1), abs(WII_2), FE, FI)

# RE1_THEORYc, RI1_THEORYc, RE2_THEORYc, RI2_THEORYc = theory_rates(abs(Aee), abs(Aee), abs(Aie), abs(Aie), abs(Aie_NL), abs(Aie_NL), abs(Aei), abs(Aei), abs(Aii), abs(Aii), s_strength, 0.)

Input_E1, Input_E2, Input_I1, Input_I2 = estimate_I(SEE, SEI, SIE, SIEL, SII, s_strength, 100) ###needs to be adjusted for when we add feedforward input to the inhibitory neurons

sEm1 = zeros(50)
sIm1 = zeros(50)
sEm2 = zeros(50)
sIm2 = zeros(50)

for i = 1:50
    es1_ = SEE[i,:][:] .+ SEI[i,:][:] .+ s_strength
    is1_ = SIE[i,:][:] .+ SIEL[i,:][:] .+ SII[i,:][:]
    es2_ = SEE[i+50,:][:] .+ SEI[i+50,:][:] .+ s_strength
    is2_ = SIE[i+50,:][:] .+ SIEL[i+50,:][:] .+ SII[i+50,:][:]
    sEm1[i] = mean(es1_)
    sIm1[i] = mean(is1_)
    sEm2[i] = mean(es2_)
    sIm2[i] = mean(is2_)
end

seL1 = zeros(50)
seL2 = zeros(50)

for i = 1:50
 seL1[i] = mean(SIEL[i,:][:])
 seL2[i] = mean(SIEL[i+50,:][:])
end

push!(SL1, mean(seL1))
push!(SL2, mean(seL2))

push!(E1mn, mean(sEm1))
push!(I1mn, mean(sIm1))
push!(E2mn, mean(sEm2))
push!(I2mn, mean(sIm2))

gti_e1 = estimated_gtiX(Input_E1)
gti_e2 = estimated_gtiX(Input_E2)
gti_i1 = estimated_gtiX(Input_I1)
gti_i2 = estimated_gtiX(Input_I2)

b1 = get_b_g(abs(WEE_1), abs(WII_1), gti_e1, gti_i1)
b2 = get_b_g(abs(WEE_2), abs(WII_2), gti_e2, gti_i2)

c1 = get_c_g(abs(WEE_1), abs(WEI_1), abs(WIE_1), abs(WII_1), gti_e1, gti_i1)
c2 = get_c_g(abs(WEE_2), abs(WEI_2), abs(WIE_2), abs(WII_2), gti_e2, gti_i2)

d1 = get_d_g(abs(WEI_1), abs(WIEL_1), gti_e1, gti_i1)
d2 = get_d_g(abs(WEI_2), abs(WIEL_2), gti_e2, gti_i2)

push!(ERT1, RE1_THEORY)
push!(ERT2, RE2_THEORY)
push!(IRT1, RI1_THEORY)
push!(IRT2, RI2_THEORY)

push!(ERS1, MER1)
push!(ERS2, MER2)
push!(IRS1, MIR1)
push!(IRS2, MIR2)

push!(b1L, b1)
push!(b2L, b2)
push!(c1L, c1)
push!(c2L, c2)
push!(d1L, d1)
push!(d2L, d2)

rex1, rix1 = theory_rates_2x2(abs(WEE_1), abs(WIE_1), abs(WEI_1), abs(WII_1), FE, FI)
rex2, rix2 = theory_rates_2x2(abs(WEE_1), abs(WIE_1), abs(WEI_1), abs(WII_1), FE, FI)
push!(REX1, rex1)
push!(REX2, rex2)
push!(RIX1, rix1)
push!(RIX2, rix2)

push!(wta, wta_ness)

# if wta_ness > .8
#     rex, rix = theory_rates_2x2(abs(WEE_1), abs(WIE_1), abs(WEI_1), abs(WII_1), FE, FI)
#     # push!(REX, rex)
#     # push!(RIX, rix)
# end


# println("\n\n\n $(MER1), $(MER2), $(MIR1), $(MIR2), $(RE1_THEORY), $(RE2_THEORY), $(RI1_THEORY), $(RI2_THEORY), $(EIG1), $(EIG2)\n\n\n")
# push!(ECV1, mean(CV_ETOP))
# push!(ECV2, mean(CV_EBOT))
# push!(ICV1, mean(CV_ITOP))
# push!(ICV2, mean(CV_IBOT))
#
# end

# Aeis = [20,30,40,50,60]
# RET .*= 1000.
# RIT .*= 1000.
# RES .*= 1000.
# RIS .*= 1000.

# plot(RET, RES, ".", ms = 20., label = "E Cells")
# plot(RIT, RIS, ".", ms = 20., label = "I Cells")
# xticks(fontsize = 16)
# yticks(fontsize = 16)
# xlabel("Theory", fontsize = 16)
# ylabel("Simulation", fontsize = 16)
# title("Scanning Aie Local", fontsize = 16)
# plot([0.1, 10.1], [0.1, 10.1], color = "r", label = "x=y")
# legend()

# I_ = zeros(N);
# for i = 1:N
#     I_[i] = mean(Input[i,:][:])
# end

# push!(wta, wta_ness)
# push!(ER1, mean(E_R_bot))
# push!(ER2, mean(E_R_top))
# push!(IR1, mean(I_R_bot))
# push!(IR2, mean(I_R_top))
# push!(ECV1, mean(CV_EBOT))
# push!(ECV2, mean(CV_ETOP))
# push!(ICV1, mean(CV_IBOT))
# push!(ICV2, mean(CV_ITOP))
# push!(means, mean(I_[1:Ne2][:]))


# function sim_2_theory(fe, fi, Aee, Aei, Aie, AieNL, Aii, re, ri)
#     WEE1 = Aee/re

# wta = []
# ER1 = []
# ER2 = []
# IR1 = []
# IR2 = []
# ECV1 = []
# ECV2 = []
# ICV1 = []
# ICV2 = []

#
# Aie_NLS = [0, 30, 60, 90]



# sEm = zeros(100)
# sIm = zeros(100)
# sEs = zeros(100)
# sIs = zeros(100)
#
# for i = 1:100
#     es_ = SEE[i,:][:] .+ SEI[i,:][:] .+ s_strength
#     is_ = SIE[i,:][:] .+ SIEL[i,:][:] .+ SII[i,:][:]
#     sEm[i] = mean(es_)
#     sEs[i] = std(es_)
#     sIm[i] = mean(is_)
#     sIs[i] = std(is_)
# end



# gti_e1 = .1
# gti_i1 = .1
#
#
# WEEz = [abs(i) for i in WEE]
# WEIz = [abs(i) for i in WEI]
# WIEz = [abs(i) for i in WIE]
# WIELz = [abs(i) for i in WIEL]
# WIIz = [abs(i) for i in WII]
#
# slopes = linspace(.1, 1., 10)
# for i in slopes
# gti_e1 = i
# gti_i1 = i
# b1 = get_b_g(WEEz, WIIz, gti_e1, gti_i1)
# c1 = get_c_g(WEEz, WEIz, WIEz, WIIz, gti_e1, gti_i1)
# d1 = get_d_g(WEIz, WIELz, gti_e1, gti_i1)
# b1 =convert(Array{Float64}, b1)
# c1 =convert(Array{Float64}, c1)
# d1 =convert(Array{Float64}, d1)
# r = b1 .+ sqrt(complex(c1 .+ d1))
# plot(WIEL, r, ".", label = "g~'=$(i)")
# end
# legend()
# xlabel("Theoretical WIE_LONG")
# ylabel("++ Eigenvalue")
# title("Slopes consistent with WTA transition")
# axvline(WIEL[6], linestyle = "dashed", alpha = .2)
#
#
# L = zeros(50)
# L2 = zeros(50)
# for i = 1:50
#     L[i] = mean(SIEL[i,:])
#     L2[i] = mean(SIEL[i+50,:])
# end

# FEz = [i-1. for i in Ss]
# FIz = [i-1. for i in Ss]
#
#
# plot(FIz, ERS1, ".", label = "pool 1")
# plot(FIz, ERS2, ".", label = "pool 2")
# xlabel("FI")
# ylabel("rate (Hz)")
# title("Finding g~': Inhibitory Neurons")
# legend()





# plot(Aie_NLS, ERS1, "b", ms = 15., label = "Sim, Pool E1")
# plot(Aie_NLS, ERS2, "r", ms = 15., label = "Sim, Pool E2")
# plot(Aie_NLS, REX1, "b|", ms = 15., label = "2x2 theory, Pool E1")
# plot(Aie_NLS, REX2, "r|", ms = 15., label = "2x2 theory, Pool E2")
# plot(Aie_NLS, ERT1, "b_", ms = 15., label = "4x4 theory, Pool E1")
# plot(Aie_NLS, ERT2, "r_", ms = 15., label = "4x4 theory, Pool E2")
#
# LR1 = WIEL .* ERS1
# LR2 = WIEL .* ERS2
#
# EX1 = []
# IX1 = []
# EX2 = []
# IX2 = []

# for i = 1:length(LR1)
#  ex1, ix1 = theory_rates_2x2(abs(WEE[1]), abs(WIE[1]), abs(WEI[1]), abs(WII[1]), FE, FI + LR1[i])
#  ex2, ix2 = theory_rates_2x2(abs(WEE[1]), abs(WIE[1]), abs(WEI[1]), abs(WII[1]), FE, FI + LR2[i])
#  push!(EX1, ex1)
#  push!(IX1, ix1)
#  push!(EX2, ex2)
#  push!(IX2, ix2)
# end

# LE1, LI1 = theory_rates_2x2(abs(WEE[1]), abs(WIE[1]), abs(WEI[1]), abs(WII[1]), FE, FI + LR2)
# LE2, LI2 = theory_rates_2x2(abs(WEE[1]), abs(WIE[1]), abs(WEI[1]), abs(WII[1]), FE, FI + LR1)
#
# plot(Aie_NLS, REX1, "b", alpha = 0.2, ms = 15., label = "Sim, Pool E1, no altered input")
# plot(Aie_NLS, LE1, "b.", ms = 15., label = "2x2 theory, Pool E1, adjusted FI")
# plot(Aie_NLS, LE2, "r.", ms = 15., label = "2x2 theory, Pool E2, adjusted FI")
# plot(Aie_NLS, ERS1, "b+", ms = 15., label = "Sim, Pool E1")
# plot(Aie_NLS, ERS2, "r+", ms = 15., label = "Sim, Pool E2")
# plot(Aie_NLS, ERT1, "b_", ms = 15., label = "4x4 theory, Pool E1")
# plot(Aie_NLS, ERT2, "r_", ms = 15., label = "4x4 theory, Pool E2")
# legend()
# xlabel("Aie_LONG")
# ylabel("Rate (Hz)")
# title("All predictions vs simulation")
# rex1, rix1 = theory_rates_2x2(abs(WEE_1), abs(WIE_1), abs(WEI_1), abs(WII_1), FE, FI)
# rex2, rix2 = theory_rates_2x2(abs(WEE_2), abs(WIE_2), abs(WEI_2), abs(WII_2), FE, FI)

# plot(E1mn, ERS1, ".", ms = 20., label = "Excitatory Pool 1")
# plot(E2mn, ERS2, ".", ms = 20., label = "Excitatory Pool 2")
# axvline(E1mn[8], linestyle = "dashed")
# axvline(E2mn[8], linestyle = "dashed")
# legend()
# title("Mean Input vs. Firing Rate: Excitatory Cells")
# ylabel("Mean Firing Rate (Hz)")
# subplot(212, sharex = sp)
# plot(I1mn, IRS1, ".", ms = 20., label = "Inhibitory Pool 1")
# plot(I2mn, IRS2, ".", ms = 20., label = "Inhibitory Pool 2")
# legend()
# ylabel("Mean Firing Rate (Hz)")
# xlabel("Mean Input (mV/ms)")
# title("Mean Input vs. Firing Rate: Inhibitory Cells")
#
# big_shot = [max(ERS1[i], ERS2[i]) for i = 1:length(ERS1)]
# lil_shot = [min(ERS1[i], ERS2[i]) for i = 1:length(ERS1)]
#
# big_one = zeros(length(ERS1))
# lil_one = zeros(length(ERS1))
#
# for i = 1:length(ERS1)
#  if ERS1[i] > ERS2[i]
#   big_one[i] = 1
#   lil_one[i] = 2
#  else
#   big_one[i] = 2
#   lil_one[i] = 1
#  end
# end
#
# bolth = zeros(2, length(ERS1))
# for i= 1:length(ERS1)
#  bolth[1, i] = E1mn[i]
#  bolth[2,i] = E2mn[i]
# end
#
# big_input = zeros(length(ERS1))
# lil_input = zeros(length(ERS1))
# for i = 1:length(ERS1)
#  big_input[i] = bolth[big_one[i], i]
#  lil_input[i] = bolth[lil_one[i], i]
# end
#
# plot(big_input, big_shot, ".", ms = 20., label = "more active pool")
# plot(lil_input, lil_shot, ".", ms = 20., label = "less active pool")
# axvline(big_input[8], linestyle = "dashed")
# axvline(lil_input[8], linestyle = "dashed")
#
# xlabel("Mean Input")
# legend()
# ylabel("Firing Rate")
# title("Mean Input vs. Firing Rate: Excitatory Cells")

# plot(Aie_NLS, ERS1, ".", ms = 20., label = "Excitatory Pool 1")
# plot(Aie_NLS, ERS2, ".", ms = 20., label = "Excitatory Pool 2")
#
# plot(Aie_NLS, IRS1, ".", ms = 20., label = "Inhibitory Pool 1")
# plot(Aie_NLS, IRS2, ".", ms = 20., label = "Inhibitory Pool 2")
#
# plot(Aie_NLS, I1mn, ".", ms = 20., label = "Inhibitory Pool 1")
# plot(Aie_NLS, I2mn, ".", ms = 20., label = "Inhibitory Pool 2")


# LR1 = WIEL1 .* ERS1
# LR2 = WIEL2 .* ERS2

# plot(SL1, LR1)
# plot(SL2, LR2)



# LE1, LI1 = theory_rates_2x2(abs(WEE1[1]), abs(WIE1[1]), abs(WEI1[1]), abs(WII1[1]), FE, FI .+ LR2)
# LE2, LI2 = theory_rates_2x2(abs(WEE2[1]), abs(WIE2[1]), abs(WEI2[1]), abs(WII2[1]), FE, FI .+ LR1)
#
# LE1, LI1 = theory_rates_2x2(abs(WEE1[1]), abs(WIE1[1]), abs(WEI1[1]), abs(WII1[1]), FE, FI .+ SL1)
# LE2, LI2 = theory_rates_2x2(abs(WEE2[1]), abs(WIE2[1]), abs(WEI2[1]), abs(WII2[1]), FE, FI .+ SL2)
# #
# plot(Aie_NLS, REX1, "b", alpha = 0.2, ms = 15., label = "2x2 theory, Pool E1, no altered input")
# plot(Aie_NLS, LE1, "b+", ms = 15., label = "2x2 theory, Pool E1, adjusted FI")
# plot(Aie_NLS, LE2, "r+", ms = 15., label = "2x2 theory, Pool E2, adjusted FI")
# plot(Aie_NLS, ERS1, "b.", ms = 15., label = "Sim, Pool E1")
# plot(Aie_NLS, ERS2, "r.", ms = 15., label = "Sim, Pool E2")
# plot(Aie_NLS, ERT1, "bx", ms = 15., label = "4x4 theory, Pool E1")
# plot(Aie_NLS, ERT2, "rx", ms = 15., label = "4x4 theory, Pool E2")
# legend(loc=2)
# xlabel("Aie_LONG")
# ylabel("Rate (Hz)")
# title("All predictions vs simulation")
# rex1, rix1 = theory_rates_2x2(abs(WEE_1), abs(WIE_1), abs(WEI_1), abs(WII_1), FE, FI)
# rex2, rix2 = theory_rates_2x2(abs(WEE_2), abs(WIE_2), abs(WEI_2), abs(WII_2), FE, FI)

# plot(Aie_NLS, IRS1, "b+", ms = 15., label = "Sim, Pool I1")
# plot(Aie_NLS, IRS2, "r+", ms = 15., label = "Sim, Pool I2")
# plot(Aie_NLS, LI1, "b.", ms = 15., label = "2x2 theory, Pool I1, adjusted FI")
# plot(Aie_NLS, LI2, "r.", ms = 15., label = "2x2 theory, Pool I2, adjusted FI")
# legend()
# xlabel("Aie_LONG")
# ylabel("Mean Firing Rate")
# title("2x2 Theory vs. 4x4 Sim: Inhibitory Cells")
# #
#
#
#
# plot(Aie_NLS, IRS1, "r+", ms = 15., label = "Sim, Pool I1")
# # plot(Aie_NLS, IRS2, "m+", ms = 15., label = "Sim, Pool I2")
# plot(Aie_NLS, ERS1, "b+", ms = 15., label = "Sim, Pool E1")
# # plot(Aie_NLS, ERS2, "c+", ms = 15., label = "Sim, Pool E2")
#
# plot(Aie_NLS, LI1, "r.", ms = 15., label = "2x2 theory, Pool I1")
# # plot(Aie_NLS, LI2, "m.", ms = 15., label = "2x2 theory, Pool I2")
# plot(Aie_NLS, LE1, "b.", ms = 15., label = "2x2 theory, Pool E1")
# # plot(Aie_NLS, LE2, "c.", ms = 15., label = "2x2 theory, Pool E2")
# legend()


# big_inh = zeros(length(ERS1))
# big_exc = zeros(length(ERS1))
# lil_inh = zeros(length(ERS1))
# lil_exc = zeros(length(ERS1))
# for i = 1:length(ERS1)
#  if ERS1[i] > ERS2[i]
#   big_inh[i] = IRS1[i]
#   big_exc[i] = ERS1[i]
#   lil_inh[i] = IRS2[i]
#   lil_exc[i] = ERS2[i]
#  else
#   big_inh[i] = IRS2[i]
#   big_exc[i] = ERS2[i]
#   lil_inh[i] = IRS1[i]
#   lil_exc[i] = ERS1[i]
#  end
# end
#
# plot(Aie_NLS, big_exc, ".", ms = 15., label = "E winner")
# plot(Aie_NLS, lil_exc, ".", ms = 15., label = "E loser")
# plot(Aie_NLS, big_inh, ".", ms = 15., label = "I winner")
# plot(Aie_NLS, lil_inh, ".", ms = 15., label = "I loser")
# legend()
#
#
# SL_loser = [max(SL1[i], SL2[i]) for i = 1:length(ERS1)] #if you got more long range input, it's because you lose
# SL_winner = [min(SL1[i], SL2[i]) for i = 1:length(ERS1)]
#
# LEL, LIL = theory_rates_2x2(abs(WEE1[1]), abs(WIE1[1]), abs(WEI1[1]), abs(WII1[1]), FE, FI .+ SL_loser)
# LEW, LIW = theory_rates_2x2(abs(WEE2[1]), abs(WIE2[1]), abs(WEI2[1]), abs(WII2[1]), FE, FI .+ SL_winner)
#
# plot(Aie_NLS, LEW, ".", ms = 15., label = "2x2 E Winner")
# plot(Aie_NLS, LIW, ".", ms = 15., label = "2x2 I Winner")
# plot(Aie_NLS, big_exc, ".", ms = 15., label = "Sim E Winner")
# plot(Aie_NLS, big_inh, ".", ms = 15., label = "Sim I Winner")
# legend()
#
# #want to isolate inhibitory inputs in dominant circuit
#
# big_I = zeros(length(ERS1))
# lil_I = zeros(length(ERS1))
#
# for i = 1:length(ERS1)
#  if ERS1[i] > ERS2[i]
#   big_I[i] = 1
#   lil_I[i] = 2
#  else
#   big_I[i] = 2
#   lil_I[i] = 1
#  end
# end
#
# bolth = zeros(2, length(ERS1))
# for i= 1:length(ERS1)
#  bolth[1, i] = I1mn[i]
#  bolth[2,i] = I2mn[i]
# end
#
# big_in = zeros(length(ERS1))
# lil_in = zeros(length(ERS1))
# for i = 1:length(ERS1)
#  big_in[i] = bolth[big_I[i], i]
#  lil_in[i] = bolth[lil_I[i], i]
# end
#
# plot(big_in, big_inh, ".", ms = 15., label = "Inhibitory Pool, Dominant Circuit")
# plot(lil_in, lil_inh, ".", ms = 15., label = "Inhibitory Pool, Supressed Circuit")
# axvline(big_in[10], linestyle = "dashed")
# axvline(lil_in[10], linestyle = "dashed")
#
#

### Threshold Adjusted Predictions
s_strength = 3.08

### Get rate, input, for dominant pool

big_inh = zeros(length(ERS1)) #inhibitory rate in dominant CIRCUIT
big_exc = zeros(length(ERS1)) #excitatory rate in dominant CIRCUIT
lil_inh = zeros(length(ERS1)) #same, for suppressed circuits
lil_exc = zeros(length(ERS1))
for i = 1:length(ERS1)
 if ERS1[i] > ERS2[i]
  big_inh[i] = IRS1[i]
  big_exc[i] = ERS1[i]
  lil_inh[i] = IRS2[i]
  lil_exc[i] = ERS2[i]
 else
  big_inh[i] = IRS2[i]
  big_exc[i] = ERS2[i]
  lil_inh[i] = IRS1[i]
  lil_exc[i] = ERS1[i]
 end
end


big_inh_input = zeros(length(ERS1))
big_exc_input = zeros(length(ERS1))
lil_inh_input = zeros(length(ERS1))
lil_exc_input = zeros(length(ERS1))
for i = 1:length(ERS1)
 if ERS1[i] > ERS2[i]
  big_inh_input[i] = I1mn[i]
  big_exc_input[i] = E1mn[i]
  lil_inh_input[i] = I2mn[i]
  lil_exc_input[i] = E2mn[i]
 else
  big_inh_input[i] = I2mn[i]
  big_exc_input[i] = E2mn[i]
  lil_inh_input[i] = I1mn[i]
  lil_exc_input[i] = E1mn[i]
 end
end

big_LR = zeros(length(ERS1))
lil_LR = zeros(length(ERS1))
for i = 1:length(ERS1)
 if ERS1[i] > ERS2[i]
  big_LR[i] = SL1[i]
  lil_LR[i] = SL2[i]
 else
  big_LR[i] = SL2[i]
  lil_LR[i] = SL1[i]
 end
end

# LE1W, LI1W = theory_rates_2x2(abs(WEE1[1]), abs(WIE1[1]), abs(WEI1[1]), abs(WII1[1]), s_strength .- big_exc_input, 0 .- big_inh_input)
# LE2L, LI2L = theory_rates_2x2(abs(WEE2[1]), abs(WIE2[1]), abs(WEI2[1]), abs(WII2[1]), s_strength .- lil_exc_input, 0 .- lil_inh_input)
#
# plot(Aie_NLS, big_exc, "bx", ms = 15, label = "E_winner Sim")
# plot(Aie_NLS, lil_exc, "rx", ms = 15, label = "E_loser Sim")
# plot(Aie_NLS, LE1W, "b.", ms = 15, label = "E_winner 2x2 theory")
# plot(Aie_NLS, LE2L, "r.", ms = 15, label = "E_loser 2x2 theory")
#
# plot(Aie_NLS, big_inh, "b+", ms = 15., label = "I_winner Sim")
# plot(Aie_NLS, lil_inh, "r+", ms = 15., label = "I_loser Sim")
# plot(Aie_NLS, LI1W, "b*", ms = 15, label = "I_winner 2x2 theory")
# plot(Aie_NLS, LI2L, "r*", ms = 15, label = "I_loser 2x2 theory")
# legend()
#
# ##########dominant only, then suppressed only
# plot(Aie_NLS, big_exc, "b+", ms = 15, label = "E_winner Sim")
# plot(Aie_NLS, LE1W, "b.", ms = 15, label = "E_winner 2x2 theory")
# plot(Aie_NLS, big_inh, "r+", ms = 15., label = "I_winner Sim")
# plot(Aie_NLS, LI1W, "r.", ms = 15, label = "I_winner 2x2 theory")
# legend()
#
#
#
# plot(Aie_NLS, lil_exc, "b+", ms = 15, label = "E_loser Sim")
# plot(Aie_NLS, LE2L, "b.", ms = 15, label = "E_loser 2x2 theory")
# plot(Aie_NLS, lil_inh, "r+", ms = 15., label = "I_loser Sim")
# plot(Aie_NLS, LI2L, "r.", ms = 15, label = "I_loser 2x2 theory")
# legend()
#
WEEz = [abs(i) for i in WEE1]
WEIz = [abs(i) for i in WEI1]
WIEz = [abs(i) for i in WIE1]
WIELz = [abs(i) for i in WIEL1]
WIIz = [abs(i) for i in WII1]


REW_4x4, RIW_4x4 = theory_rates_4x4_1C(WEEz, WIEz, WIELz, WEIz, WIIz, s_strength .- big_exc_input, 0 .- big_inh_input)
REL_4x4, RIL_4x4 = theory_rates_4x4_1C(WEEz, WIEz, WIELz, WEIz, WIIz, s_strength .- lil_exc_input, 0 .- lil_inh_input)

# plot(Aie_NLS, big_exc, "b+", ms = 15, label = "E_winner Sim")
# plot(Aie_NLS, REW_4x4, "b.", ms = 15, label = "E_winner 4x4 theory")
# plot(Aie_NLS, big_inh, "r+", ms = 15, label = "I_winner Sim")
# plot(Aie_NLS, RIW_4x4, "r.", ms = 15, label = "I_winner 4x4 theory")
# axvline(Aie_NLS[10], linestyle="dashed")
# legend()
# xlabel("Aie_LONG")
# ylabel("Rate (Hz)")
# title("Threshold Adjusted 4x4 Theory: More Active Pool")
#
#
# plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
# plot(Aie_NLS, REW_4x4, "bx", ms = 15, label = "E_winner 4x4 theory")
# plot(Aie_NLS, big_inh, "r.", ms = 15, label = "I_winner Sim")
# plot(Aie_NLS, RIW_4x4, "rx", ms = 15, label = "I_winner 4x4 theory")
# plot(Aie_NLS, LE1W, "b+", ms = 15, label = "E_winner 2x2 theory")
# plot(Aie_NLS, LI1W, "r+", ms = 15, label = "I_winner 2x2 theory")
# axvline(Aie_NLS[28], linestyle="dashed")
# legend()
# xlabel("Aie_LONG")
# ylabel("Rate (Hz)")
# title("Threshold Adjusted 4x4 Theory: More Active Pool")
#
# TE = [-.227 for i = 1:length(IRS1)]
# TI = [-.1 for i = 1:length(IRS1)]
#
# LE1W, LI1W = theory_rates_2x2(abs(WEE1[1]), abs(WIE1[1]), abs(WEI1[1]), abs(WII1[1]), s_strength .- TE, 0 .- TI)
# LE2L, LI2L = theory_rates_2x2(abs(WEE2[1]), abs(WIE2[1]), abs(WEI2[1]), abs(WII2[1]), s_strength .- lil_exc_input, 0 .- lil_inh_input)
# LE2L, LI2L = theory_rates_2x2(abs(WEE2[1]), abs(WIE2[1]), abs(WEI2[1]), abs(WII2[1]), s_strength .- TE, 0 .- TI)
#
# REW_4x4, RIW_4x4 = theory_rates_4x4_1C(WEEz, WIEz, WIELz, WEIz, WIIz, s_strength .- TE, 0 .- TI)
# REL_4x4, RIL_4x4 = theory_rates_4x4_1C(WEEz, WIEz, WIELz, WEIz, WIIz, s_strength .- lil_exc_input, 0 .- lil_inh_input)
# REL_4x4, RIL_4x4 = theory_rates_4x4_1C(WEEz, WIEz, WIELz, WEIz, WIIz, s_strength .- TE, 0 .- TI)

###Test Suppressed Circuit Rates against 2x2 theory with and without adjusting fi with long range input

# plot(Aie_NLS, lil_exc, "b.", ms = 15, label = "E_loser Sim")
# plot(Aie_NLS, REW_4x4, "bx", ms = 15, label = "E_loser 4x4 theory")
# plot(Aie_NLS, lil_inh, "r.", ms = 15, label = "I_loser Sim")
# plot(Aie_NLS, RIW_4x4, "rx", ms = 15, label = "I_loser 4x4 theory")
# plot(Aie_NLS, LE1W, "b+", ms = 15, label = "E_loser 2x2 theory")
# plot(Aie_NLS, LI1W, "r+", ms = 15, label = "I_loser 2x2 theory")
# axvline(Aie_NLS[28], linestyle="dashed")
# legend()
# xlabel("Aie_LONG")
# ylabel("Rate (Hz)")
# title("Threshold Adjusted Theory: Less Active Pool")
#
#
# #modify fi with long range input
# LE1Wi, LI1Wi = theory_rates_2x2(abs(WEE1[1]), abs(WIE1[1]), abs(WEI1[1]), abs(WII1[1]), s_strength .- TE, 0 .- TI .+ big_LR)
# LE2Li, LI2Li = theory_rates_2x2(abs(WEE2[1]), abs(WIE2[1]), abs(WEI2[1]), abs(WII2[1]), s_strength .- TE, 0 .- TI .+ lil_LR)
#
# plot(Aie_NLS, lil_exc, "b.", ms = 15, label = "E_loser Sim")
# plot(Aie_NLS, REW_4x4, "bx", ms = 15, label = "E_loser 4x4 theory")
# plot(Aie_NLS, lil_inh, "r.", ms = 15, label = "I_loser Sim")
# plot(Aie_NLS, RIW_4x4, "rx", ms = 15, label = "I_loser 4x4 theory")
# plot(Aie_NLS, LE2Li, "b+", ms = 15, label = "E_loser 2x2 theory")
# plot(Aie_NLS, LI2Li, "r+", ms = 15, label = "I_loser 2x2 theory")
# axvline(Aie_NLS[28], linestyle="dashed")
# legend()
# xlabel("Aie_LONG")
# ylabel("Rate (Hz)")
# title("Input Adjusted Theory: Less Active Pool")
# #
# plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
# plot(Aie_NLS, REW_4x4, "bx", ms = 15, label = "E_winner 4x4 theory")
# plot(Aie_NLS, big_inh, "r.", ms = 15, label = "I_winner Sim")
# plot(Aie_NLS, RIW_4x4, "rx", ms = 15, label = "I_winner 4x4 theory")
# plot(Aie_NLS, LE1Wi, "b+", ms = 15, label = "E_winner 2x2 theory")
# plot(Aie_NLS, LI1Wi, "r+", ms = 15, label = "I_winner 2x2 theory")
# axvline(Aie_NLS[28], linestyle="dashed")
# legend()
# xlabel("Aie_LONG")
# ylabel("Rate (Hz)")
# title("Threshold Adjusted 4x4 Theory: More Active Pool")



#check if we can predict the instability

# b1 = get_b_g(WEEz, WIIz, s_strength .- TE, -TI)
# c1 = get_c_g(WEEz, WEIz, WIEz, WIIz, s_strength .- TE, -TI)
# d1 = get_d_g(WEIz, WIELz, s_strength .- TE, -TI)
# b1 =convert(Array{Float64}, b1)
# c1 =convert(Array{Float64}, c1)
# d1 =convert(Array{Float64}, d1)
# rpp = b1 .+ sqrt(complex(c1 .+ d1))
# rpm = b1 .+ sqrt(complex(c1 .- d1))
# rmp = b1 .- sqrt(complex(c1 .+ d1))
# rmm = b1 .- sqrt(complex(c1 .- d1))
#
# plot(WIELz, rpp, ".", ms = 15.)
# # plot(WIELz, rpm, ".", ms = 15.)
# # plot(WIELz, rmp, ".", ms = 15.)
# # plot(WIELz, rmm, ".", ms = 15.)
#
# plot(WIELz, rpp, ".", ms = 15.)
# legend()
# xlabel("Theoretical WIE_LONG")
# ylabel("++ Eigenvalue")
# title("Predicting Instability With Adjusted Thresholds")
# axvline(WIELz[28], linestyle = "dashed")
# axhline(0,linestyle = "dashed")



###4x4 with separated feedforward input
# TE1 = [mean(E1mn) for i = 1:length(IRS1)]
# TI1 = [mean(I1mn) for i = 1:length(IRS1)]
# TE2 = [mean(E2mn) for i = 1:length(IRS1)]
# TI2 = [mean(I2mn) for i = 1:length(IRS1)]

TE1 = [-.5 for i = 1:length(IRS1)]
TI1 = [-.12 for i = 1:length(IRS1)]
TE2 = [-.5 for i = 1:length(IRS1)]
TI2 = [-.12 for i = 1:length(IRS1)]

RE1_4x4_sf, RE2_4x4_sf, RI1_4x4_sf, RI2_4x4_sf = theory_rates_4x4_sf(WEEz, WIEz, WIELz, WEIz, WIIz, fe .- TE1, 0 .- TI1, fe2 .- TE2, 0 .- TI2)

plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
plot(Aie_NLS, lil_exc, "r.", ms = 15, label = "E_loser Sim")
plot(Aie_NLS, RE2_4x4_sf, "b+", ms = 15, label = "E Theory win")
plot(Aie_NLS, RE1_4x4_sf, "r+", ms = 15, label = "E Theory los")
legend()
xlabel("Aie_LONG")
ylabel("Rate (Hz)")
title("Separate f 4x4 Theory: Excitatory Neurons")

plot(Aie_NLS, big_inh, "b.", ms = 15, label = "I_winner Sim")
plot(Aie_NLS, lil_inh, "r.", ms = 15, label = "I_loser Sim")
plot(Aie_NLS, RI2_4x4_sf, "b+", ms = 15, label = "I Theory win")
plot(Aie_NLS, RI1_4x4_sf, "r+", ms = 15, label = "I Theory los")
legend()
xlabel("Aie_LONG")
ylabel("Rate (Hz)")
title("Separate f 4x4 Theory: Inhibitory Neurons")

#standard 2x2 theory where fe = fe2
LE1W, LI1W = theory_rates_2x2(abs(WEE1[1]), abs(WIE1[1]), abs(WEI1[1]), abs(WII1[1]), fe2 .- TE1, 0 .- TI1)
LE1L, LI1L = theory_rates_2x2(abs(WEE1[1]), abs(WIE1[1]), abs(WEI1[1]), abs(WII1[1]), fe .- TE1, 0 .- TI1)

#2x2 theory with long range input absorbed into fi
LE1Wi, LI1Wi = theory_rates_2x2(abs(WEE1[1]), abs(WIE1[1]), abs(WEI1[1]), abs(WII1[1]), fe2 .- TE1, 0 .- TI1 .+ big_LR)
LE2Li, LI2Li = theory_rates_2x2(abs(WEE2[1]), abs(WIE2[1]), abs(WEI2[1]), abs(WII2[1]), fe2 .- TE1, 0 .- TI1 .+ lil_LR)
#standard 4x4 theory,
REW_4x4, RIW_4x4 = theory_rates_4x4_1C(WEEz, WIEz, WIELz, WEIz, WIIz, fe2.- TE1, 0 .- TI1)

plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
plot(Aie_NLS, REW_4x4, "bx", ms = 15, label = "E_winner 4x4 theory standard")
plot(Aie_NLS, big_inh, "r.", ms = 15, label = "I_winner Sim")
plot(Aie_NLS, RIW_4x4, "rx", ms = 15, label = "I_winner 4x4 theory standard")
plot(Aie_NLS, LE1Wi, "b+", ms = 15, label = "E_winner 2x2 theory input adjusted")
plot(Aie_NLS, LI1Wi, "r+", ms = 15, label = "I_winner 2x2 theory input adjusted")
legend()
xlabel("Aie_LONG")
ylabel("Rate (Hz)")
title("Threshold Adjusted 4x4 Theory: More Active Pool")

plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
plot(Aie_NLS, REW_4x4, "bx", ms = 15, label = "E_winner 4x4 theory standard")
plot(Aie_NLS, big_inh, "r.", ms = 15, label = "I_winner Sim")
plot(Aie_NLS, RIW_4x4, "rx", ms = 15, label = "I_winner 4x4 theory standard")
plot(Aie_NLS, LE1W, "b+", ms = 15, label = "E_winner 2x2 theory")
plot(Aie_NLS, LI1W, "r+", ms = 15, label = "I_winner 2x2 theory")
legend(loc=6)
xlabel("Aie_LONG")
ylabel("Rate (Hz)")
title("Threshold Adjusted 4x4 Theory: More Active Pool")

plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
plot(Aie_NLS, REW_4x4, "bx", ms = 15, label = "E_winner 4x4 theory standard")
plot(Aie_NLS, big_inh, "r.", ms = 15, label = "I_winner Sim")
plot(Aie_NLS, RIW_4x4, "rx", ms = 15, label = "I_winner 4x4 theory standard")
plot(Aie_NLS, RE2_4x4_sf, "b+", ms = 15, label = "E Theory win")
plot(Aie_NLS, RE1_4x4_sf, "r+", ms = 15, label = "E Theory los")
legend()
xlabel("Aie_LONG")
ylabel("Rate (Hz)")
title("Separate f 4x4 Theory: Excitatory Neurons")



plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
plot(Aie_NLS, REW_4x4, "bx", ms = 15, label = "E_winner 4x4 theory standard")
plot(Aie_NLS, RE2_4x4_sf, "b+", ms = 15, label = "E_winner SF")
plot(Aie_NLS, LE1W, "bv", ms = 15, label = "E_winner 2x2 theory")
plot(Aie_NLS, big_inh, "r.", ms = 15, label = "I_winner Sim")
plot(Aie_NLS, RIW_4x4, "rx", ms = 15, label = "I_winner 4x4 theory standard")
plot(Aie_NLS, RI2_4x4_sf, "r+", ms = 15, label = "I_winner SF")
plot(Aie_NLS, LI1W, "rv", ms = 15, label = "I_winner 2x2 theory")
legend()
xlabel("Aie_LONG")
ylabel("Rate (Hz)")
title("Theories: Winning Circuit")


plot(Aie_NLS, lil_exc, "b.", ms = 15, label = "E_loser Sim")
plot(Aie_NLS, RE1_4x4_sf, "bx", ms = 15, label = "E_loser SF")
plot(Aie_NLS, LE1L, "b+", ms = 15, label = "E_loser 2x2 theory")
plot(Aie_NLS, lil_inh, "r.", ms = 15, label = "I_loser Sim")
plot(Aie_NLS, RI1_4x4_sf, "rx", ms = 15, label = "I_loser SF")
plot(Aie_NLS, LI1L, "r+", ms = 15, label = "I_loser 2x2 theory")
legend()
xlabel("Aie_LONG")
ylabel("Rate (Hz)")
title("Theories: Loser Circuit")


plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
plot(Aie_NLS, RE2_4x4_sf, "bv", ms = 15, label = "E_winner SF")
plot(Aie_NLS, LE1W, "b+", ms = 15, label = "E_winner 2x2 theory")
plot(Aie_NLS, big_inh, "r.", ms = 15, label = "I_winner Sim")
plot(Aie_NLS, RI2_4x4_sf, "rv", ms = 15, label = "I_winner SF")
plot(Aie_NLS, LI1W, "r+", ms = 15, label = "I_winner 2x2 theory")
legend()
xlabel("Aie_LONG")
ylabel("Rate (Hz)")
title("Theories: Winning Circuit")


#
plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
plot(Aie_NLS, RE2_4x4_sf, "bv", ms = 15, label = "E_winner SF")
plot(Aie_NLS, LE1W, "b", alpha = .2, label = "E_winner 2x2 theory")
plot(Aie_NLS, REW_4x4, "bx", ms = 15, label = "E_winner 4x4 theory standard")
plot(Aie_NLS, LE1Wi, "b+", ms = 15, label = "E_winner 2x2 theory input adjusted")


# RE1_4x4_sf_dense, RE2_4x4_sf_dense, RI1_4x4_sf_dense, RI2_4x4_sf_dense = theory_rates_4x4_sf(WEEzd, WIEzd, WIELzd, WEIzd, WIIzd, fe .- TE1d, 0 .- TI1d, fe +.01 .- TE2d, 0 .- TI2d)
# #
# TE1d = [-.5 for i = 1:100]
# TE2d = [-.5 for i = 1:100]
# TI1d = [-.12 for i = 1:100]
# TI2d = [-.12 for i = 1:100]
#
# WEEzd = [.275 for i=1:100]
# WIEzd = [1.96 for i=1:100]
# WIELzd = linspace(0, 2.2, 100)
# WEIzd = [2.20 for i=1:100]
# WIIzd = [1.10 for i=1:100]
#
# plot(Aie_NLS, big_exc, "b.", ms = 15, label = "E_winner Sim")
# plot(Aie_NLS, big_inh, "r.", ms = 15, label = "I_winner Sim")
# plot(WIELzd .* 36.36, RE2_4x4_sf_dense, "bx", label = "")


# RE1_4x4_sf, RE2_4x4_sf, RI1_4x4_sf, RI2_4x4_sf = theory_rates_4x4_sf(WEEz, WIEz, WIELz, WEIz, WIIz, fe .- TE1, 0 .- TI1, fe+.2 .- TE2, 0 .- TI2)
#
# plot(WIELzd, RE1_4x4_sf_dense, "b.", ms = 15.)
# plot(WIELzd, RE2_4x4_sf_dense, "r.", ms = 15.)
