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
WEE = []
WEI = []
WIE = []
WIEL = []
WII = []
REX1 = []
RIX1 = []
REX2 = []
RIX2 = []

E1std = []
I1std = []
E1mn = []
I1mn = []

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
Aie_NLS = [0, 20, 40, 60, 65, 70, 75, 80]
for i in Aie_NLS
Aee = 5.
Aei = 40.
Aie = 70.
Aie_NL = i
Aii = 20.

#Aie = [200, 400, 600, 800, 1000]

s_strength = 3.08
fe = s_strength
fi = 0.
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
runtime = 5000 #ms
h = .1
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
netd_binsize = 50/h
end_trans = 0.
rt = ((ntotal - end_trans)/1000.)*h

W = homogenous_4x4_weights(N, IFRAC, k, Aee, Aei, Aie, Aie_NL, Aii);
CSR = sparse_rep(W, N);

@time te, re, ti, ri, SEE, SEI, SIE, SIEL, SII = euler_lif_CSR_4x4_s(h, runtime, N, IFRAC, W, CSR, fe, fi, vth, tau_m, tau_s, tau_a, g_a)

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

push!(WEE, WEE_1)
push!(WEI, WEI_1)
push!(WIE, WIE_1)
push!(WIEL, WIEL_1)
push!(WII, WII_1)

println("##PARAMETERS $(WEE_1), $(WEE_2), $(WIE_1), $(WIE_2), $(WIEL_1), $(WIEL_2), $(WEI_1), $(WEI_2), $(WII_1), $(WII_2), $(FE), $(FI)")

RE1_THEORY, RI1_THEORY, RE2_THEORY, RI2_THEORY = theory_rates(abs(WEE_1), abs(WEE_2), abs(WIE_1), abs(WIE_2), abs(WIEL_1), abs(WIEL_2), abs(WEI_1), abs(WEI_2), abs(WII_1), abs(WII_2), FE, FI)

# RE1_THEORYc, RI1_THEORYc, RE2_THEORYc, RI2_THEORYc = theory_rates(abs(Aee), abs(Aee), abs(Aie), abs(Aie), abs(Aie_NL), abs(Aie_NL), abs(Aei), abs(Aei), abs(Aii), abs(Aii), s_strength, 0.)

Input_E1, Input_E2, Input_I1, Input_I2 = estimate_I(SEE, SEI, SIE, SIEL, SII, s_strength, 100) ###needs to be adjusted for when we add feedforward input to the inhibitory neurons

sEm = zeros(100)
sIm = zeros(100)
sEs = zeros(100)
sIs = zeros(100)

for i = 1:100
    es_ = SEE[i,:][:] .+ SEI[i,:][:] .+ s_strength
    is_ = SIE[i,:][:] .+ SIEL[i,:][:] .+ SII[i,:][:]
    sEm[i] = mean(es_)
    sEs[i] = std(es_)
    sIm[i] = mean(is_)
    sIs[i] = std(is_)
end

push!(E1std, mean(sEs))
push!(E1mn, mean(sEm))
push!(I1std, mean(sIs))
push!(I1mn, mean(sIm))

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
end

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
FIz = [i-1. for i in Ss]


plot(FIz, ERS1, ".", label = "pool 1")
plot(FIz, ERS2, ".", label = "pool 2")
xlabel("FI")
ylabel("rate (Hz)")
title("Finding g~': Inhibitory Neurons")
legend()





plot(Aie_NLS, ERS1, "b", ms = 15., label = "Sim, Pool E1")
plot(Aie_NLS, ERS2, "r", ms = 15., label = "Sim, Pool E2")
plot(Aie_NLS, REX1, "b|", ms = 15., label = "2x2 theory, Pool E1")
plot(Aie_NLS, REX2, "r|", ms = 15., label = "2x2 theory, Pool E2")
plot(Aie_NLS, ERT1, "b_", ms = 15., label = "4x4 theory, Pool E1")
plot(Aie_NLS, ERT2, "r_", ms = 15., label = "4x4 theory, Pool E2")

LR1 = WIEL .* ERS1
LR2 = WIEL .* ERS2

EX1 = []
IX1 = []
EX2 = []
IX2 = []

# for i = 1:length(LR1)
#  ex1, ix1 = theory_rates_2x2(abs(WEE[1]), abs(WIE[1]), abs(WEI[1]), abs(WII[1]), FE, FI + LR1[i])
#  ex2, ix2 = theory_rates_2x2(abs(WEE[1]), abs(WIE[1]), abs(WEI[1]), abs(WII[1]), FE, FI + LR2[i])
#  push!(EX1, ex1)
#  push!(IX1, ix1)
#  push!(EX2, ex2)
#  push!(IX2, ix2)
# end

LE1, LI1 = theory_rates_2x2(abs(WEE[1]), abs(WIE[1]), abs(WEI[1]), abs(WII[1]), FE, FI + LR2)
LE2, LI2 = theory_rates_2x2(abs(WEE[1]), abs(WIE[1]), abs(WEI[1]), abs(WII[1]), FE, FI + LR1)

plot(Aie_NLS, REX1, "b", alpha = 0.2, ms = 15., label = "Sim, Pool E1, no altered input")
plot(Aie_NLS, LE1, "b.", ms = 15., label = "2x2 theory, Pool E1, adjusted FI")
plot(Aie_NLS, LE2, "r.", ms = 15., label = "2x2 theory, Pool E2, adjusted FI")
plot(Aie_NLS, ERS1, "b+", ms = 15., label = "Sim, Pool E1")
plot(Aie_NLS, ERS2, "r+", ms = 15., label = "Sim, Pool E2")
plot(Aie_NLS, ERT1, "b_", ms = 15., label = "4x4 theory, Pool E1")
plot(Aie_NLS, ERT2, "r_", ms = 15., label = "4x4 theory, Pool E2")
legend()
xlabel("Aie_LONG")
ylabel("Rate (Hz)")
title("2x2 predictions with fi = fi + long range input")
# rex1, rix1 = theory_rates_2x2(abs(WEE_1), abs(WIE_1), abs(WEI_1), abs(WII_1), FE, FI)
# rex2, rix2 = theory_rates_2x2(abs(WEE_2), abs(WIE_2), abs(WEI_2), abs(WII_2), FE, FI)
