include("4x4_equal_EI.jl")
include("Analyze.jl")

srand(4321)


# wta = []
# ER1 = []
# ER2 = []
# IR1 = []
# IR2 = []
# ECV1 = []
# ECV2 = []
# ICV1 = []
# ICV2 = []
# means = []

RET = []
RIT = []
RES = []
RIS = []
#
# Aie = [200, 300, 400, 500, 600, 700, 800, 900, 1000]
# for i in Aie
Aee = 25.
Aei = 300.
Aie = 300.
Aie_NL = 0.
Aii = 200.

#Aie = [200, 400, 600, 800, 1000]

s_strength = 3.08
p = .2

N = 16000

N_local = Int64(round(N/2)) #divide the network into two EI circuits
half = round(Int64, N_local/2)

# N2 = div(N, 2)
# NeL = div(4*N2, 5)
# NiL = N2-NeL
# Ne2 = NeL*2
# Ni2 = NiL*2

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

@time te, re, ti, ri, SEE, SEI, SIE, SIEL, SII = euler_lif_CSR_4x4_s(h, runtime, N, W, CSR, s_strength, vth, tau_m, tau_s, tau_a, g_a)

wta_ness, bias = score_analysis(re, N_local)
top_e_neurons, bot_e_neurons = Neurons_tb_ns(re, half, 10, min_e_neurons)
top_i_neurons, bot_i_neurons = Neurons_tb_ns(ri, half, 10, min_i_neurons)

CV_ETOP = CV_ISI_ALLTIME(top_e_neurons, te, re)
CV_EBOT = CV_ISI_ALLTIME(bot_e_neurons, te, re)
CV_ITOP = CV_ISI_ALLTIME(top_i_neurons, ti, ri)
CV_IBOT = CV_ISI_ALLTIME(bot_i_neurons, ti, ri)

E_R_top = [length(find(re .== i))/rt for i=half+1:N_local]
E_R_bot = [length(find(re .== i))/rt for i=1:half]
I_R_top = [length(find(ri .== i))/rt for i=half+1:N_local]
I_R_bot = [length(find(ri .== i))/rt for i=1:half]

# NSET = length(find(re .> NeL))
# NSEB = length(find(re .<= NeL))
# NSIT = length(find(ri .> NiL))
# NSIB = length(find(ri .<= NiL))

MER1 = mean(E_R_bot)*(1/1000.)
MER2 = mean(E_R_top)*(1/1000.)
MIR1 = mean(I_R_bot)*(1/1000.)
MIR2 = mean(I_R_top)*(1/1000.)

println("##RESULT $(wta_ness), $(bias), $(mean(E_R_bot)), $(mean(E_R_top)), $(mean(I_R_top)), $(mean(I_R_bot)), $(mean(CV_ETOP)), $(mean(CV_EBOT)), $(mean(CV_ITOP)), $(mean(CV_IBOT))")

WEE_1, WEE_2, WIE_1, WIE_2, WIEL_1, WIEL_2, WEI_1, WEI_2, WII_1, WII_2, FE, FI = sim_2_theory(SEE, SEI, SIE, SIEL, SII, s_strength, 0, 1., MER1, MER2, MIR1, MIR2, 100)

println("##PARAMETERS $(WEE_1), $(WEE_2), $(WIE_1), $(WIE_2), $(WIEL_1), $(WIEL_2), $(WEI_1), $(WEI_2), $(WII_1), $(WII_2), $(FE), $(FI)")

RE_THEORY, RI_THEORY = theory_rates(abs(WEE_1), abs(WEE_2), abs(WIE_1), abs(WIE_2), abs(WIEL_1), abs(WIEL_2), abs(WEI_1), abs(WEI_2), abs(WII_1), abs(WII_2), FE, FI)


push!(RET, RE_THEORY)
push!(RIT, RI_THEORY)
push!(RES, MER1)
push!(RIS, MIR1)
# end

# Aie = [200, 300, 400, 500, 600, 700, 800, 900, 1000]
#
# RET .*= 1000.
# RIT .*= 1000.
# RES .*= 1000.
# RIS .*= 1000.


# figure(1)
# plot(RET, RES, ".", ms = 20., label = "E Cells")
# plot(RIT, RIS, ".", ms = 20., label = "I Cells")
# xticks(fontsize = 16)
# yticks(fontsize = 16)
# xlabel("Theory", fontsize = 16)
# ylabel("Simulation", fontsize = 16)
# title("Scanning Aie Local", fontsize = 16)
# plot([0.1, 10.1], [0.1, 10.1], color = "r", label = "x=y")
# legend()

# subplot(211)
# plot(Aie, RET, ".", label = "Theory E cells")
# plot(Aie, RES, ".", label = "Sim E cells")
# legend()
# # xlabel("Aie")
# ylabel("Firing Rate")
# title("Theory vs. Sim: E Cells")
# subplot(212)
# plot(Aie, RIT, ".", label = "Theory I cells")
# plot(Aie, RIS, ".", label = "Sim I cells")
# legend()
# xlabel("Aie")
# ylabel("Firing Rate")
# title("Theory vs. Sim: I Cells")
# tight_layout()
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
# Aie = [200, 300, 400, 500]
# subplot(211)
# plot(Aie, ER1, ".", label = "Excitatory Pool 1")
# plot(Aie, ER2, ".", label = "Excitatory Pool 2")
# plot(Aie, IR1, ".", label = "Inhibitory Pool 1")
# plot(Aie, IR2, ".", label = "Inhibitory Pool 2")
# legend()
# ylabel("Mean Firing Rate")
# subplot(212)
# plot(Aie, ECV1, ".", label = "Excitatory Pool 1")
# plot(Aie, ECV2, ".", label = "Excitatory Pool 2")
# plot(Aie, ICV1, ".", label = "Inhibitory Pool 1")
# plot(Aie, ICV2, ".", label = "Inhibitory Pool 2")
# legend()
# xlabel("Aie NON-LOCAL")
# ylabel("CV_ISI")
