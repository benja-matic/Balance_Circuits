include("4x4_separate_synapses.jl")
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

Aee = 50.
Aei = 200.
Aie = 500.
Aie_NL = 600.
Aii = 200.

#Aie = [200, 400, 600, 800, 1000]

s_strength = 3.08
p = .2

N = 4000
N2 = div(N, 2)
NeL = div(4*N2, 5)
NiL = N2-NeL
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
end_trans = 2000
rt = ((ntotal - end_trans)/1000.)*h

W = homogenous_4x4_weights(N, p, Aee, Aei, Aie, Aie_NL, Aii);
CSR = sparse_rep(W, N);

@time te, re, ti, ri, SEE, SEI, SIE, SIEL, SII, LX = euler_lif_CSR_4x4_s(h, runtime, N, W, CSR, s_strength, vth, tau_m, tau_s, tau_a, g_a)

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

NSET = length(find(re .> NeL))
NSEB = length(find(re .<= NeL))
NSIT = length(find(ri .> NiL))
NSIB = length(find(ri .<= NiL))

MER1 = mean(E_R_bot)*h*(1/1000.)
MER2 = mean(E_R_top)*h*(1/1000.)
MIR1 = mean(I_R_bot)*h*(1/1000.)
MIR2 = mean(I_R_top)*h*(1/1000.)

println("##RESULT $(wta_ness), $(bias), $(mean(E_R_bot)), $(mean(E_R_top)), $(mean(I_R_top)), $(mean(I_R_bot)), $(mean(CV_ETOP)), $(mean(CV_EBOT)), $(mean(CV_ITOP)), $(mean(CV_IBOT))")

WEE_1, WEE_2, WIE_1, WIE_2, WIEL_1, WIEL_2, WEI_1, WEI_2, WII_1, WII_2, FE, FI = sim_2_theory(SEE, SEI, SIE, SIEL, SII, s_strength*h, 0, .1, MER1, MER2, MIR1, MIR2, 100)

println("##PARAMETERS $(WEE_1), $(WEE_2), $(WIE_1), $(WIE_2), $(WIEL_1), $(WIEL_2), $(WEI_1), $(WEI_2), $(WII_1), $(WII_2), $(FE), $(FI)")





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
