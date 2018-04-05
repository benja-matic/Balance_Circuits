include("2x2.jl")
include("Analyze.jl")

srand(4321)

RET = []
RIT = []
RES = []
RIS = []
ECV = []
ICV = []

# Aie = [700, 800, 900, 1000]
# for i in Aie
Aee = 100.
Aei = 300.
Aie = 700.
Aii = 400.

#Aie = [200, 400, 600, 800, 1000]
#
NNNN = [2000, 4000, 8000, 16000]
# N = 2000
for N in NNNN
IFRAC = 2.
Ni = Int64(round(N/IFRAC))
Ne = N - Ni

s_strength = 3.08
p = 0.1
p0 = sqrt(2000)*p/sqrt(N)

# sn = sqrt(N)
# Aee *= sn
# Aei *= sn
# Aie *= sn
# Aii *= sn

vth = 20
tau_m = 20.
tau_s = 2.

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

W = local_random_2x2(N, IFRAC, p, Aee, Aei, Aie, Aii);
# CSR = sparse_rep(W, N);

@time te, re, ti, ri, SEE, SEI, SIE, SII = euler_lif_2x2_s(h, runtime, N, IFRAC, W, s_strength, vth, tau_m, tau_s)

e_neurons = Neuron_finder(re, 5, min_e_neurons)
i_neurons = Neuron_finder(ri, 5, min_e_neurons)

CV_E = CV_ISI_ALLTIME(e_neurons, te, re)
CV_I = CV_ISI_ALLTIME(i_neurons, ti, ri)

E_R = [length(find(re .== i))/rt for i=1:Ne]
I_R = [length(find(ri .== i))/rt for i=1:Ni]

MER = mean(E_R)#*(1/1000.)
MIR = mean(I_R)#*(1/1000.)

println("##RESULT $(mean(E_R)), $(mean(I_R)), $(mean(CV_E)), $(mean(CV_I))")

WEE_1, WIE_1, WEI_1, WII_1, FE, FI = sim_2_theory_2x2(SEE, SEI, SIE, SII, s_strength, 0, 1., MER, MIR, 100)

println("##PARAMETERS $(WEE_1), $(WIE_1), $(WEI_1), $(WII_1), $(FE), $(FI)")

RE_THEORY, RI_THEORY = theory_rates_2x2(abs(WEE_1), abs(WIE_1), abs(WEI_1), abs(WII_1), FE, FI)

push!(RET, RE_THEORY)
push!(RIT, RI_THEORY)
push!(RES, MER)
push!(RIS, MIR)
push!(ECV, mean(CV_E))
push!(ICV, mean(CV_I))

end

# Aie = [700, 800, 900, 1000]

#RET .*= 1000.
#RIT .*= 1000.
#RES .*= 1000.
#RIS .*= 1000.

NNNN = [2000, 4000, 8000, 16000]
EDIFF = (RES .- RET)
IDIFF = (RIS .- RIT)
EDA = [abs(i) for i in EDIFF]
IDA = [abs(i) for i in IDIFF]


plot(log(NNNN), log(EDA), ".", ms = 20., label = "E")
plot(log(NNNN), log(IDA), ".", ms = 20., label = "I")


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
# Aie = [700, 800, 900, 1000]
# subplot(211)
# plot(Aie, RET, ".", ms = 20., label = "THEORY E CELLS")
# plot(Aie, RES, ".", ms = 20., label = "SIM E CELLS")
# legend()
# ylabel("Mean Firing Rate")
# subplot(212)
# plot(Aie, RIT, ".", ms = 20., label = "THEORY I CELLS")
# plot(Aie, RIS, ".", ms = 20., label = "SIM I CELLS")
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
