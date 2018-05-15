include("2x2_fast.jl")
include("c://Users//cohenbp//Documents//Neuroscience//Balance_Circuits//two_pools//Analyze.jl")


# RET1 = []
# RIT1 = []
# RET2 = []
# RIT2 = []
#
# RES1 = []
# RIS1 = []
# RES2 = []
# RIS2 = []

# for i in [4000, 8000, 16000, 32000]
srand(4321)

Aee = 12.5
Aie = 20#i
Aei = 50.
Aii = 50.

N = 32000#i#5000
IFRAC = 2.
Ni = Int64(round(N/IFRAC))
Ne = N - Ni

k = 800
ks = sqrt(k)
k2 = round(Int64, k/2)
Ne2 = round(Int64, Ne/2)
Ni2 = round(Int64, Ni/2)

fe1 = 3
fi1 = fe1 - .15

fe2 = fe1 + .5
fi2 = fe2 - .15

vth = 20
tau_m = 20.
tau_s = 2.

min_e_neurons = 20
min_i_neurons = 50
runtime = 10000 #ms
h = .1
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
netd_binsize = 50/h
end_trans = 0.
# rt = ((ntotal - end_trans)/1000.)*h
rt = runtime/1000
W = local_random_2x2_symmetric(N, IFRAC, k, Aee, Aei, Aie, Aii);
CSR = sparse_rep(W, N);

@time t, r = euler_lif_2x2_CSR(h, runtime, N, IFRAC, W, CSR, fe1, fi1, fe2, fi2, vth, tau_m, tau_s)

an = Set(1:length(r));
em = find(r .<= Ne);
im = collect(setdiff(an, Set(em)));
te = t[em];
re = r[em];

ti = t[im];
ri = r[im];

er1 = Set(find(re .< Ne2));
er2 = setdiff(Set(1:length(re)), er1);
ir1 = Set(find(ri .< Ne + Ni2));
ir2 = setdiff(Set(1:length(ri)), ir1);

ME = length(re)/(Ne*rt);
MI = length(ri)/(Ni*rt);
MER1 = length(er1)/(Ne2*rt);
MER2 = length(er2)/(Ne2*rt);
MIR1 = length(ir1)/(Ni2*rt);
MIR2 = length(ir2)/(Ni2*rt);

# figure(1)
# plot(t,r,".")

s_inputs = zeros(16);
compartments = [1, Ne2, Ne, Ne+Ni2]

#this loop assumes Ne=Ni for convenience
for i = 1:Ne2
  weights_lookup(i, re, [1, 3, 5, 7])
  weights_lookup(i+Ne2, re, [2, 4, 6, 8])
  weights_lookup(i + Ne, ri, [9, 11, 13, 15])
  weights_lookup(i+Ne+Ne2, ri, [10, 12, 14, 16])
end

normz = Ne2 * runtime #note that when you pretend this is 2x2, if you sum over all 4 inputs of a given type, you will get double what you expect since you normalized by Ne2 here, not Ne

input_strings = ["EE11", "EE22", "EE21", "EE12", "IE11", "IE22", "IE21", "IE12", "EI11", "EI22", "EI21", "EI12", "II11", "II22", "II21", "II12",]

s_inf  = s_inputs ./ normz
# figure(2)
# bar(1:16, s_inf, tick_label = input_strings)
# xticks(rotation = 45, fontsize = 24)
# yticks(fontsize = 24)
# title("Characterizing Inputs in Brent's Network", fontsize = 36)

#1.
WEE1 = (s_inf[1])/MER1 #output from pool1 excitatory neurons to other excitatory neurons
WEE2 = (s_inf[2])/MER2 #same, but pool 2
WIE1 = (s_inf[5])/MER1 #from E to I
WIE2 = (s_inf[6])/MER2 #from E to I
WEI1 = (s_inf[9])/MIR1 #I to E
WEI2 = (s_inf[10])/MIR2 #I to E
WII1 = (s_inf[13])/MIR1 #I to I
WII2 = (s_inf[14])/MIR2 #I to E

#2.
EL12 = (s_inputs[4] + s_inputs[12]) / (Ne2 * runtime)
EL21 = (s_inputs[3] + s_inputs[11]) / (Ne2 * runtime)
IL12 = (s_inputs[8] + s_inputs[16]) / (Ni2 * runtime)
IL21 = (s_inputs[7] + s_inputs[15]) / (Ni2 * runtime)

#3': just calculate fe, let it be negative if it's negative
fe_1 = fe1 + EL12
fi_1 = fi1 + IL12

fe_2 = fe2 + EL21
fi_2 = fi2 + IL21
#OK, turns out it's negative for one pool and positive for the other...which is what we'd expect

RE1T, RI1T = theory_rates_2x2(abs(WEE1), abs(WEI1), abs(WIE1), abs(WII1), fe_1, fi_1)
RE2T, RI2T = theory_rates_2x2(abs(WEE2), abs(WEI2), abs(WIE2), abs(WII1), fe_2, fi_2)

# push!(RET1, RE1T)
# push!(RIT1, RI1T)
# push!(RET2, RE2T)
# push!(RIT2, RI2T)
#
# push!(RES1, MER1)
# push!(RES2, MER2)
# push!(RIS1, MIR1)
# push!(RIS2, MIR2)
# end
#
#
# plot([4000, 8000, 16000, 32000], RES2, ".", ms = 10., label = "E2 simulation")
# plot([4000, 8000, 16000, 32000], RET2, ".", ms = 10., label = "E2 theory")
#
# plot([4000, 8000, 16000, 32000], RIS2, ".", ms = 10., label = "I2 simulation")
# plot([4000, 8000, 16000, 32000], RIT2, ".", ms = 10., label = "I2 theory")
# legend()
# xlabel("N", fontsize = 24)
# ylabel("Rate (Hz)")
# title("Dominant Pool Firing Rate vs. Empirical 2x2 Theory")
#
# plot([4000, 8000, 16000, 32000], RES1, ".", ms = 10., label = "E1 simulation")
# plot([4000, 8000, 16000, 32000], RET1, ".", ms = 10., label = "E1 theory")
#
# plot([4000, 8000, 16000, 32000], RIS1, ".", ms = 10., label = "I1 simulation")
# plot([4000, 8000, 16000, 32000], RIT1, ".", ms = 10., label = "I1 theory")
# xlabel("N", fontsize = 24)
# ylabel("Rate (Hz)")
# title("Suppressed Pool Firing Rate vs. Empirical 2x2 Theory")
