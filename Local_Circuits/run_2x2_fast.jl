include("2x2_fast.jl")
include("c://Users//cohenbp//Documents//Neuroscience//Balance_Circuits//two_pools//Analyze.jl")


RET = []
RIT = []

RES = []
RIS = []

# for i in [200, 400, 600, 800, 1000, 1200]
srand(4321)

Aee = 12.5
Aie = 20#i
Aei = 50.
Aii = 50.

N = 5000
IFRAC = 2.
Ni = Int64(round(N/IFRAC))
Ne = N - Ni

k = i#800
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
figure(2)
bar(1:16, s_inf, tick_label = input_strings)
xticks(rotation = 45, fontsize = 24)
yticks(fontsize = 24)
title("Characterizing Inputs in Brent's Network", fontsize = 36)


s_inf[2] + s_inf[3] + s_inf[10] + s_inf[11] + fe2 #input to E pool 2 (dominant)
s_inf[1] + s_inf[4] + s_inf[9] + s_inf[12] + fe1 #input to E pool 1 (suppressed)

s_inf[6] + s_inf[7] + s_inf[14] + s_inf[15] + fi2#input to I pool 2 (dominant)
s_inf[5] + s_inf[8] + s_inf[13] + s_inf[16] + fi1#input to I pool 1 (suppressed)

WEE1 = (s_inf[1] + s_inf[3])/MER1 #output from pool1 excitatory neurons to other excitatory neurons
WEE2 = (s_inf[2] + s_inf[4])/MER2 #same, but pool 2
WIE1 = (s_inf[5] + s_inf[7])/MER1 #from E to I
WIE2 = (s_inf[6] + s_inf[8])/MER2 #from E to I
WEI1 = (s_inf[9] + s_inf[11])/MIR1 #I to E
WEI2 = (s_inf[10] + s_inf[12])/MIR2 #I to E
WII1 = (s_inf[13] + s_inf[15])/MIR1 #I to I
WII2 = (s_inf[14] + s_inf[16])/MIR2 #I to E

#when inputs are the same
s_2x2 = s_inputs ./ (Ne * runtime)
WEE2x = sum(s_2x2[1:4])/ME
WIE2x = sum(s_2x2[5:8])/ME
WEI2x = sum(s_2x2[9:12])/MI
WII2x = sum(s_2x2[13:16])/MI

##OK, these all checked out during a normaliation sim
##and a WTA sim

##now we can start predicting firing rates
# RE, RI = theory_rates_2x2(abs(WEE2x), abs(WIE2x), abs(WEI2x), abs(WII2x), fe1, fi1)

# push!(RET, RE)
# push!(RIT, RI)
# push!(RES, ME)
# push!(RIS, MI)
# end

###OK, now try getting the empirical 2x2.
#for a first try, fe += sum(EE_long, EI_long), fi += sum(IE_long, II_long)
#OK, maybe not now that I think about it because the input is net negative, but then you're just going to use a different threshold to cheat anyway
#in that case, it doesn't matter what f is, you're just going to find threshold such that one guy is off while the other is on
#might still be worth showing what happens

#calculate net excitatory input to both E guys
EE1 = s_inputs[1] + s_inputs[4]
EE2 = s_inputs[2] + s_inputs[3]
#inhibitory to excitatory
EI1 = s_inputs[9] + s_inputs[12]
EI2 = s_inputs[10] + s_inputs[11]
#excitatory to inhibitory
IE1 = s_inputs[5] + s_inputs[8]
IE2 = s_inputs[6] + s_inputs[7]
#inhibitory to inhibitory
II1 = s_inputs[13] + s_inputs[16]
II2 = s_inputs[14] + s_inputs[15]
###ok so the net inputs are identical and the feedforward inputs are the only things that differ
#everything dominant person gives to themselves, they give to suppressed pool too
#Input = r * W, solving for W we get W = Input/r

#calculate just long range inputs
EL12 = (s_inputs[4] + s_inputs[12]) / (Ne2 * runtime)
EL21 = (s_inputs[3] + s_inputs[11]) / (Ne2 * runtime)
IL12 = (s_inputs[8] + s_inputs[16]) / (Ni2 * runtime)
IL21 = (s_inputs[7] + s_inputs[15]) / (Ni2 * runtime)


# TE = [0, .2, .4, .6, .8, 1.]
# for i = 1:length(TE)
RE1, RI1 = theory_rates_2x2(abs(WEE1), abs(WIE1), abs(WEI1), abs(WII1), fe1 - EL1 - TE[i], fi1 - IL1 - TE[i])
RE2, RI2 = theory_rates_2x2(abs(WEE2), abs(WIE2), abs(WEI2), abs(WII2), fe2 - EL2 - TE[i], fi2 - IL2 - TE[i])
# plot(RE1, MER1, ".", label = "E1")
# plot(RE2, MER2, ".", label = "E1")
# plot(RI1, MIR1, ".", label = "E1")
# plot(RI2, MIR2, ".", label = "E1")
# end
# legend()

TE = .3
TI = .6
RE1, RI1 = theory_rates_2x2(abs(WEE1), abs(WIE1), abs(WEI1), abs(WII1), fe1 + EL12 - TE, fi1 + IL12 - TI)
RE2, RI2 = theory_rates_2x2(abs(WEE2), abs(WIE2), abs(WEI2), abs(WII2), fe2 + EL21 - TE, fi2 + IL21 - TI)


# RE, RI = theory_rates_2x2(abs(WEE2x), abs(WIE2x), abs(WEI2x), abs(WII2x), fe1, fi1)



#
#
# # plot([0,3], [0,3], "r")
# # plot(RET, RES, ".", ms = 25., label = "Excitatory Neurons")
# # plot(RIT, RIS, ".", ms = 25., label = "Inhibitory Neurons")
# # xlabel("Theory Rate (Hz)", fontsize = 24)
# # ylabel("Sim Rate (Hz)", fontsize = 24)
#
# subplot(211)
# plot([20, 30, 40, 50, 60], RET, ".", ms = 25., label = "E Theory")
# plot([20, 30, 40, 50, 60], RES, ".", ms = 25., label = "E Sim")
# legend()
# xticks([])
# yticks(fontsize = 18)
# ylabel("Rate (Hz)", fontsize = 24)
#
# subplot(212)
# plot([20, 30, 40, 50, 60], RIT, ".", ms = 25., label = "I Theory")
# plot([20, 30, 40, 50, 60], RIS, ".", ms = 25., label = "I Sim")
# legend()
# yticks(fontsize = 18)
# xticks(fontsize = 18)
# ylabel("Rate (Hz)", fontsize = 24)
#
# xlabel("Aie", fontsize = 36)
#
#for i in [200, 400, 600, 800, 1000, 1200]

# subplot(211)
# plot([200, 400, 600, 800, 1000, 1200], RET, ".", ms = 25., label = "E Theory")
# plot([200, 400, 600, 800, 1000, 1200], RES, ".", ms = 25., label = "E Sim")
# legend()
# xticks([])
# yticks(fontsize = 18)
# ylabel("Rate (Hz)", fontsize = 24)
#
# subplot(212)
# plot([200, 400, 600, 800, 1000, 1200], RIT, ".", ms = 25., label = "I Theory")
# plot([200, 400, 600, 800, 1000, 1200], RIS, ".", ms = 25., label = "I Sim")
# legend()
# yticks(fontsize = 18)
# xticks(fontsize = 18)
# ylabel("Rate (Hz)", fontsize = 24)
#
# xlabel("K", fontsize = 36)
# tight_layout()
#
#
# plot([0,25], [0,25], "r")
# plot(RET, RES, ".", ms = 25., label = "Excitatory Neurons")
# plot(RIT, RIS, ".", ms = 25., label = "Inhibitory Neurons")
# xlabel("Theory Rate (Hz)", fontsize = 24)
# ylabel("Sim Rate (Hz)", fontsize = 24)
