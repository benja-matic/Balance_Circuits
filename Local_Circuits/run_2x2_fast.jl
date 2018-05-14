include("2x2_fast.jl")
include("c://Users//cohenbp//Documents//Neuroscience//Balance_Circuits//two_pools//Analyze.jl")

srand(4321)

Aee = 12.5
Aie = 20.
Aei = 50.
Aii = 50.

N = 5000
IFRAC = 2.
Ni = Int64(round(N/IFRAC))
Ne = N - Ni

k = 800
ks = sqrt(k)
k2 = round(Int64, k/2)
Ne2 = round(Int64, Ne/2)
Ni2 = round(Int64, Ni/2)

fe1 = 3.
fi1 = fe1 - .15

fe2 = fe1 #+ .5
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

figure(1)
plot(t,r,".")

s_inputs = zeros(16);
compartments = [1, Ne2, Ne, Ne+Ni2]

function weights_lookup(i, raster, j_inds)
  kx = find(raster .== i)
  if length(ks) > 0
    for j = 1:4
      ws = sum(W[compartments[j]:compartments[j] + Ne2, i])
      output = ws * tau_s * length(kx)
      s_inputs[j_inds[j]] += output
    end
  end
end

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

#
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

s_2x2 = s_inputs ./ (Ne * runtime)

WEE2x = sum(s_2x2[1:4])/ME
WIE2x = sum(s_2x2[5:8])/ME
WEI2x = sum(s_2x2[9:12])/MI
WII2x = sum(s_2x2[13:16])/MI

##OK, these all checked out during a normaliation sim
##and a WTA sim

##now we can start predicting firing rates
RE = theory_rates_2x2(WEE2x, WIE2x, WEI2x, WII2x, fe1, fi1)
