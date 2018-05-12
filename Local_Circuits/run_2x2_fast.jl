include("2x2_fast.jl")
# include("Analyze.jl")

srand(4321)

Aee = 12.5
Aie = 50.
Aei = 20.
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
MER1 = length(er1)/(Ne2*rt);
MER2 = length(er2)/(Ne2*rt);

figure(1)
plot(t,r,".")

###################
#IE means E to I, 12 means from 2 to 1
#EE11, EE22, EE21, EE12,
#IE11, IE22, IE21, IE12,
#EI11, EI22, EI21, EI12,
#II11, II22, II21, II12
s_inputs = zeros(16);
#  sum over A__/sqrt(k) * k_out * tau_s, divide by number of neurons, time

#loop over neurons
#get firing rate
#loop over 4 compartments of matrix (local vs non-local X 2 connection types) ... (e.g. local EE, nonlocal EE, Local IE, nonlocal IE)
#

compartments = [1, Ne2, Ne, Ne+Ni2]


function weights_lookup(i, raster, j_inds)
  kx = find(raster .== i)
  if length(ks) > 0
    for j = 1:4
      # wxr = find(W[compartments[j]:compartments[j] + Ne2, i])
      # A = W[compartments[j]:compartments[j] + Ne2, i][wxr[1]]
      # output = A * tau_s * length(wxr) * length(kx)#integral over synaptic waveform * number_outputs * number of spikes *
      ### fails to account for tupled connections
      ### should really be sum over column compartment * tau_s * number of spikes
      ws = sum(W[compartments[j]:compartments[j] + Ne2, i])
      output = ws * tau_s * length(kx)
      s_inputs[j_inds[j]] += output
    end
  end
end

#

#EE11, EE22, EE21, EE12,
#IE11, IE22, IE21, IE12,
#EI11, EI22, EI21, EI12,
#II11, II22, II21, II12
#compartments = [1, Ne2, Ne, Ne+Ni2]


#this loop assumes Ne=Ni for convenience
for i = 1:Ne2
  weights_lookup(i, re, [1, 3, 5, 7])
  weights_lookup(i+Ne2, re, [2, 4, 6, 8])
  weights_lookup(i + Ne, ri, [9, 11, 13, 15])
  weights_lookup(i+Ne+Ne2, ri, [10, 12, 14, 16])
end

normz = Ne2 * runtime

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


###checking for bugs
#first issue is that you can have repeat indices
#to resolve this we can compute sums over every row and column and compare to what I think we should get
Jee = Aee/ks
Jei = -Aei/ks
Jie = Aie/ks
Jii = -Aii/ks

e_row_t = Jee*k + Jei*k
i_row_t = Jie*k + Jii*k
# e_col_t =

# c_out = zeros(N);
# c_outx = zeros(N);
# c_in = zeros(N);
# c_inx = zeros(N);
#
# for i = 1:N
#   kx = sum(W[:,i])
#   c_out[i] = kx
#   kx2 = sum(W[i,:])
#   c_outx[i] = length(find(W[:,i]))
#   c_in[i] = kx2
#   c_inx[i] = length(find(W[i,:]))
# end
# figure(3)
# plot(c_in)
# axhline(e_row_t, linestyle = "dashed")
# axhline(i_row_t, linestyle = "dashed")
#that's a good lookin' check

#plot(c_out)
#outputs are not constrained, so we expect some heterogeneity
#Each time we decide who supplies inputs, we choose a neuron with probability k2/(Ne2=Ni2) = k/(Ne=Ni)
#This experiment is repeated N times, thus we expect each neuron to have N*k/Ne outputs = 2k
#plot(c_outx) #looks like the number of UNIQUE outgoing connections is a bit less than 2k...approximately correct
#is it really 2k - (N*2k/N^2) - (N*2k/N^3) - ...? whenever you get doubled, someone else gets left out
# #we know past N, the terms are virtually vanishing, just subtract out N^2 terms
# p0 = 2k/N
# p1 = 2k/N^2
# #probability of being chosen is 2k/N for each row, so probability of not being chosen is 1 - 2k/n ~ .68
# figure(4)
# plot(sort(c_inx), sort(c_outx))
# plot([1200,1600], [1200, 1600], "r")
#this plot seems to say that there are a higher number of unique outgoing connections than incoming connections


##try checking some sums over Ne2 sized pieces of columns
ees1 = zeros(Ne2);
ees2 = zeros(Ne2);
eis1 = zeros(Ni2);
eis2 = zeros(Ni2);
eis3 = zeros(Ni2);#just spot check pool 2 E neurons
eis4 = zeros(Ni2);
ies1 = zeros(Ne2);
ies2 = zeros(Ne2);
iis1 = zeros(Ni2);
iis2 = zeros(Ni2);

for i = 1:Ne2
  e1 = sum(W[1:Ne2, i])
  e2 = sum(W[1:Ne2, i+Ne2])
  ees1[i] = e1
  ees2[i] = e2
end

for i = 1:Ne2
  e1 = sum(W[1:Ne2, i+Ne])
  e2 = sum(W[1:Ne2, i+Ne+Ni2])
  eis1[i] = e1
  eis2[i] = e2
end

for i = 1:Ne2
  e1 = sum(W[Ne+1:Ne+Ni2, i])
  e2 = sum(W[Ne+1:Ne+Ni2, i+Ni2])
  ies1[i] = e1
  ies2[i] = e2
end

for i = 1:Ne2
  e1 = sum(W[Ne+1:Ne+Ni2, i+Ne])
  e2 = sum(W[Ne+1:Ne+Ni2, i+Ne+Ni2])
  iis1[i] = e1
  iis2[i] = e2
end

for i = 1:Ne2
  e1 = sum(W[Ne2+1:Ne, i+Ne])
  e2 = sum(W[Ne2+1:Ne, i+Ne+Ni2])
  eis3[i] = e1
  eis4[i] = e2
end

#verified matrix is not biased, issue must be in my measures

#


#c_in should have 0 heterogeneity, should just be 2k every time

# Wx = zeros(N,N);
# e_inds = rand(1:Ne2, k2);
# i_inds = rand(Ne+1:Ne+Ni2, k2);
#
# for j in eachindex(e_inds)
#   Wx[1, e_inds[j]] += 1
#   Wx[1, e_inds[j] + Ne2] += 1
#   Wx[1, i_inds[j]] += 1
#   Wx[1, i_inds[j] + Ni2] += 1
# end
