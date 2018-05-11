include("2x2_fast.jl")
# include("Analyze.jl")

srand(4321)

Aee = 12.5
Aei = 20.
Aie = 50.
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
MER1 = length(er1)/(Ne2*rt);
MER2 = length(er2)/(Ne2*rt);



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


function weights_lookup(i, raster, modu, j_inds)
  kx = find(raster .== i)
  if length(ks) > 0
    for j = 1:4
      wxr = find(W[compartments[j]:compartments[j] + Ne2, i + modu])
      A = W[compartment:compartment + Ne2][wxr[1]]
      output = A * length(wxr) * length(kx)
      s_inputs[j_inds[j]] += output
    end
  end
end

#



#this loop assumes Ne=Ni for convenience
for i = 1:Ne2
  weights_lookup(i, re, 0, [1, 3, ])



i = 2000
kx = find(re .== i);
wx1 = length(find(W[1:Ne2, i]));
wx2 = length(find(W[Ne2+1:Ne, i]));
wx3 = length(find(W[Ne+1:Ne + Ne2, i]));
wx4 = length(find(W[Ne+Ne2:N, i]));




#
