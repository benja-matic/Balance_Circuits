include("simplest_network.jl")

N = 5000
Ni = Int64(round(N/5))
Ne = N - Ni
width = .35
p = .2

Aee = 10.
Aei = 350.
Aie = 150.0
Aie_NL = 250.
Aii = 200.

W = homogenous_weights(N, p, Aee, Aei, Aie, Aie_NL, Aii, width);
CSR = sparse_rep(W, N);

runtime = 10000
rt = runtime /1000.
h = 0.1
total = runtime/h
s1 = s2 = 2.
vth = 20.
tau_m = 20.
tau_s = 2.

@time t, r = Simple_Network_CSR(h, runtime, CSR, W, N, s1, s2, vth, tau_m, tau_s);

ref = find(r .<= Ne)
rif = find(r .> Ne)

re = r[ref]
te = t[ref]
ri = r[rif]
ti = t[rif]
tem = find(te.> end_trans)
te_pt = te[tem]
re_pt = re[tem]
tim = find(ti.> end_trans)
ti_pt = ti[tim]
ri_pt = ri[tim]

min_e_neurons = 200
min_i_neurons = 50
runtime = 14000 #ms
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
end_trans = 40000
rt = ((ntotal - end_trans)/1000.)*h
half_e = div(Ne, 2)
half_i = div(Ni, 2)

wta_ness, bias = score_analysis(re_pt, Ne)
top_e_neurons, bot_e_neurons = Neurons_tb_ns(re_pt, half_e, 10, min_e_neurons)#22
if ((top_e_neurons == -5) & (bot_e_neurons == -5))
  println("##RESULT -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, $(kee), $(kei), $(kie), $(kii), $(Aee*p), $(Aei*p), $(Aie*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
  #Set stuff equal to error code
  #println("##RESULT stuff")
else
I_Neurons = Neuron_finder(ri_pt, 10, min_i_neurons)
#rates
E_Rates = [length(find(re .== i))/rt for i=1:Ne]
I_Rates = [length(find(ri .== i))/rt for i=1:Ni]
E_rate = (length(re_pt)/rt)/Ne
I_rate = (length(ri_pt)/rt)/Ni

#counts
E_count_top = count_train(fbinsize, te_pt, re_pt, top_e_neurons, length(top_e_neurons))
E_count_bot = count_train(fbinsize, te_pt, re_pt, bot_e_neurons, length(bot_e_neurons))
I_count_all = count_train(fbinsize, ti_pt, ri_pt, I_Neurons, length(I_Neurons))
#FF
ntfa = network_fano(te_pt, re_pt, cbinsize, ntotal)
E_FANO_mean_top, E_FANO_median_top, E_FANO_std_top = fano_train(E_count_top, -5)
E_FANO_mean_bot, E_FANO_median_bot, E_FANO_std_bot = fano_train(E_count_bot, -5)
I_FANO_mean, I_FANO_median, I_FANO_std = fano_train(I_count_all, -5)
#cv
E_CV_mean_top, E_CV_median_top, E_CV_STD_top = CV_ISI_ALLTIME(top_e_neurons, te_pt, re_pt)
E_CV_mean_bot, E_CV_median_bot, E_CV_STD_bot = CV_ISI_ALLTIME(bot_e_neurons, te_pt, re_pt)
I_CV_mean, I_CV_median, I_CV_STD = CV_ISI_ALLTIME(I_Neurons, ti_pt, ri_pt)
#synchrony
E_spike_correlation_top = rand_pair_cor(cbinsize, te_pt, re_pt, top_e_neurons, 1000)
E_spike_correlation_bot = rand_pair_cor(cbinsize, te_pt, re_pt, bot_e_neurons, 1000)
I_spike_correlation = rand_pair_cor(cbinsize, ti_pt, ri_pt, I_Neurons, 500)
#inputs
# etm, etstd, etbal, etskew, etkurt = moments(e_top_pt)
# ebm, ebstd, ebbal, ebskew, ebkurt = moments(e_bot_pt)
# igm, igstd, igbal, igskew, igkurt = moments(i_guy_pt)
#normality analysis
println("##RESULT $(E_rate), $(I_rate), $(wta_ness), $(E_FANO_mean_top), $(E_FANO_mean_bot), $(I_FANO_mean), $(E_CV_mean_top), $(E_CV_mean_bot), $(I_CV_mean)")
# println("##RESULT $(E_rate), $(I_rate), $(wta_ness), $(bias), $(ntfa), $(E_FANO_mean_top), $(E_FANO_mean_bot), $(I_FANO_mean), $(E_CV_mean_top), $(E_CV_mean_bot), $(I_CV_mean), $(mean(E_spike_correlation_top)), $(mean(E_spike_correlation_bot)), $(I_spike_correlation)")
end


figure(1)
plot(te, re, "g.")
figure(2)
plot(ti, ri, "b.")
figure(3)
a=imshow(W)
colorbar(a)
