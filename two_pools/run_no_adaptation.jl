#println("GOT THIS FAR")
#include("measure_balance.jl")
include("euler_WTA.jl")

srand(1234)
function sd_2_k(sd)
  return 1./(2*pi*((sd)^2))
end

kee = .26#parse(Float64, ARGS[1])
kei = .93#parse(Float64, ARGS[2])
kie = .97#parse(Float64, ARGS[3])
kii = .5#parse(Float64, ARGS[4])

Aee = 84#parse(Float64, ARGS[5])
Aei = 314#parse(Float64, ARGS[6])
Aie = 1319#parse(Float64, ARGS[7])
Aii= 689#parse(Float64, ARGS[8])

s_strength = 3.08#parse(Float64, ARGS[9])
p = .34#parse(Float64, ARGS[10])
Ne = 3200#parse(Int64, ARGS[11])
Ni = div(Ne, 4)
h = .1
vth = 20
tau_m = 20.
tau_ee = 2
tau_ei = 2
tau_ie = 2
tau_ii = 2


kee = sd_2_k(kee)
kei = sd_2_k(kei)
kie = sd_2_k(kie)
kii = sd_2_k(kii)

# Aee = 0.#80.19552235
# Aei = 366.75109066
# Aie_L = 648.52219365
# Aii = 400.#130.916487007

# p = 0.245403843441
# s_strength = 2.84452373606
s1 = s_strength
s2 = s_strength

# Ne = 3200
# Ni = 800
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

Aee /= p
Aei /= p
Aie /= p
Aii /= p
#tic()
c=0


wee,wei,wie,wii = weights(Ne,Ni,kee, kei, kie, kii, Aee, Aei, Aie, Aii, p,p,p,p);
# for gamma = .0001:0.005:.05
te,re,ti,ri,kill_flag, e_top, e_bot, i_guy = euler_lif(h,runtime,Ne,Ni,wee,wei,wie,wii,s1,s2, vth, tau_m, tau_ee, tau_ei, tau_ie, tau_ii);

###Kill transient
e_top_pt = e_top[end_trans:end]
e_bot_pt = e_bot[end_trans:end]
i_guy_pt = i_guy[end_trans:end]

tem = find(te.> end_trans)
te_pt = te[tem]
re_pt = re[tem]
tim = find(ti.> end_trans)
ti_pt = ti[tim]
ri_pt = ri[tim]

wta_ness, bias = score_analysis(re_pt, Ne)
top_e_neurons, bot_e_neurons = Neurons_tb_ns(re_pt, half_e, 10, min_e_neurons)#22
if ((top_e_neurons == -5) & (bot_e_neurons == -5))
  println("##RESULT -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, $(kee), $(kei), $(kie), $(kii), $(Aee*p), $(Aei*p), $(Aie*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
  #Set stuff equal to error code
  #println("##RESULT stuff")
else
I_Neurons = Neuron_finder(ri_pt, 10, min_i_neurons)
#rates
E_rate = (length(re_pt)/rt)/Ne
I_rate = (length(ri_pt)rt)/Ni

#counts
E_count_top = count_train(fbinsize, te_pt, re_pt, top_e_neurons, length(top_e_neurons))
E_count_bot = count_train(fbinsize, te_pt, re_pt, bot_e_neurons, length(bot_e_neurons))
I_count_all = count_train(fbinsize, ti_pt, ri_pt, I_Neurons, length(I_Neurons))
#FF
ntf = network_fano(te_pt, re_pt, cbinsize)
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
etm, etstd, etbal, etskew, etkurt = moments(e_top_pt)
ebm, ebstd, ebbal, ebskew, ebkurt = moments(e_bot_pt)
igm, igstd, igbal, igskew, igkurt = moments(i_guy_pt)
#normality analysis

println("##RESULT $(E_rate), $(I_rate), $(wta_ness), $(bias), $(ntf), $(E_FANO_mean_top), $(E_FANO_mean_bot), $(I_FANO_mean), $(E_CV_mean_top), $(E_CV_mean_bot), $(I_CV_mean), $(E_spike_correlation_top), $(E_spike_correlation_bot), $(I_spike_correlation), $(etm), $(etstd), $(etbal), $(etskew), $(etkurt), $(ebm), $(ebstd), $(ebbal), $(ebskew), $(ebkurt), $(igm), $(igstd), $(igbal), $(igskew), $(igkurt), $(kee), $(kei), $(kie), $(kii), $(Aee*p), $(Aei*p), $(Aie*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
end
