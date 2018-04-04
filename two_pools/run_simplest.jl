include("simplest_network_in.jl")
include("Analyze.jl")
N = 5000
Ni = Int64(round(N/5))
Ne = N - Ni
width = .35
p = .2

Aee = parse(Float64, ARGS[1])#10.
Aei = 500.
Aie = parse(Float64, ARGS[2])#2050.0
Aie_NL = parse(Float64, ARGS[3])#0.#150.
Aii = 500.
# Aie = [150, 550, 1050, 1550, 2050]
# Aie_NL = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1500, 3000]

W = homogenous_weights(N, p, Aee, Aei, Aie, Aie_NL, Aii, width);
CSR = sparse_rep(W, N);

runtime = 17000
rt = runtime /1000.
end_trans = 2000
h = 0.1
total = runtime/h
s1 = s2 = 2.
vth = 20.
tau_m = 20.
tau_s = 2.

min_e_neurons = 200
min_i_neurons = 50
ntotal = round(runtime/h) #time points
fbinsize = 400/h
cbinsize = 100/h
# end_trans = 40000
rt = ((ntotal - end_trans)/1000.)*h
half_e = div(Ne, 2)
half_i = div(Ni, 2)

@time t, r, Input = Simple_Network_CSR(h, runtime, CSR, W, N, s1, s2, vth, tau_m, tau_s);

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

wta_ness, bias = score_analysis(re_pt, Ne)
top_e_neurons, bot_e_neurons = Neurons_tb_ns(re_pt, half_e, 10, min_e_neurons)#22
I_Neurons = Neuron_finder(ri_pt, 10, min_i_neurons)
#rates
E_R_bot = [length(find(re_pt .== i))/rt for i=1:half_e]
E_R_top = [length(find(re_pt .== i))/rt for i=half_e+1:Ne]
I_R = [length(find(ri_pt .== i))/rt for i=Ne+1:N]
E_R = [length(find(re_pt .== i))/rt for i=1:Ne]

#cv
CV_TOP = CV_ISI_ALLTIME(top_e_neurons, te_pt, re_pt)
CV_BOT = CV_ISI_ALLTIME(bot_e_neurons, te_pt, re_pt)
CV_INH = CV_ISI_ALLTIME(I_Neurons, ti_pt, ri_pt)


println("##RESULT $(wta_ness), $(bias), $(E_R_bot), $(E_R_top), $(I_R), $(E_R), $(mean(CV_TOP)), $(mean(CV_BOT)), $(mean(CV_INH)), $(mean(Input[1,:][:])), $(mean(Input[2,:][:])), $(mean(Input[3,:][:])), $(mean(Input[4,:][:]))")
# figure(1)
# plot(te, re, "g.", ms = 1.)
#
# push!(IR, mean(I_R))
# push!(ECVT, mean(CV_TOP))
# push!(ECVB, mean(CV_BOT))
# push!(ER, mean(E_R))
# figure(2)
# plot(ti, ri, "b.")
# figure(3)
# a=imshow(W)
# colorbar(a)

# Aie_NL = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1500, 3000]
# subplot(211)
# plot(Aie_NL, IR, ".")
# axvline(400, linestyle = "dashed")
# ylabel("Firing Rate")
# title("Inhibitory Firing Rate")
# subplot(212)
# title("Excitatory CV_ISI")
# plot(Aie_NL, ECVT, ".", label = "E Pool 1")
# plot(Aie_NL, ECVB, ".", label = "E Pool 2")
# ylabel("CV_ISI")
# axvline(400, linestyle = "dashed")
# legend()
# Aie = [150, 300, 450, 600, 750, 900 ]
#
# Aie = [150, 550, 1050, 1550, 2050]
# subplot(211)
# plot(Aie, IR, ".", label = "Inhibitory")
# plot(Aie, ER, ".", label = "Excitatory")
# ylabel("Rate")
# legend()
# subplot(212)
# plot(Aie, ECVT, ".", label = "Excitatory Pool 1")
# plot(Aie, ECVB, ".", label = "Excitatory Pool 2")
# ylabel("CV_ISI")
# xlabel("Aie LOCAL")
# legend()
