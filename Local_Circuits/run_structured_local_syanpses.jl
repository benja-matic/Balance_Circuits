#run structured local synapses
include("structured_local_synapse.jl")

function sd_2_k(sd)
  return 1./(2*pi*((sd)^2))
end

kee = .9#parse(Float64, ARGS[1])
kei = .9#parse(Float64, ARGS[2])
kie_L = .9#parse(Float64, ARGS[3])
kii = .9#parse(Float64, ARGS[4])

Aee = 50#parse(Float64, ARGS[5])
Aei = 300#parse(Float64, ARGS[6])
Aie_L = 800#parse(Float64, ARGS[7])
Aii= 300#parse(Float64, ARGS[8])

kee = sd_2_k(kee)
kei = sd_2_k(kei)
kie_L = sd_2_k(kie_L)
kii = sd_2_k(kii)

sE = 2.
sI = 0.

Ne = 3200
Ni = 800
N = 4000

p = 0.2
pee = .2
pei = .2
pie = .2
pii = .2

h = .1
vth = 20
tau_m = 20.
tau_ee = 2
tau_ei = 2
tau_ie = 2
tau_ii = 2
tau_ae = 550.
tau_ai = 1000.
g_e = 0.
g_i = 0.

Aee /= p
Aei /= p
Aie_L /= p
Aii /= p

runtime = 10000 #ms
total = round(runtime/h) #time points

wee, wei, wie, wii = weights(Ne,Ni,kee,kei,kie_L,kii,Aee,Aei,Aie_L,Aii,pee,pei,pie,pii)
@time te, re, ti, ri, kill_flag = euler_lif(h, total, Ne, Ni, wee, wei, wie, wii, sE, sI, vth, tau_m, tau_ee, tau_ei, tau_ie, tau_ii, tau_ae, tau_ai, g_e, g_i)

bal = (Aie_L/Aii)*(Ne/Ni)
ral = (length(ti)/Ni)/(length(te)/Ne)

println("Predicted Rates Ratio = $(bal)")
println("Observed Rates Ratio = $(ral)")
