#Analyze sims with julia
using DataFrames
data = readtable("DFS_results.csv");
data2 = readtable("DC_results.csv");
N32 = readtable("N32K_f01_results.csv");
N64 = readtable("N64K_f01_results.csv");

# cd("c://Users\\cohenbp\\Documents\\Neuroscience\\Balance_Circuits\\two_pools\\")
cd("c://Users\\cohenbp\\Documents\\Neuroscience\\data\\theory")
include("c://Users\\cohenbp\\Documents\\Neuroscience\\Balance_Circuits\\two_pools\\Analyze.jl")


#println("##RESULT $(wta_ness), $(bias), $(mean(E_R_bot)), $(mean(E_R_top)), $(mean(I_R_top)), $(mean(I_R_bot)), $(mean(CV_ETOP)), $(mean(CV_EBOT)), $(mean(CV_ITOP)), $(mean(CV_IBOT)), $(WEE_1), $(WEE_2), $(WIE_1), $(WIE_2), $(WIEL_1), $(WIEL_2), $(WEI_1), $(WEI_2), $(WII_1), $(WII_2), $(FE1), $(FE2), $(Aee), $(Aei), $(Aie), $(Aie_NL), $(Aii), $(fe), $(fe2), $(fi), $(fi2), $(delcon)")


#RE1_4x4_sf_dense, RE2_4x4_sf_dense, RI1_4x4_sf_dense, RI2_4x4_sf_dense = theory_rates_4x4_sf(WEEzd, WIEzd, WIELzd, WEIzd, WIIzd, fe .- TE1d, 0 .- TI1d, fe +.01 .- TE2d, 0 .- TI2d)

win_e_rate = zeros(99)
los_e_rate = zeros(99)
for i = 1:99
    if data[i,3] > data[i,4]
        win_e_rate[i] = data[i,3]
        los_e_rate[i] = data[i,4]
    else
        win_e_rate[i] = data[i,4]
        los_e_rate[i] = data[i,3]
    end
end

WEE = [abs(data[i, 11]) for i = 1:99]
WIE = [abs(data[i, 13]) for i = 1:99]
WIEL = [abs(data[i, 15]) for i = 1:99]
WEI = [abs(data[i, 17]) for i = 1:99]
WII = [abs(data[i, 19]) for i = 1:99]

RE1_4x4_sf_dense, RE2_4x4_sf_dense, RI1_4x4_sf_dense, RI2_4x4_sf_dense = theory_rates_4x4_sf(WEE, WIE, WIEL, WEI, WII, -1. + data[:,28] .- -.5, data[:,29][:] .- -.12, -1 .+ data[:,30][:] .- -.5, data[:,31][:] .- -.12)

plot(data[:,end-6][:], RE2_4x4_sf_dense, ".", ms = 5., label = "SF 4x4 Exc Prediction")
plot(data[:,end-6][:], win_e_rate, ".", ms = 5., label = "winner firing rate")
plot(data[:,end-6][:], los_e_rate, ".", ms = 5., label = "loser firing rate")
legend()
xlabel("Long Range Strength")
ylabel("Firing Rate (Hz)")
title("4x4 SF Theory vs. Sim")


RE1_4x4_sf_dense2, RE2_4x4_sf_dense2, RI1_4x4_sf_dense2, RI2_4x4_sf_dense2 = theory_rates_4x4_sf(data2[:,11][:], data2[:,13][:], data2[:,15][:], data2[:,17][:], data2[:,19][:], data2[:,28] .- -.5, data2[:,29][:] .- -.12, data2[:,30][:] .- -.5, data2[:,31][:] .- -.12)

RE1_4x4_sf_dense2b, RE2_4x4_sf_dense2b, RI1_4x4_sf_dense2b, RI2_4x4_sf_dense2b = theory_rates_4x4_sf(data2[:,11][:], data2[:,13][:], data2[:,15][:], data2[:,18][:], data2[:,19][:], data2[:,28] .- -.5, data2[:,29][:] .- -.12, data2[:,30][:] .- -.5, data2[:,31][:] .- -.12)

win_e_rate2 = max(data[:,3][:], data[:,4][:])
los_e_rate2 = min(data[:,3][:], data[:,4][:])

plot(data2[:,end-6][:], RE1_4x4_sf_dense2, ".", ms = 5., label = "SF 4x4 Exc Prediction")
plot(data2[:,end-6][:], RE2_4x4_sf_dense2, ".", ms = 5., label = "SF 4x4 Exc Prediction")

plot(data2[:,end-6][:], RE1_4x4_sf_dense2b, ".", ms = 5., label = "SF 4x4 Exc Prediction")
plot(data2[:,end-6][:], RE2_4x4_sf_dense2b, ".", ms = 5., label = "SF 4x4 Exc Prediction")

plot(data2[:,end-6][:], win_e_rate2, ".", ms = 5., label = "winner firing rate")
plot(data2[:,end-6][:], los_e_rate2, ".", ms = 5., label = "loser firing rate")
legend()
xlabel("Long Range Strength")
ylabel("Firing Rate (Hz)")
title("4x4 SF Theory vs. Sim")




#println("##RESULT $(wta_ness), $(bias), $(mean(E_R_bot)), $(mean(E_R_top)), $(mean(I_R_top)), $(mean(I_R_bot)), $(mean(CV_ETOP)), $(mean(CV_EBOT)), $(mean(CV_ITOP)), $(mean(CV_IBOT)), $(WEE_1), $(WEE_2), $(WIE_1), $(WIE_2), $(WIEL_1), $(WIEL_2), $(WEI_1), $(WEI_2), $(WII_1), $(WII_2), $(FE1), $(FE2), $(Aee), $(Aei), $(Aie), $(Aie_NL), $(Aii), $(fe), $(fe2), $(fi), $(fi2), $(delcon)")


#we know A_ij *sqrt(k)*tau_s/1000. ~ theoretical W (slightly above due to lx correction)
#you have to divide Wie and WieL by two since E to I is half WieL and half Wie

function s2ta(A, k, tau_s)
    return A*sqrt(k)*tau_s/1000.
end
#$(Aee), $(Aei), $(Aie), $(Aie_NL), $(Aii)
Wee=s2ta(data[:,end-9][:], 800, 2)
Wei=s2ta(data[:,end-8][:], 800, 2)
Wie=s2ta(data[:,end-7][:], 800, 2)/2.
WieL=s2ta(data[:,end-6][:], 800, 2)/2.
Wii=s2ta(data[:,end-5][:], 800, 2)
#WEE_1, WIE_1, WIEL_1, WEI_1, WII_1, FE1, FI1, FE2, FI2

RE1_4x4_sf_dense, RE2_4x4_sf_dense, RI1_4x4_sf_dense, RI2_4x4_sf_dense = theory_rates_4x4_sf(Wee, Wie, WieL, Wei, Wii, 3.08, 0., 3.09, 0.)


# plot(data[:, end-6], RE1_4x4_sf_dense, ".", ms = 8., label = "4x4sf theory prediction")
plot(data[:, end-6], RE2_4x4_sf_dense, ".", ms = 5., label = "4x4sf theory prediction")
plot(data[:,end-6][:], win_e_rate, ".", ms = 5., label = "winner firing rate")
plot(data[:,end-6][:], los_e_rate, ".", ms = 5., label = "loser firing rate")
legend()
xlabel("Long Range Strength")
ylabel("Firing Rate (Hz)")
title("4x4 SF Theory vs. Sim")



#
Wee_N32=s2ta(N32[:,end-9][:], 800, 2)
Wei_N32=s2ta(N32[:,end-8][:], 800, 2)
Wie_N32=s2ta(N32[:,end-7][:], 800, 2)/2.
WieL_N32=s2ta(N32[:,end-6][:], 800, 2)/2.
Wii_N32=s2ta(N32[:,end-5][:], 800, 2)

RE1_4x4_sf_dense_N32, RE2_4x4_sf_dense_N32, RI1_4x4_sf_dense_N32, RI2_4x4_sf_dense_N32 = theory_rates_4x4_sf(Wee_N32, Wie_N32, WieL_N32, Wei_N32, Wii_N32, 3.08, 0., 3.09, 0.)

Wee_N64=s2ta(N64[:,end-9][:], 800, 2)
Wei_N64=s2ta(N64[:,end-8][:], 800, 2)
Wie_N64=s2ta(N64[:,end-7][:], 800, 2)/2.
WieL_N64=s2ta(N64[:,end-6][:], 800, 2)/2.
Wii_N64=s2ta(N64[:,end-5][:], 800, 2)

RE1_4x4_sf_dense_N64, RE2_4x4_sf_dense_N64, RI1_4x4_sf_dense_N64, RI2_4x4_sf_dense_N64 = theory_rates_4x4_sf(Wee_N64, Wie_N64, WieL_N64, Wei_N64, Wii_N64, 3.08, 0., 3.09, 0.)

win_e_n32 = max(N32[:,3], N32[:,4]);
los_e_n32 = min(N32[:,3], N32[:,4]);

win_e_n64 = max(N64[:,3], N64[:,4]);
los_e_n64 = min(N64[:,3], N64[:,4]);

#verified that the predictions are always the same in the theory
plot(data[:, end-6], RE2_4x4_sf_dense, ".", ms = 5., label = "4x4sf theory prediction")
plot(data[:, end-6], RE2_4x4_sf_dense_N32, ".", ms = 5., label = "4x4sf theory prediction")
plot(data[:, end-6], RE2_4x4_sf_dense_N64, ".", ms = 5., label = "4x4sf theory prediction")

plot(data[:, end-6], RE2_4x4_sf_dense, ".", ms = 5., label = "4x4sf theory prediction")
plot(data[:,end-6][:], win_e_rate, ".", ms = 5., label = "N=16k")
plot(N32[:, end-6], win_e_n32, ".", ms = 5., label = "N=32k")
plot(N64[:, end-6], win_e_n64, ".", ms = 5., label = "N=64k")
legend()

#
