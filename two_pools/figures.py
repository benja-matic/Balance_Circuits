import numpy as np
import matplotlib.pyplot as plt

def parse_inputs(f):
    with open(f,'r') as x:
        a = x.read()
    b = a.split('\n')
    c = np.array([float(b[i]) for i in range(len(b)-1)])
    return c

###Normalization Data
NTIM = parse_inputs("N_TOP_IM.txt")
NBIM = parse_inputs("N_BOT_IM.txt")
NTRSC = parse_inputs("N_TOP_RSC.txt")
NBRSC = parse_inputs("N_BOT_RSC.txt")
NTFF = parse_inputs("N_TOP_FF.txt")
NBFF =parse_inputs("N_BOT_FF.txt")
NTCV = parse_inputs("N_TOP_CV.txt")
NBCV = parse_inputs("N_BOT_CV.txt")

###WTA Data
WTIM = parse_inputs("W_TOP_IM.txt")
WBIM = parse_inputs("W_BOT_IM.txt")
WTRSC = parse_inputs("W_TOP_RSC.txt")
WBRSC = parse_inputs("W_BOT_RSC.txt")
WTFF = parse_inputs("W_TOP_FF.txt")
WBFF =parse_inputs("W_BOT_FF.txt")
WTCV = parse_inputs("W_TOP_CV.txt")
WBCV = parse_inputs("W_BOT_CV.txt")

###Rivalry Data
#up-states
RT_W_IM = parse_inputs("R_TOP_W_IM.txt")
RB_W_IM = parse_inputs("R_BOT_W_IM.txt")
RT_W_RSC = parse_inputs("R_TOP_W_RSC.txt")
RB_W_RSC = parse_inputs("R_BOT_W_RSC.txt")
RT_W_FF = parse_inputs("R_TOP_W_FF.txt")
RB_W_FF =parse_inputs("R_BOT_W_FF.txt")
RT_W_CV = parse_inputs("R_TOP_W_CV.txt")
RB_W_CV = parse_inputs("R_BOT_W_CV.txt")
#down-states
RT_L_IM = parse_inputs("R_TOP_L_IM.txt")
RB_L_IM = parse_inputs("R_BOT_L_IM.txt")
RT_L_RSC = parse_inputs("R_TOP_L_RSC.txt")
RB_L_RSC = parse_inputs("R_BOT_L_RSC.txt")
RT_L_FF = parse_inputs("R_TOP_L_FF.txt")
RB_L_FF =parse_inputs("R_BOT_L_FF.txt")
RT_L_CV = parse_inputs("R_TOP_L_CV.txt")
RB_L_CV = parse_inputs("R_BOT_L_CV.txt")

#
# los = parse_inputs("E_BOT_LOS.txt") #time series loser wta
# win = parse_inputs("E_TOP_WIN.txt") #time series winner wta
# nrm = parse_inputs("e_top_norm1.txt") #time series normalization
# nrm2 = parse_inputs("e_bot_norm1.txt") #time series other normalization
# riv = parse_inputs("e_riv1_input.txt") #time series raw rivalry
# winr = parse_inputs("wta_sample_raster_W.txt") #spikes winner
# losr = parse_inputs("wta_sample_raster_L.txt") #spikes loser
# nrmr1 = parse_inputs("norm_sample_raster_1.txt") #spikes norm1
# nrmr2 = parse_inputs("norm_sample_raster_2.txt") #spikes norm2
# rivr = parse_inputs("e_riv1_raster.txt") #spikes rivalry
# doms = parse_inputs("dominance_times.txt") #dominance times
#
# #calculated means and stds for all neurons in norm, wta
bm = parse_inputs("wta_bot_means.txt")
tm = parse_inputs("wta_top_means.txt")
bs = parse_inputs("wta_bot_stds.txt")
ts = parse_inputs("wta_top_stds.txt")
nm = parse_inputs("nmz_means.txt")
ns = parse_inputs("nmz_stds.txt")
#rivalry
etwm = parse_inputs("riv_top_win_means.txt")
ebwm = parse_inputs("riv_bot_win_means.txt")
etlm = parse_inputs("riv_top_los_means.txt")
eblm = parse_inputs("riv_bot_los_means.txt")
etws = parse_inputs("riv_top_win_stds.txt")
ebws = parse_inputs("riv_bot_win_stds.txt")
etls = parse_inputs("riv_top_los_stds.txt")
ebls = parse_inputs("riv_bot_los_stds.txt")
all_2624_R = parse_inputs("riv_2624_rast.txt")
dom_2624_R = parse_inputs("riv_2624_dom_rast.txt")
sup_2624_R = parse_inputs("riv_2624_sup_rast.txt")
all_2624_I = parse_inputs("riv_2624_input.txt")
dom_2624_I = parse_inputs("riv_2624_IU.txt")
sup_2624_I = parse_inputs("riv_2624_ID.txt")

# fig = plt.figure()
#
# ax1 = fig.add_subplot(5,2,1)
# ax1.set_title("Distribution of Input Means Across Pools")
# ax1.hist(bm, 50)
# ax1.axvline(np.mean(bm), linestyle = "dashed", color="k")
# ax1.axvline(0.1, linestyle = "dashed", color="r")
# ax1.set_ylabel("Loser")
# # ax1.set_xticklabels([])
#
# ax3 = fig.add_subplot(5,2,3, sharex = ax1)
# ax3.hist(nm, 50)
# ax3.axvline(np.mean(nm), linestyle = "dashed", color="k")
# ax3.axvline(0.1, linestyle = "dashed", color="r")
# ax3.set_ylabel("Normalization")
# # ax3.set_xticklabels([])
#
#
# ax5 = fig.add_subplot(5,2,5, sharex = ax1)
# ax5.hist(tm, 50)
# ax5.axvline(np.mean(tm), linestyle = "dashed", color="k")
# ax5.axvline(0.1, linestyle = "dashed", color="r")
# ax5.set_ylabel("Winner")
# # ax5.set_xticklabels([])
#
# ax7 = fig.add_subplot(5,2,7, sharex = ax1)
# ax7.hist(etwm, 50)
# ax7.axvline(np.mean(etwm), linestyle = "dashed", color="k")
# ax7.axvline(0.1, linestyle = "dashed", color="r")
# ax7.set_ylabel("Up_state")
# # ax7.set_xticklabels([])
#
# ax9 = fig.add_subplot(5,2,9, sharex = ax1)
# ax9.hist(etlm, 50)
# ax9.axvline(np.mean(etlm), linestyle = "dashed", color="k")
# ax9.axvline(0.1, linestyle = "dashed", color="r")
# ax9.set_ylabel("Down_state")
# ax9.set_xlabel("Mean Input to Neuron")
# # ax9.set_xticklabels([-0.15, -0.1, -0.05, 0, 0.05, 0.1])
#
#
# ax2 = fig.add_subplot(5,2,2)
# ax2.set_title("Distribution of Input Standard Deviations Across Pools")
# ax2.hist(bs, 50)
# ax2.axvline(np.mean(bs), linestyle = "dashed", color="k")
#
# ax4 = fig.add_subplot(5,2,4, sharex = ax2)
# ax4.hist(ns, 50)
# ax4.axvline(np.mean(ns), linestyle = "dashed", color="k")
#
# ax6 = fig.add_subplot(5,2,6, sharex = ax2)
# ax6.hist(ts, 50)
# ax6.axvline(np.mean(ts), linestyle = "dashed", color="k")
#
# ax8 = fig.add_subplot(5,2,8, sharex = ax2)
# ax8.hist(etws, 50)
# ax8.axvline(np.mean(etws), linestyle = "dashed", color="k")
#
# ax10 = fig.add_subplot(5,2,10, sharex = ax2)
# ax10.hist(etls, 50)
# ax10.axvline(np.mean(etls), linestyle = "dashed", color="k")
# ax10.set_xlabel("Standard Deviation of Input to Neuron")
#
# axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
# for i in axes:
#     plt.setp(i.get_xticklabels(), visible = False)
#

#

fig = plt.figure()

ax1 = fig.add_subplot(5,1,1)
ax1.set_title("Distribution of Input Means Across Pools")
ax1.hist(bm, 50)
ax1.axvline(np.mean(bm), linestyle = "dashed", color="k")
ax1.axvline(0.1, linestyle = "dashed", color="r")
ax1.set_ylabel("Loser")
# ax1.set_xticklabels([])

ax3 = fig.add_subplot(5,1,2, sharex = ax1)
ax3.hist(nm, 50)
ax3.axvline(np.mean(nm), linestyle = "dashed", color="k")
ax3.axvline(0.1, linestyle = "dashed", color="r")
ax3.set_ylabel("Normalization")
# ax3.set_xticklabels([])


ax5 = fig.add_subplot(5,1,3, sharex = ax1)
ax5.hist(tm, 50)
ax5.axvline(np.mean(tm), linestyle = "dashed", color="k")
ax5.axvline(0.1, linestyle = "dashed", color="r")
ax5.set_ylabel("Winner")
# ax5.set_xticklabels([])

ax7 = fig.add_subplot(5,1,4, sharex = ax1)
ax7.hist(etwm, 50)
ax7.axvline(np.mean(etwm), linestyle = "dashed", color="k")
ax7.axvline(0.1, linestyle = "dashed", color="r")
ax7.set_ylabel("Up_state")
# ax7.set_xticklabels([])

ax9 = fig.add_subplot(5,1,5, sharex = ax1)
ax9.hist(etlm, 50)
ax9.axvline(np.mean(etlm), linestyle = "dashed", color="k")
ax9.axvline(0.1, linestyle = "dashed", color="r")
ax9.set_ylabel("Down_state")
ax9.set_xlabel("Mean Input to Neuron")
# ax9.set_xticklabels([-0.15, -0.1, -0.05, 0, 0.05, 0.1])

axes = [ax1, ax3, ax5, ax7]
for i in axes:
    plt.setp(i.get_xticklabels(), visible = False)















#
# fig = plt.figure()
#
# ax1 = fig.add_subplot(3,2,1)
# ax1.set_title("Distribution of Input Means Across Pools")
# ax1.hist(bm, 50)
# ax1.axvline(np.mean(bm), linestyle = "dashed", color="k")
# ax1.axvline(0.1, linestyle = "dashed", color="r")
# ax1.set_ylabel("Loser")
# # ax1.set_xticklabels([])
#
#
# ax3 = fig.add_subplot(3,2,3, sharex = ax1)
# ax3.hist(nm, 50)
# ax3.axvline(np.mean(nm), linestyle = "dashed", color="k")
# ax3.axvline(0.1, linestyle = "dashed", color="r")
# ax3.set_ylabel("Normalization")
# # ax3.set_xticklabels([])
#
#
# ax5 = fig.add_subplot(3,2,5, sharex = ax1)
# ax5.hist(tm, 50)
# ax5.axvline(np.mean(tm), linestyle = "dashed", color="k")
# ax5.axvline(0.1, linestyle = "dashed", color="r")
# ax5.set_ylabel("Winner")
# ax5.set_xlabel("Mean Input to Neuron")
# # ax5.set_xticklabels([])
#
#
# ax2 = fig.add_subplot(3,2,2)
# ax2.set_title("Distribution of Input Standard Deviations Across Pools")
# ax2.hist(bs, 50)
# ax2.axvline(np.mean(bs), linestyle = "dashed", color="k")
#
# ax4 = fig.add_subplot(3,2,4, sharex = ax2)
# ax4.hist(ns, 50)
# ax4.axvline(np.mean(ns), linestyle = "dashed", color="k")
#
# ax6 = fig.add_subplot(3,2,6, sharex = ax2)
# ax6.hist(ts, 50)
# ax6.axvline(np.mean(ts), linestyle = "dashed", color="k")
# ax6.set_xlabel("Standard Deviation of Input to Neuron")










# fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True)
# ax1.set_title("Histogrammed Input to Normalization Neuron")
# ax1.hist(nrm, 100)
# ax1.axvline(np.mean(nrm), color = 'k', linestyle = 'dashed', linewidth = 2)
# ax1.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
#
# ax3.hist(win, 100)
# ax3.axvline(np.mean(win), color = 'k', linestyle = 'dashed', linewidth = 2)
# ax3.set_title("Histogrammed Input to Winning Neuron")
# ax3.set_xlabel("Current (mV/dt)")
# ax3.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
#
# ax2.hist(los, 100)
# ax2.axvline(np.mean(los), color = 'k', linestyle = 'dashed', linewidth = 2)
# ax2.set_title("Histogrammed Input to Losing Neuron")
# ax2.set_ylabel("Frequency")
# ax2.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
#
# plt.show()

#set order to be loser, normalization, winner

####################################################### INPUT & RASTER #########################
'''
fig = plt.figure()
ax1 = fig.add_subplot(5, 2, 1)
ax3 = fig.add_subplot(5, 2, 3, sharex = ax1)
ax5 = fig.add_subplot(5, 2, 5, sharex = ax1)
ax7 = fig.add_subplot(5, 2, 7, sharex = ax1)
ax9 = fig.add_subplot(5, 2, 9, sharex = ax1)
ax2 = fig.add_subplot(5, 2, 2)
ax4 = fig.add_subplot(5, 2, 4, sharex = ax2)
ax6 = fig.add_subplot(5, 2, 6, sharex = ax2)
ax8 = fig.add_subplot(5, 2, 8, sharex = ax2)
ax10 = fig.add_subplot(5, 2, 10, sharex = ax2)

ax1.set_title("Histogrammed Input")
ax1.hist(los, 100)
ax1.axvline(np.mean(los), color = 'k', linestyle = 'dashed', linewidth = 2)
ax1.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
ax1.set_ylabel("Loser")

ax5.hist(win, 100)
ax5.axvline(np.mean(win), color = 'k', linestyle = 'dashed', linewidth = 2)
# ax5.set_title("Histogrammed Input to Winning Neuron")
ax5.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
ax5.set_ylabel("Winner")

ax3.hist(nrm, 100)
ax3.axvline(np.mean(nrm), color = 'k', linestyle = 'dashed', linewidth = 2)
# ax3.set_title("Histogrammed Input to Losing Neuron")
ax3.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
ax3.set_ylabel("Normalization")
# all_2624_R = parse_inputs("riv_2624_rast.txt")
# dom_2624_R = parse_inputs("riv_2624_dom_rast.txt")
# sup_2624_R = parse_inputs("riv_2624_sup_rast.txt")
# all_2624_I = parse_inputs("riv_2624_input.txt")
# dom_2624_I = parse_inputs("riv_2624_IU.txt")
# sup_2624_I = parse_inputs("riv_2624_ID.txt")
ax7.hist(dom_2624_I, 100)
ax7.axvline(np.mean(dom_2624_I), color = 'k', linestyle = 'dashed', linewidth = 2)
ax7.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
ax7.set_ylabel("Up-State")

ax9.hist(sup_2624_I, 100)
ax9.axvline(np.mean(sup_2624_I), color = 'k', linestyle = 'dashed', linewidth = 2)
ax9.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
ax9.set_ylabel("Down-State")

# ax7.hist(riv, 100)
# ax7.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
# ax7.axvline(np.mean(riv), color = 'k', linestyle = 'dashed', linewidth = 2)
# ax7.set_xlabel("Current (mV/dt)")
# ax7.set_ylabel("Rivalry")

ax2.plot(losr*.1/1000., np.zeros(len(losr)), "|", ms = 5.)
ax2.set_title("Corresponding Rasters")
ax2.set_yticklabels([])
ax6.plot(winr*.1/1000., np.zeros(len(winr)), "|", ms = 5.)
ax4.set_yticklabels([])
ax4.plot(nrmr1*.1/1000., np.zeros(len(nrmr1)), "|", ms = 5.)
ax6.set_yticklabels([])
# ax8.plot(rivr*.1/1000., np.zeros(len(rivr)), "|", ms = 5.)
# ax8.set_xlabel("Time (s)")
ax8.plot(dom_2624_R*1./10000., np.zeros(len(dom_2624_R)), "|", ms = 5.)
ax8.set_yticklabels([])
ax10.plot((sup_2624_R*1./10000) -8, np.zeros(len(sup_2624_R)), "|", ms = 5.)

axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
for i in axes:
    plt.setp(i.get_xticklabels(), visible = False)



plt.show()
'''
# f, (ax1, ax2) = plt.subplots(1, 2)
# ax1.hist(riv, 100)
# ax2.plot(rivr, np.zeros(len(rivr)), "|", ms = 5.)
# ax1.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
# ax1.axvline(np.mean(riv), color = 'k', linestyle = 'dashed', linewidth = 2)
#
# #make a plot with histogram and raster
# #beneath it show up state hist and raster
# #beneath it show down state hist and raster
# #using e_bot so sequence goes 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0; 0= lose, 1 = win
# dflg = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
# def ligate_inputs(t, flg, inputs):
#     xw = np.zeros(1)
#     xl = np.zeros(1)
#     for i in range(len(t)-1):
#         print i
#         if flg[i] == 0:
#             xwi = inputs[t[i]:t[i+1]]
#             xw = np.concatenate([xw, xwi])
#         if flg[i] == 1:
#             xli = inputs[t[i]:t[i+1]]
#             xl = np.concatenate([xl, xli])
#     return xw, xl
#
# xw, xl = ligate_inputs(doms, dflg, riv)
#
# fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex = True)
# ax2.hist(riv, 100)
# ax2.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
# ax2.axvline(np.mean(riv), color = 'k', linestyle = 'dashed', linewidth = 2)
# ax2.set_yticklabels([])
# ax2.set_ylabel("All")
#
# ax3.hist(xw, 100)
# ax3.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
# ax3.axvline(np.mean(xw), color = 'k', linestyle = 'dashed', linewidth = 2)
# ax3.set_yticklabels([])
# ax3.set_ylabel("Up-State")
#
# ax1.set_title("Distribution of Inputs During Rivalry")
# ax1.hist(xl, 100)
# ax1.axvline(0.1, color = 'r', linestyle = 'dashed', linewidth = 2)
# ax1.axvline(np.mean(xl), color = 'k', linestyle = 'dashed', linewidth = 2)
# ax1.set_yticklabels([])
# ax1.set_ylabel("Down-State")
#
# plt.show()

# ws = np.sort(win)
# ls = np.sort(los)
# ns = np.sort(nrm)
#
#
#
# plt.plot(ns, ns, "g")
# plt.plot(ws, ns, "r")
# plt.plot(ls, ns, "b")
# plt.title("QQ Plot: Nomralization vs Winning and Losing Samples")
# plt.show()
#

#


#








#histogram of means in normalization, winning, and losing pools
