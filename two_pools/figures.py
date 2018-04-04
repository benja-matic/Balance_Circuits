import numpy as np
import matplotlib.pyplot as plt

def parse_inputs(f):
    with open(f,'r') as x:
        a = x.read()
    b = a.split('\n')
    c = np.array([float(b[i]) for i in range(len(b)-1)])
    return c

def read_raster(x):
    with open(x, 'r') as f:
        a = f.read()
    b = a.split('\n')
    b.pop()
    c = [i.split(',') for i in b]
    te = np.array([float(i[0]) for i in c])
    re = np.array([float(i[1]) for i in c])
    return te, re

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




fig, ax = plt.subplots(3, sharex = True)
ax[0].set_title("Distribution of Mean Inputs")
ax[0].hist(NTIM, 50)
ax[0].set_ylabel("Normalization")
ax[0].axvline(np.mean(NTIM), linestyle="dashed", color = "g")
ax[0].axvline(0.1, linestyle="dashed", color = "r")
ax[1].hist(WBIM, 50)
ax[1].set_ylabel("Winner")
ax[1].axvline(np.mean(WBIM), linestyle="dashed", color = "g")
ax[1].axvline(0.1, linestyle="dashed", color = "r")
ax[2].hist(WTIM, 50)
ax[2].set_ylabel("Loser")
ax[2].axvline(np.mean(WTIM), linestyle="dashed", color = "g")
ax[2].axvline(0.1, linestyle="dashed", color = "r")

fig, ax = plt.subplots(1,2, sharex=True)
ax[0].hist(NTFF, 50)
ax[1].hist(NTCV, 50)
ax[0].set_title("Fano Factor")
ax[1].set_title('$CV_{isi}$')
ax[1].axvline(np.mean(NTCV), linestyle="dashed", color = "g")
ax[1].axvline(1., linestyle="dashed", color = "r")
ax[0].axvline(np.mean(NTFF), linestyle="dashed", color = "g")
ax[0].axvline(1., linestyle="dashed", color = "r")


fig, ax = plt.subplots(1)
ax.hist(NTRSC, 50)
ax.set_title('$r_{sc}$')
ax.axvline(0., linestyle = "dashed", color = "r", alpha = .5, lw = 3.)
ax.axvline(np.mean(NTRSC), linestyle = "dashed", color = "g")

fig, ax = plt.subplots(1,3)
ax[0].hist(NTFF, 50)
ax[1].hist(NTCV, 50)
ax[0].set_title("Fano Factor")
ax[1].set_title('$CV_{isi}$')
ax[1].axvline(np.mean(NTCV), linestyle="dashed", color = "g")
ax[1].axvline(1., linestyle="dashed", color = "r")
ax[0].axvline(np.mean(NTFF), linestyle="dashed", color = "g")
ax[0].axvline(1., linestyle="dashed", color = "r")
ax[2].hist(NTRSC, 50)
ax[2].set_title('$r_{sc}$')
ax[2].axvline(0., linestyle = "dashed", color = "r", alpha = .5, lw = 3.)
ax[2].axvline(np.mean(NTRSC), linestyle = "dashed", color = "g")



fig, ax = plt.subplots()
ax.axis('off')
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_xticks([])
ax.set_yticks([])

ax0 = fig.add_subplot(2,2,1)
ax0.hist(WTCV, 50)
ax0.set_title('$CV_{isi}$')
ax0.axvline(np.mean(WTCV), linestyle="dashed", color = "g")
ax0.axvline(1., linestyle="dashed", color = "r")
ax0.set_ylabel("Loser")

ax1 = fig.add_subplot(2,2,2)
ax1.hist(WTRSC, 50)
ax1.set_title('$r_{sc}$')
ax1.axvline(np.mean(WTRSC), linestyle="dashed", color = "g")
ax1.axvline(0., linestyle="dashed", color = "r")

ax2 = fig.add_subplot(2,2,3, sharex = ax0)
ax2.hist(WBCV, 50)
ax2.axvline(np.mean(WBCV), linestyle="dashed", color = "g")
ax2.axvline(1., linestyle="dashed", color = "r")
ax2.set_ylabel("Winner")

ax3 = fig.add_subplot(2,2,4, sharex = ax1)
ax3.hist(WBRSC, 50)
ax3.axvline(np.mean(WBRSC), linestyle="dashed", color = "g")
ax3.axvline(0., linestyle="dashed", color = "r")



fig, ax = plt.subplots()
ax.axis('off')
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_xticks([])
ax.set_yticks([])

ax0 = fig.add_subplot(2,2,1)
ax0.hist(RT_W_CV, 50)
ax0.set_title('$cv_{isi}$')
ax0.axvline(1., linestyle="dashed", color = "r", lw = 3, alpha = .5)
ax0.axvline(np.mean(RT_W_CV), linestyle="dashed", color = "g")
ax0.set_ylabel("Up-State")

ax1 = fig.add_subplot(2,2,2)
ax1.hist(RT_W_RSC, 50)
ax1.set_title('$r_{sc}$')
ax1.axvline(np.mean(WTRSC), linestyle="dashed", color = "g")
ax1.axvline(0., linestyle="dashed", color = "r")

ax2 = fig.add_subplot(2,2,3, sharex = ax0)
ax2.hist(RT_L_CV, 50)
ax2.axvline(np.mean(RT_L_CV), linestyle="dashed", color = "g")
ax2.axvline(1., linestyle="dashed", color = "r")
ax2.set_ylabel("Down-State")

ax3 = fig.add_subplot(2,2,4, sharex = ax1)
ax3.hist(RT_L_RSC, 50)
ax3.axvline(np.mean(RT_L_RSC), linestyle="dashed", color = "g")
ax3.axvline(0., linestyle="dashed", color = "r")







longs = ['AR_long_rast.txt', 'Balance_long_rast.txt', 'Brent_long_rast.txt']

exec(open('C:\\Users\\cohenbp\\Documents\\Neuroscience\\Realistic_Simulations\\Hodgkin_Huxley_Adaptation\\read_rasters.py'))
ntotal = 5000000
tr, rr = read_raster('AR_long_rast.txt')
tb, rb = read_raster('Balance_long_rast.txt')
ts, rs = read_raster('Brent_long_rast.txt')
m = np.where(rs <= 2500)
ts = ts[m]
rs = rs[m]

k = np.where(tr < 100*10000)[0]
k0 = k[-1]
flags, times = WLD_01(sr, -1./3, 1./3)
tx = times * 1500./10000.
kt = np.where(tx < 100)[0]
kt0 = kt[-1]
plt.plot(tr[:k0]/10000., rr[:k0], "g.", ms = .1)
plt.plot(tx[:kt0], np.ones(len(tx[:kt0]))*1600., "r.")

k = np.where(tb < 100*10000)[0]
k0 = k[-1]
flags, times = WLD_01(sb, -1./3, 1./3)
tx = times * 1500./10000.
kt = np.where(tx < 100)[0]
kt0 = kt[-1]
plt.plot(tb[:k0]/10000., rb[:k0], "g.", ms = 1.)
plt.plot(tx[:kt0], np.ones(len(tx[:kt0]))*1600., "r.")

k = np.where(ts < 100*10000)[0]
k0 = k[-1]
flags, times = WLD_01(SS, -1./3, 1./3)
tx = times * 1500./10000.
kt = np.where(tx < 100)[0]
kt0 = kt[-1]
plt.plot(ts[:k0]/10000., rs[:k0], "g.", ms = 1.)
plt.plot(tx[:kt0], np.ones(len(tx[:kt0]))*1250., "r.")


half = 1600
netd_binsize = 1500

def d_get(sig):
    flags, times = WLD_01(sig, -1./3, 1./3)
    tx = times * 1500.
    dx = np.diff(tx)
    CV0 = CV(dx)
    lp = 3000. #300ms
    dlp = [i for i in dx if i > lp]
    return np.array(dlp)

sr = nt_diff_H(tr, rr, 5000000, 1600, 1500)
sb = nt_diff_H(tb, rb, 5000000, 1600, 1500)
SS = nt_diff_H(ts, rs, 5000000, 1250, 1500)


flags, times = WLD_01(sb, -1./3, 1./3)
dlp_B = d_get(sb)


dlp_R = d_get(sr)/10000.
dlp_B = d_get(sb)/10000.
dlp_S = d_get(SS)/10000.
fig, ax = plt.subplots(1,3, sharex = True)
ax[0].hist(dlp_S, 40)
ax[1].hist(dlp_R, 40)
ax[2].hist(dlp_B, 40)
plt.show()



SS = np.load("Brent_Alternations_Long.npy")
sr = np.load("Reg_alternations_long.npy")
sb = np.load("Bal_alternations_long.npy")












t, r = read_raster("Brent_Basic_Rast.txt")
s1 = parse_inputs("Brent_S1.txt")
s2 = parse_inputs("Brent_S2.txt")

m = np.where(r < 2500)[0]
te = t[m]
re = r[m]

fig, ax = plt.subplots(2,1, sharex = True)
ax[0].plot(te, re, "g.", ms = 1.)
ax[1].plot(s1, label = "Input 1")
ax[1].plot(s2, label = "Input 2")





DA = nt_diff_H(te, re, 200000, 1250, 1500)
S_diff = s1 - s2
s_sum = s1 + s2
S = S_diff/s_sum

























###BW rasters

contents = os.listdir(os.getcwd())
files = []
for i in contents:
    if i[-3:] == 'txt':
        files.append(i)

sim2 = [i for i in files if i[3] == '2']
sim3 = [i for i in files if i[3] == '3']
sim5 = [i for i in files if i[3] == '5']
data = np.zeros((12, len(sim5)))
data3 = np.zeros((12, len(sim3)))
data2 = np.zeros((12, len(sim2)))


for i in range(len(sim5)):
    t, r = read_raster(sim5[i])
    data[:, i] = get_stats(t, r, ntotal, half, netd_binsize, fbinsize, cbinsize)


for i in range(len(sim3)):
    t, r = read_raster(sim3[i])
    data3[:, i] = get_stats(t, r, ntotal, half, netd_binsize, fbinsize, cbinsize)


for i in range(len(sim2)):
    t, r = read_raster(sim2[i])
    data2[:, i] = get_stats(t, r, ntotal, half, netd_binsize, fbinsize, cbinsize)


s = np.array([float(i.split('_')[1]) for i in sim5])
s2 = np.array([float(i.split('_')[1]) for i in sim2])
s3 = np.array([float(i.split('_')[1]) for i in sim3])


fig, ax = plt.subplots(2, sharex = True)
ax[0].plot(s3, data3[5,:], "g.")
ax[1].plot(s3, data3[2,:], ".", label = "Dominances")
# ax[1].plot(s3, data3[8,:], ".", label = "ISI")
fit1 = np.polyfit(s3, data3[5,:], deg = 1)
fit2 = np.polyfit(s3, data3[2,:], deg = 1)
ax[0].plot(s3, (fit1[0]*s3) + fit1[1], 'r')
ax[1].plot(s3, (fit2[0]*s3) + fit2[1], 'r')
ax[1].set_xlabel("Input Strength")
ax[0].set_ylabel("Alternation Rate")
ax[1].set_ylabel("CV")

fig, ax = plt.subplots(1,2)
ax[0].plot(s3, data3[2,:], ".", label = "Dominances")
ax[0].plot(s3, data3[8,:], ".", label = "ISI")
fit1 = np.polyfit(s3, data3[2,:], deg = 1)
fit2 = np.polyfit(s3, data3[8,:], deg = 1)
ax[0].plot(s3, (fit1[0]*s3) + fit1[1], 'r')
ax[0].plot(s3, (fit2[0]*s3) + fit2[1], 'r')
ax[0].set_xlabel("Input Strength")
ax[0].set_ylabel("CV")
ax[0].legend()
ax[1].plot(data3[2,:], data3[8,:], ".")
ax[1].set_xlabel("CVD")
ax[1].set_ylabel("CVISI")
fit3 = np.polyfit(data3[2,:], data3[8,:], deg = 1)
ax[1].plot(data3[2,:], (fit3[0]*data3[2,:]) + fit3[1], 'r')
plt.tight_layout()

fig, ax = plt.subplots(2, sharex = True)
ax[0].plot(s, data[5,:], "g.")
ax[1].plot(s, data[2,:], ".", label = "Dominances")
# ax[1].plot(s, data3[8,:], ".", label = "ISI")
fit1 = np.polyfit(s, data[5,:], deg = 1)
fit2 = np.polyfit(s, data[2,:], deg = 1)
ax[0].plot(s, (fit1[0]*s) + fit1[1], 'r')
ax[1].plot(s, (fit2[0]*s) + fit2[1], 'r')
ax[1].set_xlabel("Input Strength")
ax[0].set_ylabel("Alternation Rate")
ax[1].set_ylabel("CV")
# ax[1].axhline(data[2,0], linestyle = "dashed")






te, re = read_raster(sim5[19])
sig = nt_diff_H(te, re, ntotal, half, netd_binsize)
flags, times = WLD_01(sig, -1./3, 1./3)
tx = times * netd_binsize
dx = np.diff(tx)
dlp = [i for i in dx if i > 30000]
t2, f2 = splice_reversions(flags, times)
d = np.diff(np.array(netd_binsize)*t2)/100000.














###########################L2


ntotal = 10000000.
import os
contents = os.listdir(os.getcwd())
files = []
for i in contents:
    if i[-3:] == 'txt':
        files.append(i)


sim5 = [i for i in files if i[3+3] == '5']
s = np.array([float(i.split('_')[2]) for i in sim5])


data = np.zeros((6, len(sim5)))

for i in range(len(sim5)):
    print i
    te, ri = read_raster(sim5[i])
    data[:,i] = get_L2(t, r, ntotal, half, netd_binsize, fbinsize, cbinsize)
    # res = get_L2(t, r, ntotal, half, netd_binsize, fbinsize, cbinsize)
    # for j in range(len(res)):
    #     data[j,i] = res[j]




# [MDT, MDB, len(times), MD2, MDT2, MDB2]
plt.plot(s, data[0,:], "r.", label = "Pool 1")
plt.plot(s, data[1,:], "b.", label = "Pool 2")
plt.plot(s, data[2,:], "g.", label = "Alternations")
