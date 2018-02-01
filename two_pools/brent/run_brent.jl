#run brent small


N = 5000
half = round(Int64, N/2)
quar = round(Int64, N/4)

runtime = 200000
rt = runtime /1000.
h = 0.1
rot = round(Int64, runtime/h)
rot2 = rot*2

# BW = Brent_W(N, 1/8.)
# CSR = sparse_rep(BW, N)

scan = linspace(.2, .8, 20)
MD = []
alt = []

stdev = .1*10
tau_n = 100*10
mu = 0.2

R1 = rand(Normal(0, stdev), rot2)
R2 = rand(Normal(0, stdev), rot2)

s1 = OU_Model(R1, tau_n, 0.1)
s2 = OU_Model(R2, tau_n, 0.1)

s1 .+= mu
s2 .+= mu

s1 = s1[rot:end]
s2 = s2[rot:end]

STD = mean([std(s1), std(s2)])
stdev/(sqrt(2)*sqrt(tau_n))

A = s1./(s1 .+ s2)

#stdev_resulting = stdev/(sqrt(2)*sqrt(tau_n))

@time t, r = Brent_Network_Euler_CSR(h, runtime, CSR, BW, N, s1, s2, 20., 20., 2.)

e_m = find(r .<= half)
i_m = find(r .> half)

te_pt = t[e_m]
re_pt = r[e_m]
ti_pt = t[i_m]
ri_pt = r[i_m]

ntd, nts = nt_diff(te_pt, re_pt, rot, quar, 25/h)
s = ntd ./ nts

# push!(sm, mean(s))
# push!(ss, std(s))

tl = .3
th = .7

flags, times = WLD_01(s, tl, th)
top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times)
m = (rot - tnmz)/(length(top) + length(bot))
MDT = tdom/length(top)
MDB = bdom/length(bot)
ad = []
for i in eachindex(top)
    push!(ad, top[i][2]-top[i][1])
end
for i in eachindex(bot)
    push!(ad, bot[i][2]-bot[i][1])
end
ad = convert(Array{Float64}, ad)/10000.
adx = find(ad .> 0.2)
adz = ad[adx]
cvd = cv(ad)
println("##RESULT $(MDT), $(MDB), $(m), $(cvd), $(length(times))")
push!(alt, length(times))
push!(MD, mean(ad))

plt[:hist](adz, 20)


#
#
# if flags == ["lose", "end"]
#     println("##RESULT $(0), $(rot), $(rot), $(0)")
# elseif flags == ["win", "end"]
#     println("##RESULT $(rot), $(0), $(rot), $(0)")
# elseif flags == ["draw", "end"]
#     println("##RESULT $(0), $(0), $(0), $(0)")
# else
#     t2, f2 = splice_reversions(flags, times)
#     top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times)
#
#     # d = convert(Array{Float64}, diff(250*t2))
#     # m = mean(d)
#     # cvd = cv(d)
#     # fw = find(f2.=="win")
#     # fl = find(f2.=="lose")
#     # MDT = mean(d[fw])
#     # MDB = mean(d[fl])
#     m = (rot - tnmz)/(length(top) + length(bot))
#     MDT = tdom/length(top)
#     MDB = bdom/length(bot)
#     ad = []
#     for i in eachindex(top)
#         push!(ad, top[i][2]-top[i][1])
#     end
#     for i in eachindex(bot)
#         push!(ad, bot[i][2]-bot[i][1])
#     end
#     ad = convert(Array{Float64}, ad)/10000.
#     cvd = cv(ad)
#
#     println("##RESULT $(MDT), $(MDB), $(m), $(cvd)")
# end

# subplot(311)
# ax1 = gca()
# title("Simulation of Brent's Network")
# plot(te_pt./10000, re_pt, "g.", ms = 1.)
# plot(t2*25/1000, fill(quar, length(t2)), "k.", ms = 10.)
# subplot(312)
# ax2 = gca()
# plot(linspace(0, rt, length(s1)), s1)
# plot(linspace(0, rt, length(s1)), s2)
# ylabel("Inputs")
# subplot(313)
# plot(linspace(0, rt, length(s)), s)
# xlabel("Time (s)")
# ylabel("Activity")
# axhline(0.5, color = "r", linestyle = "dashed")
# plot(t2*25/1000., fill(0.5, length(t2)), "k.", ms = 10.)
# plot(linspace(0, rt, length(s)), s)
# xlabel("Time (s)")
# ylabel("Activity")
# title("Normalized Activity Over Time")
# axhline(0.7, linestyle = "dashed", color="r")
# axhline(0.3, linestyle = "dashed", color="r")
# axhline(0.5, linestyle = "dashed", color="g")
# plot([i[1]/(1000/h) for i in top], fill(0.7, length(top)), "k.", ms = 10.)
# plot([i[1]/(1000/h) for i in bot], fill(0.3, length(bot)), "k.", ms = 10.)
# plot([i[1]/(1000/h) for i in nmz], fill(0.5, length(nmz)), "k.", ms = 10.)
#

# plot(t2[fw]*25/1000., fill(0.6, length(fw)), "k.", ms = 10.)
# plot(t2[fl]*25/1000., fill(0.4, length(fl)), "k.", ms = 10.)
