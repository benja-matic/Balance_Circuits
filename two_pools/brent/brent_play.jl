#Brent_Play
include("brent.jl")

runtime = 10000
rt = runtime /1000.
h = 0.1
rot = round(Int64, runtime/h)
rot2 = rot*2

stdev = .1
tau_n = 100
mu = .2

R1 = rand(Normal(0, stdev), rot2)
R2 = rand(Normal(0, stdev), rot2)

s1 = OU_Model(R1, tau_n, 0.1)
s2 = OU_Model(R2, tau_n, 0.1)

s1 .+= mu
s2 .+= mu

s1 = s1[rot:end]
s2 = s2[rot:end]

# STD = mean([std(s1), std(s2)])
# stdev/(sqrt(2)*sqrt(tau_n))
s3 = s1 .+ s2
A1 = s1./(s3)
A2 = s2./(s3)

BW = Brent_W(N, 1/8.)
CSR = sparse_rep(BW, N)

@time t, r = Brent_Network_Euler_CSR(h, runtime, CSR, BW, N, s1, s2, 20., 20., 2.)

e_m = find(r .<= half)
i_m = find(r .> half)

te_pt = t[e_m]
re_pt = r[e_m]
ti_pt = t[i_m]
ri_pt = r[i_m]

ntd, nts = nt_diff(te_pt, re_pt, rot, quar, 25/h)
s = ntd ./ nts

ma = mean(A2)
A2A = A2 .- ma

r = std(s)/std(A2)
A2A *= r
A2A += ma

xtime = linspace(0, rot/10000., length(A2A))

# subplot(311)
# title("Random Network Tracks Competitive Inputs")
# plot(te_pt./10000, re_pt, "g.", ms = 1.)
# xticks([])
# yticks([625, 1875.], ["Pool 1", "Pool 2"])
# subplot(312)
# ax1 = gca()
# plot(s)
# axhline(0.5, linestyle = "dashed", color = "k")
# xticks([])
# yticks([])
# yticks([])
# ylabel("Network Activity")
# subplot(313)
# plot(xtime, A2)
# axhline(0.5, linestyle = "dashed", color = "k")
# yticks([])
# ylabel("Input Activity")
#
#
#
# println("variance of A: $(var(A2))")
# println("variance of s: $(var(s))")
# println("var(s)/var(a): $(var(s)/var(A2))")
# println("sqrt(runtime): $(sqrt(length(s)))")



function zscore(x)
    m = mean(x)
    s = std(x)
    x1 = x .- m
    x2 = x1 ./s
    return x2
end

A2z = zscore(A2)
s2z = zscore(s)
subplot(211)
plot(s2z)
subplot(212)
plot(A2z)
