#run brent small

N = 5000
half = round(Int64, N/2)
quar = round(Int64, N/4)

runtime = 100000
rt = runtime /1000.
h = 0.1
rot = round(Int64, runtime/h)
rot2 = rot*2

stdev = 0.2
tau_n = 100

R1 = rand(Normal(0, stdev), rot2)
R2 = rand(Normal(0, stdev), rot2)

s1 = OU_Model(R1, tau_n, 0.1)
s2 = OU_Model(R2, tau_n, 0.1)

s1 .+= parse(Float64, ARGS[1])
s2 .+= 0.3

s1 = s1[rot:end]
s2 = s2[rot:end]

BW = Brent_W(5000, 1/8.)
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

tl = .4
th = .6

flags, times = WLD_01(s, tl, th)

if flags == ["lose", "end"]
    println("##RESULT $(0), $(rot), $(rot), $(0)")
elseif flags == ["win", "end"]
    println("##RESULT $(rot), $(0), $(rot), $(0)")
elseif flags == ["draw", "end"]
    println("##RESULT $(0), $(0), $(0), $(0)")
else
    t2, f2 = splice_reversions(flags, times)
    top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times)

    d = convert(Array{Float64}, diff(250*t2))
    m = mean(d)
    cvd = cv(d)
    fw = find(f2.=="win")
    fl = find(f2.=="lose")
    MDT = mean(d[fw])
    MDB = mean(d[fl])

    println("##RESULT $(MDT), $(MDB), $(m), $(cvd)")
end
