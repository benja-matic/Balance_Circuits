#mcb_r3_scratch
function von_mises_dist(x, k, mu, N)
  a = exp(k*cos(x-mu))/(N*besseli(0, k))
end

function weights(Ne,Ni,kee,kei,kie_L,kii,Aee,Aei,Aie_L,Aii,pee,pei,pie,pii)
# Construct weight functions

wee = zeros(Ne,Ne)
wei = zeros(Ne,Ni)
wie = zeros(Ni,Ne)
wii = zeros(Ni,Ni)

we = 2*collect(0:Ne-1)*pi/Ne
wi = 2*collect(0:Ni-1)*pi/Ni

for i = 1:Ne
    wee[i,:]=Aee*circshift(von_mises_dist(we, kee, 0, Ne), i-1)
    wei[i,:]=Aei*circshift(von_mises_dist(wi, kei, 0, Ni),div(Ni*i,Ne)-1) #
end
for i = 1:Ni
    wie[i,:]=Aie_L*circshift(von_mises_dist(we, kie_L, 0, Ne),div(Ne*i,Ni)-1)# + Aie_D*circshift(von_mises_dist(we, kie_D, pi, Ne),div(Ne*i,Ni)-1)
    wii[i,:]=Aii*circshift(von_mises_dist(wi, kii, 0, Ni),i-1)
end

# keep connections with probability pab
wee = wee.*(rand(Ne,Ne).<pee)
wei = wei.*(rand(Ne,Ni).<pei)
wie = wie.*(rand(Ni,Ne).<pie)
wii = wii.*(rand(Ni,Ni).<pii)
#println("weights complete")

return wee, -wei, wie, -wii
end

function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  ###x = dv/dt
  ###vth = v2 - xt
  ###(vth - v2)/x = t
  return t
end

function euler_lif(h, total, Ne, Ni, wee, wei, wie, wii, s1, s2, vth, tau_m, tau_ee, tau_ei, tau_ie, tau_ii, tau_ae, tau_ai, g_e, g_i)

  ntotal = round(Int,total/h)
  srand(1234)

  #initial conditions
  ve = rand(Ne)*vth
  vi = rand(Ni)*vth
  ve_buff = ve
  vi_buff = vi
  #input = zeros(Ne, ntotal)
  e_top = zeros(ntotal)
  e_bot = zeros(ntotal)
  # ev = zeros(ntotal)
  adapt_top = zeros(ntotal)
  adapt_bot = zeros(ntotal)

  #storage space
  see = zeros(Ne)#synapse
  sei = zeros(Ne)
  sie = zeros(Ni)
  sii = zeros(Ni)
  ae = zeros(Ne)#adaptation
  ai = zeros(Ni)
  #raster
  time_e = Float64[0]
  raster_e = Float64[0]
  time_i = Float64[0]
  raster_i =Float64[0]
  #Calculate drive and leak constants ahead of time
  m_leak = h/tau_m
  ee_leak = h/tau_ee
  ei_leak = h/tau_ei
  ie_leak = h/tau_ie
  ii_leak = h/tau_ii
  ae_leak = h/tau_ae
  ai_leak = h/tau_ai
  #Set up drive to two pools special for competitive network
  drive = zeros(Ne)
  FPe = div(Ne,5) #each pool is 2/5 of the network
  P1s = round(Int,1+Ne/4-FPe)
  P1e = round(Int,1+Ne/4+FPe)
  P2s = round(Int,3*Ne/4-FPe)
  P2e = round(Int,3*Ne/4+FPe)
  drive[P1s:P1e] = s1*h
  drive[P2s:P2e] = s2*h
  half = div(Ne, 2)
  quarter = div(half, 2)
  rt = half + quarter
  rb = quarter

  kill_flag = false
  for iter = 1:ntotal
    #administer drive, leak, synapse, and adaptation to voltage
    ve += drive + (h*see) + (h*sei) - (ve*m_leak) - (ae*g_e)
    vi += (h*sie) + (h*sii) -(vi*m_leak) - (ai*g_i)
    #ev[iter] = ve[rt]
    #input[:,iter] = drive + (h*see[:]) + (h*sei[:]) - (ae[:]*g_e)
    #e_top[iter] = drive[rt] + (h*see[rt]) + (h*sei[rt]) - (ve[rt]*m_leak) - (ae[rt]*g_e)
    #e_bot[iter] = drive[rb] + (h*see[rb]) + (h*sei[rb]) - (ve[rb]*m_leak) - (ae[rb]*g_e)
    e_top[iter] = drive[rt] + (h*see[rt]) + (h*sei[rt]) - (ae[rt]*g_e)
    e_bot[iter] = drive[rb] + (h*see[rb]) + (h*sei[rb]) - (ae[rb]*g_e)
    adapt_top[iter] = sum(ae[P2s:P2e])
    adapt_bot[iter] = sum(ae[P1s:P1e])
    #administer leak to synapse
    see -= (see*ee_leak)
    sei -= (sei*ei_leak)
    sie -= (sie*ie_leak)
    sii -= (sii*ii_leak)
    #decay of adaptation
    ae -= (ae*ae_leak)
    ai -= (ai*ai_leak)
    #check for spikes
    ves = (ve.>vth)
    vis = (vi.>vth)
    vesum = sum(ves)
    visum = sum(vis)
    #update excitatory synapses and adaptation
    if vesum > 0
      spe = find(ves)
      for j = 1:vesum
        #interpolate spike time
        delta_h = interpolate_spike(ve[spe[j]], ve_buff[spe[j]], vth) #time since the spike (units of h)
        lee = exp(delta_h/tau_ee)
        lie = exp(delta_h/tau_ie)
        see += wee[:, spe[j]]*lee #modify amplitude of synapse by decay since estimated time of spike
        sie += wie[:, spe[j]]*lie
        ae[spe[j]] +=1
        push!(raster_e,spe[j])
  	    push!(time_e,iter-delta_h)
      end
    end
    #update inhibitory synapses and adaptation
    if visum > 0
      spi = find(vis)
      for j = 1:visum
        #interpolate spike time
        delta_h = interpolate_spike(vi[spi[j]], vi_buff[spi[j]], vth)
        lei = exp(delta_h/tau_ei)
        lii = exp(delta_h/tau_ii)
        sei += wei[:,spi[j]]*lei #modify amplitude of synapse by decay since estimated time of spike
        sii += wii[:,spi[j]]*lii
        ai[spi[j]] +=1
        push!(raster_i,spi[j])
  	    push!(time_i,iter-delta_h)
      end
    end
    #reset after spike
    ve -= vth*ves
    vi -= vth*vis
    #
    ve_buff = ve
    vi_buff = vi
    if iter == 1000
      if (length(raster_e) > 200*Ne*h*iter*.5*(1/1000)) | (length(raster_i) > 200*Ne*h*iter*.5*(1/1000))
        kill_flag = true
        break
      end
    end
  end



  return time_e, raster_e, time_i, raster_i, kill_flag, e_top, e_bot, adapt_top, adapt_bot
end

function cv(isi)
  SD = std(isi)
  return SD/mean(isi)
end

function FANO(aRAY)
  return var(aRAY)/mean(aRAY)
end

function tscore(mu, sigma)
  tsc = (20. - mu)/sigma
  #20 is voltage threshold
  return tsc
end

function emptiness(x, funk, error_code) #checks if an empty set got populated
  if length(x) > 0
    return funk(x)
  else
    return error_code
  end
end

function zeroness(x, funk, error_code) #checks if an initialized array of 0s got populated
  if Set(x) == Set([0])
    return error_code
  else
    return funk(x)
  end
end

#NTF = Network Fano Factor
#Discretize time into bins, loop over bins and get counts of spikes in a set of neurons in each time bin
#Calculate the variability (i.e. Fano Factor = variance/mean) of this count series
function ntf(t, r, ntotal, half, netd_binsize)

  t = t[2:end]
  r = r[2:end]
  netd_bins = collect(1:netd_binsize:ntotal)
  ntd = zeros(length(netd_bins)-1)
  nts = zeros(length(netd_bins)-1)

  for j = 2:length(netd_bins)
    tf = find(netd_bins[j-1] .<= t .< netd_bins[j])
    T = sum(r[tf] .> half)
    B = sum(r[tf] .<= half)
    NTD = T - B #################Differences in Spikes
    ntd[j-1] = NTD
    nts[j-1] = T+B
  end

  return ntd, nts
end

#Discretize time into bins, loop over bins and return (A)/(A+B)
#This returns a number between 0 and 1
#Alternatively, Heeger's technique is (A-B)/(A+B) which returns a number between -1 and 1
#Later use a shmitt trigger olgorithm on this quantity to infer dominance times offline
function nt_diff(t, r, ntotal, half, netd_binsize)

  t = t[2:end]
  r = r[2:end]
  netd_bins = collect(1:netd_binsize:ntotal)
  ntd = zeros(length(netd_bins)-1)
  nts = zeros(length(netd_bins)-1)

  for j = 2:length(netd_bins)
    tf = find(netd_bins[j-1] .<= t .< netd_bins[j])
    T = sum(r[tf] .> half)
    #B = sum(r[tf] .<= half)
    #NTD = T - B #################Differences in Spikes
    ntd[j-1] = T
    nts[j-1] = length(tf)
    #println(T+B==length(tf))
  end

  return ntd, nts
end

function nt_diff_H(t, r, ntotal, half, netd_binsize)

  netd_bins = collect(1:netd_binsize:ntotal)
  ntd = zeros(length(netd_bins)-1)
  nts = zeros(length(netd_bins)-1)

  for j = 2:length(netd_bins)
    tf = find(netd_bins[j-1] .<= t .< netd_bins[j])
    T = sum(r[tf] .> half)
    B = sum(r[tf] .<= half)
    NTD = T - B #################Differences in Spikes
    ntd[j-1] = NTD
    nts[j-1] = length(tf)
    #println(T+B==length(tf))
  end

  return ntd, nts
end

#Given a an array with "win" and "lose" lables, and the corresponding time points...
#Combine the losses and wins so you can easily analyze up-state or down-state activity for a single pool
function splice_exons(intdices, tnb)
  shmitt_exons = [[Int64(intdices[i]), Int64(intdices[i+1])] for i = 1:length(intdices)-1]
  top_time = []
  bot_time = []
  tnbs = tnb[1:end-1]
  dt = diff(intdices)
  tdom = emptiness(dt[tnbs], sum, 0)
  bdom = emptiness(dt[!tnbs], sum, 0)

  for i in eachindex(shmitt_exons)
    if tnb[i] == true
      push!(top_time, [intdices[i], intdices[i+1]])
    else
      push!(bot_time, [intdices[i], intdices[i+1]])
    end
  end
  return top_time, tdom, bot_time, bdom
end

#Part of the shmitt trigger
#given a number and two thresholds, determine which region the number is in
function comp_01(x, th, tl)
  if x > th
    return "win"
  elseif tl <= x <= th
    return "draw"
  elseif x < tl
    return "lose"
  else
    return "weird"
  end
end

#s is A/(A+B) from nt_diff above
#loop over each time bin and get "win", "lose", or "draw" labels from the comp_01 function above
#Thresholds of .3 and .7 are hard-coded here
function WLD_01(s)
  if maximum(s) < .3
    return ["lose", "end"], [1, length(s)]
  elseif minimum(s) >= .7
    return ["win", "end"], [1, length(s)]
  elseif (.3 <= maximum(s) <= .7) & (.3 <= minimum(s) <= .7)
    return ["draw", "end"], [1, length(s)]
  end
  times = []
  flags = []
  if s[1] >=.7
    flag = "win"
  elseif .3 <= s[1] <= .7
    flag = "draw"
  elseif s[1] < .3
    flag = "lose"
  end

  push!(times, 1)
  push!(flags, flag)

  s2 = s[2:end]
  for i in eachindex(s2)
    f1 = comp_01(s2[i], .7, .3)
    if f1 != flag
      flag = f1
      push!(times, i+1)
      push!(flags, flag)
    end
  end
  push!(times, length(s))
  push!(flags, "end")
  return flags, times
end

function WLD_01(s, tl, th)
  if maximum(s) < tl
    return ["lose", "end"], [1, length(s)]
  elseif minimum(s) >= th
    return ["win", "end"], [1, length(s)]
  elseif (tl <= maximum(s) <= th) & (tl <= minimum(s) <= th)
    return ["draw", "end"], [1, length(s)]
  end
  times = []
  flags = []
  if s[1] >= th
    flag = "win"
  elseif tl <= s[1] <= th
    flag = "draw"
  elseif s[1] < tl
    flag = "lose"
  end

  push!(times, 1)
  push!(flags, flag)

  s2 = s[2:end]
  for i in eachindex(s2)
    f1 = comp_01(s2[i], tl, th)
    if f1 != flag
      flag = f1
      push!(times, i+1)
      push!(flags, flag)
    end
  end
  push!(times, length(s))
  push!(flags, "end")
  return flags, times
end

function comp_01(x, tl, th)
  if x > th
    return "win"
  elseif tl <= x <= th
    return "draw"
  elseif x < tl
    return "lose"
  else
    return "weird"
  end
end

#WLD_01 but wihtout hard coding anything
function WLD_A(s, th, tl)
  if maximum(s) < tl
    return ["lose", "end"], [1, length(s)]
  elseif minimum(s) >= th
    return ["win", "end"], [1, length(s)]
  elseif (tl <= maximum(s) <= th) & (tl <= minimum(s) <= th)
    return ["draw", "end"], [1, length(s)]
  end
  times = []
  flags = []
  if s[1] >=th
    flag = "win"
  elseif tl <= s[1] <= th
    flag = "draw"
  elseif s[1] < tl
    flag = "lose"
  end

  push!(times, 1)
  push!(flags, flag)

  s2 = s[2:end]
  for i in eachindex(s2)
    f1 = comp_01(s2[i], th, tl)
    if f1 != flag
      flag = f1
      push!(times, i+1)
      push!(flags, flag)
    end
  end
  push!(times, length(s))
  push!(flags, "end")
  return flags, times
end

#get the total time network spends in each state (pool 1 dominates, pool 2 dominates, normalization)
function dom_time(flags, times)
  mixed = []
  win = []
  lose = []
  for i=1:length(times)-1
    if flags[i] == "draw"
      if (flags[i+1] == "win") | (flags[i+1] == "lose")
        push!(mixed, times[i+1] - times[i])
      end
    elseif flags[i] == "win"
      push!(win, times[i+1] - times[i])
    elseif flags[i] == "lose"
      push!(lose, times[i+1] - times[i])
    end
  end
  mixed = convert(Array{Float64}, mixed)
  win = convert(Array{Float64}, win)
  lose = convert(Array{Float64}, lose)
  return mixed, win, lose
end

#If the system is clearly rivaling and any reversions (going from dominance to normalization and back again) are short...
#Then loop over our labeled time series, get rid of reversions, and count it as dominance until the other pool dominates
#You can do the same thing by changing the thresholds, but you'd lose information about normalization (i.e. when is neither pool dominant)
#This way, we preserve information about normalization, and then we can throw it out if we decide it isn't important
#This is preferable to never collecting normalization information in the first place.
function splice_reversions(flags, times)
  w0 = findfirst(flags, "win")
  l0 = findfirst(flags, "lose")
  #assume rivalry, and that any reversions are short and fail to persist
  empezar = min(w0, l0)
  nf = [flags[empezar]]
  a = ["win", "lose"]
  t = [empezar]
  if start == w0
    f = 2
  else
    f = 1
  end
  flag = a[f]
  for i = empezar:length(flags)
    if flags[i] == a[f]
      push!(t, times[i])
      push!(nf, flags[i])
      a = circshift(a, 1)
    end
  end
  push!(t, times[end])
  return t, nf
end

#Not sure if I finished this function...can't remember
#I think this loops over a labeled time series and deletes dominances that are too short
#Problem with this is that then each time you delete, you have to start over and re-analyze the whole time series
#Otherwise the labels get mangled
function diff_eater(t,f, mint)
  tc = copy(t)
  fc = copy(f)
  f1 = f[1]
  c1 = t[1]
  wflag = false
  for i=1:length(t)-1
    td = (tc[i+1] - tc[i])
    if td < mint

      wflag = true
      deleteat!(tc, i+1)
      deleteat!(fc, i+1)
      if i != 1
        deleteat!(tc, i)
        deleteat!(fc, i)
      end
      break
    end
  end
  return tc, fc, wflag
end

#So I think this was designed to run a while-loop and call diff_eater, and then start over every time you eat
#At some point I think I decided this was more trouble than it was worth
function low_pass_DD(flags, times, mint)
  wflag = true
  tn, fn, wflag = diff_eater(times, flags, mint)
  while wflag == true
    tn, fn, wflag = diff_eater(tn, fn, mint)
  end
  return tn, fn
end

#Another way of combining wins, losses, and draws and calculating time spent in each regime
function splice_flags(flags, times)
  l = find(flags.=="lose")
  w = find(flags.=="win")
  d = find(flags.=="draw")
  bot = [[((times[i]-1)*250)+1, ((times[i+1]-1)*250)] for i in l]
  top = [[((times[i]-1)*250)+1, ((times[i+1]-1)*250)] for i in w]
  nmz = [[((times[i]-1)*250)+1, ((times[i+1]-1)*250)] for i in d]
  # bot = [[times[i], times[i+1]] for i in l]*250
  # top = [[times[i], times[i+1]] for i in w]*250
  # nmz = [[times[i], times[i+1]] for i in d]*250
  bdom = emptiness([(i[2]-i[1]) for i in bot], sum, 0)
  tdom = emptiness([(i[2]-i[1]) for i in top], sum, 0)
  tnmz = emptiness([(i[2]-i[1]) for i in nmz], sum, 0)
  # bdom = emptiness([250*(i[2]-i[1]) for i in bot], sum, 0)
  # tdom = emptiness([250*(i[2]-i[1]) for i in top], sum, 0)
  # tnmz = emptiness([250*(i[2]-i[1]) for i in nmz], sum, 0)
  return top, tdom, bot, bdom, nmz, tnmz
end


#Get the set of neurons who fired more than ns times in a simulation
#Divide them up according to which pool they came from
#Then you can just calculate statistics from these neurons
#You don't calculate statistics from neurons with too few spikes because the error on those measures is too high
function Neurons_tb_ns(r, half, ns, mini)
  Neurons = collect(Set{Float64}(r))
  f1, f2 = find(Neurons .> half), find(Neurons .<= half)
  N1, N2 = Neurons[f1], Neurons[f2]
  TN, BN = [], [] #top neurons and bottom neurons; who spiked > 50 times
  for i in N1
    if length(find(r .== i)) > ns
      push!(TN, i)
    end
  end

  for i in N2
    if length(find(r .== i)) > ns
      push!(BN, i)
    end
  end


  if ((length(TN) < mini) & (length(BN) < mini))
    return -5, -5
  elseif ((length(TN) < mini) & (length(BN) > mini))
    return -5, convert(Array{Float64}, BN)
  elseif ((length(BN) < mini) & (length(TN) > mini))
    return convert(Array{Float64}, TN), -5
  else
  tn = convert(Array{Float64}, TN)
  bn = convert(Array{Float64}, BN)
  return tn, bn
end
end

#Takes a raster, a set of neurons from above function
#CV and shuff are user options.
#When true, shuff shuffles the neurons so you get a random sampling of stats for nn of them
#CV asks you to calculate coefficient of variation for the neurons as you go
#Main job o fthis function is to bin time and get spike counts for individual neurons in each bin so we can calcualte stats on it (namely FF and correlations)
function count_train(t, r, Neurons, nn, binsize, shuff, CV)
  #t, r, and neurons are determined up front!
  if shuff == true
    tp = shuffle(Neurons)  #in order to randomly sample neurons
  else
    tp = Neurons
  end
  tps = tp[1:nn]
  if CV == true
    cvs = [Float64[0] for i in tps]
  end
  bins = collect(t[1]:binsize:t[end])
  ##############################################################################You fool, you indexed count1 at j, and iterated length(bins)-1 times.
  ##############################################################################common set of 0s is biasing the calculation?
  count1 = zeros(length(this_pool), length(bins)-1)
  for i in eachindex(tps)
    INT = t[find(r.==tps[i])] #time points when neuron i fired
    if CV == true
      disi = diff(INT)
      if length(disi) > 0
        push!(cvs[i], d)
      end
    end
    for j = 1:length(bins)-1
      count1[i,j] = length(find(bins[j] .<= INT .< bins[j+1]))
    end
  end
  if CV == true
    return count1, cvs
  else
    return count1
  end
end

#First you splice flags and return exons
#Then you use ligase to combine into a contiguous time series
function ligase(exons, tom, t, r, Neurons)
  #decide top and bottom business up front, including DDs_1
  # fbomb = zeros(length(Neurons), length(collect(1:fbin:tom))) #fano bins
  # cbomb = zeros(length(Neurons), length(collect(1:cbin:tom))) #corr bins
  #exons*=250
  n1 = minimum(Neurons)
  n2 = maximum(Neurons)
  c = 0
  tf, rf = Float64[0], Float64[0]
  for i in eachindex(exons)
    if i >1
      space = (exons[i][1] - exons[i-1][2])-1
      c+=space
    end
    m = find(exons[i][1] .<= t .<= exons[i][2]) #all activity during this intron
    tm = t[m] - c
    rm = r[m]
    m2 = find(n1 .<= rm .<= n2) #localized to this pool
    rmm = rm[m2]
    tmm = tm[m2]
    tf, rf = cat(1, tf, tmm), cat(1, rf, rmm)
  end
  return tf[2:end], rf[2:end]
end

#Randomly sample n pairwise pearson correlations of count-trains from above (i.e. get a sample of the correlation matrix)
function rand_pair_cor(bin, lt, lr, Neurons, n)
  bins = collect(lt[1]:bin:lt[end])
  ya = [] #if you already looked at this neuron
  neurons = shuffle(Neurons) #
  lank = length(Neurons)
  shank = zeros(lank, length(bins)-1) #store count trains as you go
  cor_store = []
  #loop through n random pairwise correlations
  for i = 1:n
    #pick two random neurons in this pool
    x1 = rand(1:lank)
    x2 = rand(1:lank)
    if x1 == x2 #don't correlate you with yourself (don't worry about the main diagonal of the correlation matrix)
      while x1 == x2
        x2 = rand(1:lank)
      end
    end
    #conditional statements just ensure that you don't waste time re-calculating count trains that you already have
    if ((x1 in ya) & (x2 in ya)) #already calculated count trains for both neurons
      c = cor(vec(shank[x1,:]), vec(shank[x2,:])) #go ahead and correlate them
      if isnan(c) == false
        push!(cor_store, c)
      end
    elseif  ((x1 in ya) & ((x2 in ya) == false)) #if you already did one but not the other
      #update shank for the new neuron
      INT2 = lt[find(lr.==neurons[x2])]
      for j = 1:length(bins)-1
        shank[x2,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      push!(ya, x2)
      if isnan(c) == false
        push!(cor_store, c)
      end
    elseif ((x2 in ya) & ((x1 in ya) == false)) #if you already did the other but not the one
      #update shank for the one neuron
      INT2 = lt[find(lr.==neurons[x1])]
      for j = 1:length(bins)-1
        shank[x1,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      push!(ya, x1)
      if isnan(c) == false
        push!(cor_store, c)
      end
    else #you don't have either
      #get count trains
      INT1 = lt[find(lr.==neurons[x1])]
      INT2 = lt[find(lr.==neurons[x2])]
      #add them to shank so you don't have to recalculate
      for j = 1:length(bins)-1
        shank[x1,j] = length(find(bins[j] .<= INT1 .< bins[j+1]))
        shank[x2,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      #update the list of counts you already have
      push!(ya, x1)
      push!(ya, x2)
    end
  end
  a = convert(Array{Float64}, cor_store)
  #now return the mean of these correlations
  return a
end

#Get count trains within introns or exons
function count_train_intron(bin, lt, lr, Neurons, n, shuff)
  if shuff == true
    neurons = shuffle(Neurons)[1:n]
  else
    neurons = Neurons[1:n]
  end
  bins = collect(lt[1]:bin:lt[end])
  count1 = zeros(n, length(bins)-1)
  for i in eachindex(neurons)
    INT = lt[find(lr.==neurons[i])]
    for j = 1:length(bins)-1
      count1[i,j] = length(find(bins[j] .<= INT .< bins[j+1]))
    end
  end
  return count1
end
#Or just get a raw count train
function count_train(t, r, Neurons, nn, binsize, shuff, CV)
  #t, r, and neurons are determined up front!
  if shuff == true
    tp = shuffle(Neurons)  #in order to randomly sample neurons
  else
    tp = Neurons
  end
  tps = tp[1:nn]
  if CV == true
    cvs = [Float64[0] for i in tps]
  end
  bins = collect(t[1]:binsize:t[end])
  ##############################################################################You fool, you indexed count1 at j, and iterated length(bins)-1 times.
  ##############################################################################common set of 0s is biasing the calculation?
  count1 = zeros(length(this_pool), length(bins)-1)
  for i in eachindex(tps)
    INT = t[find(r.==tps[i])] #time points when neuron i fired
    if CV == true
      disi = diff(INT)
      if length(disi) > 0
        push!(cvs[i], d)
      end
    end
    for j = 1:length(bins)-1
      count1[i,j] = length(find(bins[j] .<= INT .< bins[j+1]))
    end
  end
  if CV == true
    return count1, cvs
  else
    return count1
  end
end

#Pass in an nxc matrix where n is neurons and c is number of time bins
#Get fano factor for each neuron, return descriptive statistics for distribution of Fano factors observed
function fano_train(Fcount, error_code)
  x = []
  for i = 1:size(Fcount)[1]
    f1 = FANO(Fcount[i,:][:])
    if isnan(f1) == false
      push!(x, f1)
    end
  end
  #println(x)
  a = convert(Array{Float64}, x)
  #println(a)
  meanf = emptiness(a, mean, error_code)
  medif = emptiness(a, median, error_code)
  stdvf = emptiness(a, std, error_code)
  return meanf, medif, stdvf
end
#If the simulation measured inputs coming into neurons and the network is rivaling, partition the inputs according to dominance times
function ligate_inputs(introns, input1, input2)
  i1, i2 = Float64[0], Float64[0]
  for i in introns
    s1 = input1[i[1]:i[2]]
    s2 = input2[i[1]:i[2]]
    i1, i2 = cat(1, i1, s1), cat(1, i2, s2)
  end
  s1, s2 = i1[2:end], i2[2:end]
  return s1, s2
end
#Measure inputs to neurons and partition according to dominances when rivaling
#rheobase is hard-coded as 1.
#bal is the t-score meausre: (mean-threshold)/standard deviation
function input_analysis(introns,input1, input2)
  i1, i2 = Float64[0], Float64[0]
  for i in introns
    s1 = input1[i[1]:i[2]]
    s2 = input2[i[1]:i[2]]
    i1, i2 = cat(1, i1, s1), cat(1, i2, s2) #this ends up being length(introns) + 1 time steps longer than it should be
  end
  s1, s2 = i1[2:end], i2[2:end]
  m1, m2 = mean(s1), mean(s2)
  std1, std2 = std(s1), std(s2)

  vs = zeros()

  bal1 = (m1-1.)/std1
  bal2 = (m2-1.)/std2
  return m1, std1, m2, std2, bal1, bal2, s1, s2
end
#Loop over a set of neurons, get the spike times for each, and record the coefficient of variation of inter-spike-intervals
#Store cv_isi for all neurons and return descriptive statistics for distribution of observed cv_isis
function CV_ISI(introns, Neurons, t, r)
  CVS = []
  for i in eachindex(Neurons)
    cvs = []
    INT = t[find(r.==Neurons[i])] #times when neuron i spiked
    for j in eachindex(introns)
      isi = diff(INT[find(introns[j][1] .<= INT .<= introns[j][2])])
      if length(isi) > 0
        for d in isi
          push!(cvs, d)
        end
      end
    end
    cvc = convert(Array{Float64}, cvs)
    cov = emptiness(cvc, cv, -5)
    if (isnan(cov) == false) & (cov != -5)
      push!(CVS, cov)
    end
  end
  a = convert(Array{Float64}, CVS)
  mcv = emptiness(a, mean, -5)
  medcv = emptiness(a, median, -5)
  stdcv = emptiness(a, std, -5)
  return mcv, medcv, stdcv
end
#Same, but instead of returning descriptive stats (as needed for parameter scanning)
#Just return the whole distribution of cv_isi so you can study it in a local simulation
function CV_ISI_D(introns, Neurons, t, r)
  CVS = []
  for i in eachindex(Neurons)
    cvs = []
    INT = t[find(r.==Neurons[i])] #times when neuron i spiked
    for j in eachindex(introns)
      isi = diff(INT[find(introns[j][1] .<= INT .<= introns[j][2])])
      if length(isi) > 0
        for d in isi
          push!(cvs, d)
        end
      end
    end
    cvc = convert(Array{Float64}, cvs)
    cov = emptiness(cvc, cv, -5)
    if (isnan(cov) == false) & (cov != -5)
      push!(CVS, cov)
    end
  end
  a = convert(Array{Float64}, CVS)
  return a
end
#same thing, but don't worry about isolating dominances. Just get raw cvisi
function CV_ISI_ALLTIME(Neurons, t, r)
  CVS = []
  for i in eachindex(Neurons)
    INT = t[find(r.==Neurons[i])]
    isi = diff(INT)
    c = cv(isi)
    if isnan(c) == false
      push!(CVS, c)
    end
  end
  mcv = emptiness(CVS, mean, -5)
  x = convert(Array{Float64}, CVS)
  return CVS
end

#Get mean, std, skew, kurtosis of a signal x
function moments(x)
  n = length(x)
  #correction = sqrt(n*(n-1))/(n-2)
  m = mean(x)
  me = median(x)
  skew = zeros(n)
  kurt = zeros(n)
  vary = zeros(n)
  for i = 1:n
    d = m-x[i]
    skew[i] = d^3
    kurt[i] = d^4
    vary[i] = d^2
  end
  ms = mean(skew)
  mk = mean(kurt)
  mv = mean(vary)
  cmv = sum(vary)/(n-1)
  skewness = ms/(mv^1.5)
  #sample_skewness = skewness*correction
  kurtosis = mk/(mv^2)
  return m,cmv, skewness, kurtosis
end
