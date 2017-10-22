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

function euler_lif(h, total, Ne, Ni, wee, wei, wie, wii, sE, sI, vth, tau_m, tau_ee, tau_ei, tau_ie, tau_ii, tau_ae, tau_ai, g_e, g_i)

  ntotal = round(Int,total/h)
  srand(1234)

  #initial conditions
  ve = rand(Ne)*vth
  vi = rand(Ni)*vth
  ve_buff = ve
  vi_buff = vi
  ishuff = shuffle!(collect(1:Ni))
  eshuff = shuffle!(collect(1:Ne))

  e_input = zeros(Ne, ntotal)
  e_half = div(Ne, 2)
  i_half = div(Ni, 2)

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
  drive = fill(sE, Ne)
  drive_i = fill(sI, Ni)

  kill_flag = false
  for iter = 1:ntotal
    #administer drive, leak, synapse, and adaptation to voltage
    ve += drive + (h*see) + (h*sei) - (ve*m_leak) - (ae*g_e)
    vi += drive_i + (h*sie) + (h*sii) -(vi*m_leak) - (ai*g_i)
    e_input[:,iter] = drive + (h*see) + (h*sei) - (ve*m_leak)#

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
        #push!(time_e, iter)
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
        #push!(time_i, iter)
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
      if (length(raster_e)*(1./Ne) > 20) | (length(raster_i)*(1./Ni) > 20) #Assume h=.1ms; 20 spikes per neuron in 100 ms
        kill_flag = true
        break
      end
    end
  end



  return time_e, raster_e, time_i, raster_i, kill_flag, e_input
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

function Neuron_finder(r, ns, mini)
  Neurons = collect(Set{Float64}(r))
  N = []
  S = []
  for i in Neurons
    a = find(r .== i)
    if length(a) > ns
      push!(N, i)
      push!(S, length(a))
    end
  end
  if length(N) < mini
    NF = -5
    SF = -5
  else
    NF = convert(Array{Float64}, N)
    SF = convert(Array{Int64}, S)
  end
  return NF, SF
end


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

function heeger_dominance(ntd, ntotal)
  shmitt_thresh = .7
  ntd_post_transient = ntd[10:end]# this kills initial transient
  intdices = []
  if maximum(ntd) < -.3
    return -1, -1
  elseif minimum(ntd) > .3
    return -1, -1
  elseif (-.7 <= maximum(ntd) <= .7) & (-.7 <= minimum(ntd) <= .7)
    return -2, -2
  else
    i1t = find(ntd_post_transient .> shmitt_thresh)
    i1b = find(ntd_post_transient .< -shmitt_thresh)
    # println(length(i1t))
    # println(length(i1b))
    if isempty(i1t) & isempty(i1b) #crosses 0, but never crosses threshold. This can be due to 2 enourmouse synchronous fluctuations elevating the threshold, for example.
      return -3, -3
    elseif isempty(i1t) # crosses threshold, and one is empty but the other isn't. Definitely WTA.
      i1s = i1b[1]
      c_thresh = shmitt_thresh
      flag = false
      u_score = ntd[i1s:end]
      u_sum = ntd[i1s:end]
    elseif isempty(i1b) #crosses threshold, and one is empty but the other isn't. Definitely WTA.
      i1s = i1t[1]
      c_thresh = -shmitt_thresh
      flag = true
      u_score = ntd[i1s:end]
      u_sum = ntd[i1s:end]
    else
      i1s = minimum([i1t[1], i1b[1]])#who gets the first dominance duration
      u_score = ntd[i1s:end] #start measuring here
      u_sum = ntd[i1s:end]
      if i1s == i1t[1]
        c_thresh = -shmitt_thresh
        flag = true
      else
        c_thresh = shmitt_thresh
        flag = false
      end
    end
    top_n_bottom = [flag]
    empezamos = (i1s*250)
    push!(intdices, empezamos)
    #push!(top_n_bottom, flag)
    for i in eachindex(u_score)
      if (u_score[i] < c_thresh) == flag
        c_thresh *= -1
        flag = !flag
        push!(intdices, empezamos + (i*250)) #this is the time point of a transition starting after the transient_end time
        push!(top_n_bottom, flag)
      end
    end

  if length(intdices) == 0
    return 0, 0
  else
    push!(intdices, ntotal)
    push!(top_n_bottom, !flag)
  end #times of transitions, who is winning for each time, and indices in ntd for transitions
    return intdices, top_n_bottom
  end
end

function comparison(flag, x)
  if flag == "win"
    c = x >= .3
    if c == true
      return flag
    else
      return "draw"
    end
  elseif flag == "lose"
    c = x <= -.3
    if c == true
      return flag
    else
      return "draw"
    end
  elseif flag == "draw"
    c = -.3 <= x <= .3
    if x >.7
      return "win"
    elseif x < -.7
      return "lose"
    else
      return flag
    end
  elseif flag == "almost"
    if x > .7
      return "win"
    elseif x < -.7
      return "lose"
    elseif -.3 <= x <= .3
      return "draw"
    else
      return flag
    end
  else
    if x > .7
      return "win"
    elseif x < -.7
      return "lose"
    elseif -.3 <= x <= .3
      return "draw"
    elseif .3 < abs(s[1]) < .7
      flag = "almost"
    else
      return flag
    end
  end
end

function comp_01(x)
  if x > .7
    return "win"
  elseif .3 <= x <= .7
    return "draw"
  elseif x < .3
    return "lose"
  else
    return "weird"
  end
end

function win_lose_draw(s)
  if maximum(s) < -.3
    #return -1, -1 ############################################################## highly asymmetric=wta
    return ["lose", "end"], [1, 799]
  elseif minimum(s) > .3
    #return -1, -1 ############################################################## highly asymmetric=wta
    return ["win", "end"], [1, 799]
  elseif (-.7 <= maximum(s) <= .7) & (-.7 <= minimum(s) <= .7)
    #return -2, -2 ############################################################## normalization
    return ["draw", "end"], [1, 799]
  end
  times = []
  flags = []
  if s[1] >.7
    flag = "win"
  elseif -3. <= s[1] <= .3
    flag = "draw"
  elseif s[1] < -.7
    flag = "lose"
  elseif .3 < abs(s[1]) < .7
    flag = "almost"
  else
    flag = "weird"
  end
  push!(times, 1)
  push!(flags, flag)
  s2 = s[2:end]
  for i in eachindex(s2)
    f1 = comparison(flag, s2[i])
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

function WLD_01(s)
  if maximum(s) < .3
    return ["lose", "end"], [1, 799]
  elseif minimum(s) >= .7
    return ["win", "end"], [1, 799]
  elseif (.3 <= maximum(s) <= .7) & (.3 <= minimum(s) <= .7)
    return ["draw", "end"], [1, 799]
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
    f1 = comp_01(s2[i])
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




function fatigue_switches(aet, aeb, tauadapt, h)
  aeb_pt = aeb[Int64(tauadapt*(3/h)):end];
  aet_pt = aet[Int64(tauadapt*(3/h)):end];
  if aeb_pt[1] > aet_pt[1]
    flag = false
  else aet_pt[1] > aeb_pt[1] #assuming fatigue variables are never exactly equal in else statement
    flag = true
  end
  transition_times = []
  for i = 2:length(aeb_pt)
    now = aet_pt[i] > aeb_pt[i] #true if top more fatigued, else false
    if now != flag # if
      flag = !flag
      push!(transition_times, i)
    end
  end
  a = convert(Array{Float64}, transition_times)
  m = mean(aet_pt[a])
  return a, m
end

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

function local_peaks(signal1, signal2, boundaries)
  bnds = [1]
  for i in boundaries
    push!(bnds, i)
  end
  z = length(boundaries)
  s1ma, s1mi = zeros(z), zeros(z)
  s2ma, s2mi = zeros(z), zeros(z)
  for i = 1:z
    s1 = signal1[bnds[i]:bnds[i+1]]
    s2 = signal2[bnds[i]:bnds[i+1]]
    s1ma[i], s1mi[i] = indmax(s1)+bnds[i], indmin(s1)+bnds[i]
    s2ma[i], s2mi[i] = indmax(s2)+bnds[i], indmin(s2)+bnds[i]
  end
  return s1ma, s1mi, s2ma, s2mi
end

function really_rivaling(intdices, tnb, tdom, bdom, ntotal, s)
  a = diff(intdices)
  x = []
  inds = div(intdices, 250)
  if length(intdices) > 32
    return 0
  end
  # if inds[end] == length(s+1)
  #   println("last index rectified")
  #   inds[end] = length(s)
  # end
  # if inds[1] == 0
  #   prinln("first index rectified")
  #   inds[1] +=1
  # end
  inds[end] -=1
  asym = zeros(length(inds)-1)
  for i = 1:length(inds) -1
    asym[i] = abs(mean(s[inds[i]:inds[i+1]]))
  end
  if maximum(a) < 6000
    return 0
  elseif maximum(asym) < .3
    return 0
  else
    for i in eachindex(a)
      if ((a[i]>6000) & (asym[i] > .5))
        push!(x, tnb[i])
      end
    end
  end
  if length(Set(x)) == 1
    if minimum([bdom, tdom])/ntotal < .2
      return 2
    else
      return 1
    end
  else
    return 1
  end
end

function shmitt_on_ad(ntd, range_frac)
  range_score = maximum(ntd) - minimum(ntd)
  shmitt_thresh = range_score/range_frac
  m = mean(ntd)
  thresh = [shmitt_thresh+m, shmitt_thresh-m]
  println("$(m), is the mean of this signal")
  intdices = []
  i1t = find(ntd .> m+shmitt_thresh)
  i1b = find(ntd .< m-shmitt_thresh)
  if isempty(i1t) & isempty(i1b)
    return -3, -3
  elseif isempty(i1t)
    i1s = i1b[1]
    ci = 1
    c_thresh = thresh[ci]
    flag = false
    u_score = ntd[i1s:end]
  elseif isempty(i1b)
    i1s = i1t[1]
    ci=2
    c_thresh = thresh[ci]
    flag = true
    u_score = ntd[i1s:end]
  else
    i1s = minimum([i1t[1], i1b[1]])
    u_score = ntd[i1s:end]
    if i1s == i1t[1]
      ci = 1
      c_thresh = thresh[ci]
      flag = true
    else
      ci = 2
      c_thresh = thresh[ci]
      flag = false
    end
  end
  top_n_bottom = [flag]
  empezamos = i1s
  push!(intdices, empezamos)
  #push!(top_n_bottom, flag)
  for i in eachindex(u_score)
    if (u_score[i] < c_thresh) == flag
      thresh = circshift(thresh, 1)
      c_thresh = thresh[ci]
      flag = !flag
      push!(intdices, empezamos + i) #this is the time point of a transition starting after the transient_end time
      push!(top_n_bottom, flag)
    end
  end

  if length(intdices) == 0
    return 0, 0
  end
  if ntotal - intdices[end] > 10000 #longer than 1 second
    flag = !flag
    push!(intdices, ntotal)
    push!(top_n_bottom, flag)
  end #times of transitions, who is winning for each time, and indices in ntd for transitions
    return intdices, top_n_bottom
  end

function interpolate_asymmetry(DDs_1, t, r)
  dspike = zeros(length(DDs_1)-1)
  for i = 1:length(DDs_1)-1
    x = find(DDs_1[i] .<= t .< DDs_1[i+1]) #mask for spike times within dominance
    #tx = te[x]
    rx = re[x] #Neurons who spiked during dominance
    s1 = sum(rx .> half)
    s2 = sum(rx .<= half)
    wta_ness = max(s1, s2)/length(rx)
    dspike[i] = wta_ness
  end
  dtimes = diff(DDs_1)
  return dspike, dtimes
end

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

function rand_pair_cor(bin, lt, lr, Neurons, n)
  bins = collect(lt[1]:bin:lt[end])
  count1 = zeros(n, length(bins)-1)
  ya = []
  neurons = shuffle(Neurons)
  lank = length(Neurons)
  shank = zeros(lank, length(bins)-1)
  cor_store = []
  #loop through however many sample correlations you want
  for i = 1:n
    #pick two random neurons in this pool
    x1 = rand(1:lank)
    x2 = rand(1:lank)
    #conditional statements just ensure that you don't waste time re-calculating count trains that you already have
    if ((x1 in ya) & (x2 in ya))
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      if isnan(c) == false
        push!(cor_store, c)
      end
    elseif  ((x1 in ya) & ((x2 in ya) == false))
      INT2 = lt[find(lr.==neurons[x2])]
      for j = 1:length(bins)-1
        shank[x2,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      push!(ya, x2)
      if isnan(c) == false
        push!(cor_store, c)
      end
    elseif ((x2 in ya) & ((x1 in ya) == false))
      INT2 = lt[find(lr.==neurons[x1])]
      for j = 1:length(bins)-1
        shank[x1,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      push!(ya, x1)
      if isnan(c) == false
        push!(cor_store, c)
      end
    else
      #get count trains
      INT1 = lt[find(lr.==neurons[x1])]
      INT2 = lt[find(lr.==neurons[x2])]
      #add them to the counts matrix so you don't have to recalculate
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
  return mean(a)
end


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

function correlate_between(count1, count2)
  #println(size(count1))
  #println(size(count2))
  bc = zeros(size(count1)[1]*size(count2)[1])
  c=0
  p1 = size(count1)[1]
  p2 = size(count2)[1]
  #println(p1)
  #println(p2)
  for i=1:p1
    for j = 1:p2
      c1, c2 = collect(count1[i,:]), collect(count2[j,:])
      corre = cor(c1, c2)
      c+=1
      bc[c] = corre
    end
  end
  ############################################################################Handle NaN issues properly
  cb = []
  for i in bc
    if isnan(i) == false
      push!(cb, i)
    end
  end
  cor_between = emptiness(cb, mean, error_code)
  return cor_between
end

function fano_train(Fcount, error_code)
  x = []
  for i = 1:size(Fcount)[1]
    f1 = FANO(Fcount[i,:])
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

function input_measures(s1, s2, vbinsize)
  m1 = mean(s1)
  m2 = mean(s2)
  std1 = std(s1)
  std2 = std(s2)
  bn = div(length(s1), vbinsize)
  mvs = zeros(6, bn)
  c=1
  for i = 1:bn
    chunk1 = s1[c:c+vbinsize]
    chunk2 = s2[c:c+vbinsize]
    ms1 = mean(chunk1)
    ms2 = mean(chunk2)
    vs1 = var(chunk1)
    vs2 = var(chunk2)
    ts1 = tscore(ms1, std(chunk1))
    ts2 = tscore(ms2, std(chunk2))
    mvs[1,i] = mean(chunk1)
    mvs[2,i] = mean(chunk2)
    mvs[3,i] = var(chunk1)
    mvs[4,i] = var(chunk2)
    mvs[5,i] = ts1
    mvs[6,i] = ts2
    c+=vbinsize
  end
  #fano factor
  #VOMIT
  #variance of the variance
  #tscore in time windows
  mtsc1 = mean(mvs[5,:]) #mean tscore across windows
  mtsc2 = mean(mvs[6,:]) #same
  VOMIT1 = var(mvs[1,:]) #variane of the mean
  VOMIT2 = var(mvs[2,:]) #variance of the mean
  f1 = VOMIT1/m1
  f2 = VOMIT2/m2
  MOVIT1 = mean(mvs[3,:])
  MOVIT2 = mean(mvs[4,:])
  VOVIT1 = var(mvs[3,:])
  VOVIT2 = var(mvs[4,:])
  return mtsc1, mtsc2, VOMIT1, VOMIT2, f1, f2, MOVIT1, MOVIT2, VOVIT1, VOVIT2, m1, m2, std1, std2
end





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

  bal1 = (20.-m1)/std1
  bal2 = (20.-m2)/std2
  return m1, std1, m2, std2, bal1, bal2, s1, s2
end

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
  CVS = convert(Array{Float64}, CVS)
  mean_cv = emptiness(CVS, mean, -5)
  median_cv = emptiness(CVS, median, -5)
  std_cv = emptiness(CVS, std, -5)
  return mean_cv, median_cv, std_cv
end

function correlate_within(counts, error_code)
  l = size(counts)[1]
  tamano = binomial(l, 2)
  cores = zeros(tamano)
  c=0
  for i=1:l -1
    for j =i+1:l
      c1, c2 = collect(counts[i,:]), collect(counts[j,:])
      #println(cor(c1, c2))
      c+=1
      cores[c] = cor(c1, c2)
    end
  end
  cn = []
  for i in cores
    if isnan(i) == false
      push!(cn, i)
    end
  end
  cor_within = emptiness(cn, mean, error_code)
  return cor_within
end

function correlate_between(count1, count2)
  #println(size(count1))
  #println(size(count2))
  bc = zeros(size(count1)[1]*size(count2)[1])
  c=0
  p1 = size(count1)[1]
  p2 = size(count2)[1]
  #println(p1)
  #println(p2)
  for i=1:p1
    for j = 1:p2
      c1, c2 = collect(count1[i,:]), collect(count2[j,:])
      corre = cor(c1, c2)
      c+=1
      bc[c] = corre
    end
  end
  ############################################################################Handle NaN issues properly
  cb = []
  for i in bc
    if isnan(i) == false
      push!(cb, i)
    end
  end
  cor_between = emptiness(cb, mean, -5)
  return cor_between
end

function real_dominances(dspike, dtimes, tnb)
  RDT = []
  tbi = []
  for i in eachindex(dspike)
    if ((dtimes[i] > 6000) & (dspike[i] > .6))
      push!(RDT, i)
      push!(tbi, tnb[i])
    end
  end
  return RDT, tbi
end

function circ_cor(s1, s2)
  if length(s1) != length(s2)
    return "signals not same length"
  else
    s1 = reshape(s1, length(s1))
    s2 = reshape(s2, length(s2))
    corz = zeros(length(s1))
    for i = 1:length(s1)
      sn = circshift(s1, i)
      c = cor(vec(s2), vec(sn))
      corz[i] = c
    end
    return corz
  end
end


function sliding_correlation(s1, s2)
  if length(s1) != length(s2)
    return "signals not same length"
  else
    #2 through half, half+1 through end
    # s1 = convert(Array{Float64}, s1)
    # s2 = convert(Array{Float64},(length(s2),) s2)
    largo = length(s1)
    l2 = div(length(s1), 2) - 1
    cpos = zeros(l2)
    cneg = zeros(l2)
    println("looping")
    for i = 1:l2
      sp = circshift(s1, i)
      sn = circshift(s1, largo - i)
      corp = dot(vec(s2), vec(sp))
      corn = dot(vec(s2), vec(sn))
      cpos[i] = corp
      cneg[i] = corn
    end
    aa = cat(1, cneg, 1)
    a = cat(1, aa, cpos)
    return a
  end
end

# function sliding_dot(s1, s2)
#   z = minimum(length(s1), length(s2))

function moments(x)
  n = length(x)
  correction = sqrt(n*(n-1))/(n-2)
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
  sample_skewness = skewness*correction
  kurtosis = mk/(mv^2)
  return m, me, cmv, skewness, sample_skewness, kurtosis
end

function basic_moments(x)
  n = length(x)
  correction = sqrt(n*(n-1))/(n-2)
  m = mean(x)
  #me = median(x)
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
  stdv = sqrt(sum(vary)/(n-1))
  skewness = ms/(mv^1.5)
  sample_skewness = skewness*correction
  kurtosis = mk/(mv^2)
  return m, stdv, skewness, kurtosis
end


function fraction_quiescent(N, r)
  q = []
  for i =1:N
    f = find(r.==i)
    if length(f) < 1
      push!(q, i)
    end
  end
  return q
end

function rate_dist(r)
  s = Set(r)
  a = zeros(length(s))
  for i = 1:length(s)
    f = find(r.==s[i])
    a[i] = length(f)
  end
  return a
end


function network_sps(ttf)
  s = ttf[1]
  xt = ttf[end] - s
  xti = Int64(xt)
  nsps = zeros(xti)
  for i = 1:xti
    a = length(find(ttf.==i + s -1))
    nsps[i] = a
  end
  return nsps
end

function running_mean(s, window, shift)
  x = []
  c=1
  w = window-1
  while c < length(s) - window
    a = mean(s[c:c+w])
    push!(x, a)
    c+=shift
  end
  return x
end
