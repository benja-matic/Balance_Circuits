function von_mises_dist(x, k, mu, N)
  a = exp(k*cos(x-mu))/(N*besseli(0, k))
end

#

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

function euler_lif_NA(h, total, Ne, Ni, wee, wei, wie, wii, s1, s2, vth, tau_m, tau_ee, tau_ei, tau_ie, tau_ii)

  ntotal = round(Int,total/h)
  srand(1234)

  #initial conditions
  ve = rand(Ne)*vth
  vi = rand(Ni)*vth
  ve_buff = ve
  vi_buff = vi
  e_top = zeros(ntotal)
  e_bot = zeros(ntotal)
  i_guy = zeros(ntotal)
  # ev = zeros(ntotal)

  #storage space
  see = zeros(Ne)#synapse
  sei = zeros(Ne)
  sie = zeros(Ni)
  sii = zeros(Ni)
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
  #Set up drive to two pools special for competitive network
  drive = zeros(Ne)
  FPe = div(Ne,5) #each pool is 2/5 of the network
  P1s = round(Int,Ne/4-FPe)
  P1e = round(Int,Ne/4+FPe)
  P2s = round(Int,3*Ne/4-FPe)
  P2e = round(Int,3*Ne/4+FPe)
  drive[P1s:P1e] = s1*h
  drive[P2s:P2e] = s2*h
  half = div(Ne, 2)
  quarter = div(half, 2)
  rt = half + quarter
  rb = quarter
  ig = div(Ni, 2)

  kill_flag = false
  for iter = 1:ntotal
    #administer drive, leak, synapse, and adaptation to voltage
    ve += drive + (h*see) + (h*sei) - (ve*m_leak) #- (ae*g_e)
    vi += (h*sie) + (h*sii) -(vi*m_leak) #- (ai*g_i)
    #ev[iter] = ve[rt]
    e_top[iter] = (h*see[rt]) + (h*sei[rt])
    e_bot[iter] = (h*see[rb]) + (h*sei[rb])
    #e_top[iter] = drive[rt] + (h*see[rt]) + (h*sei[rt]) #- (ve[rt]*m_leak) #- (ae[rt]*g_e)
    #e_bot[iter] = drive[rb] + (h*see[rb]) + (h*sei[rb]) #- (ve[rb]*m_leak) #- (ae[rb]*g_e)
    i_guy[iter] = (h*sie[ig]) + (h*sii[ig]) #-(vi[ig]*m_leak)
    #administer leak to synapse
    see -= (see*ee_leak)
    sei -= (sei*ei_leak)
    sie -= (sie*ie_leak)
    sii -= (sii*ii_leak)
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
      if (length(raster_e)/Ne > 200*h*iter*(1/1000)) | (length(raster_i)/Ni > 200*h*iter*(1/1000))
        kill_flag = true
        break
      end
    end
  end



  return time_e, raster_e, time_i, raster_i, kill_flag, e_top, e_bot, i_guy
end

function write_result(wta, etm, etv, ebm, ebv, cor_within, cor_between, ntf, mif, mcv, lte, id)
  println("##RESULT $(wta), $(etm), $(etv), $(ebm), $(ebv), $(cor_within), $(cor_between), $(ntf), $(mif), $(mcv), $(lte), $(id)")
end

function emptiness(x, funk, error_code) #checks if an empty set got populated
  if length(x) > 0
    return funk(x)
  else
    return error_code
  end
end

function CV_ISI_ALLTIME(Neurons, t, r)
  if Neurons == -5
    return -5, -5, -5
  else
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
end

function score_analysis(re, N)
  #1 = top, 2 = bot
  half = div(N, 2)
  a1 = length(find(re .<= half))
  a2 = length(find(re .> half))
  biggie = max(a1, a2)
  spike_asymmetry = biggie/length(re)
  if a1 > a2
    bias = 2
  else
    bias = 1
  end
  return spike_asymmetry, bias
end

function network_fano(t, r, bin_size, ntotal)

  ntf_bins = collect(minimum(t):bin_size:ntotal)
  ntf_bin_num = length(ntf_bins)
  ntf_nt = zeros(ntf_bin_num-1)
  nt = zeros(length(ntf_nt))

  for j = 2:ntf_bin_num
    nt[j-1] = sum(ntf_bins[j-1] .<= t .< ntf_bins[j]) #network counts
  end

  fano = FANO(nt)
  return fano
end

function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  ###x = dv/dt
  ###vth = v2 - xt
  ###(vth - v2)/x = t
  return t
end

function cv(isi)
  SD = std(isi)
  return SD/mean(isi)
end

function FANO(aRAY)
  return var(aRAY)/mean(aRAY)
end

function analyze_inputs(s)
  m = mean(s)
  SD = std(s)
  t = (m-1.)/SD
  return m, SD, t
end


function Neuron_finder(r, ns, mini)
  Neurons = collect(Set{Float64}(r))
  N = []
  #S = []
  for i in Neurons
    a = find(r .== i)
    if length(a) > ns
      push!(N, i)
      #push!(S, length(a))
    end
  end
  if length(N) < mini
    NF = -5
    #SF = -5
  else
    NF = convert(Array{Float64}, N)
    #SF = convert(Array{Int64}, S)
  end
  return NF
end

function Neurons_tb_ns(r, half, ns, mini)
  Neurons = collect(Set{Float64}(r))
  f1, f2 = find(Neurons .> half), find(Neurons .<= half)
  N1, N2 = Neurons[f1], Neurons[f2]
  TN, BN = [], [] #top neurons and bottom neurons; who spiked > ns times
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

function count_train(bin, lt, lr, Neurons, n)
  if Neurons == -5
    return -5
  else
  bins = collect(lt[1]:bin:lt[end])
  count1 = zeros(n, length(bins)-1)
  for i in eachindex(Neurons)
    INT = lt[find(lr.==Neurons[i])]
    for j = 1:length(bins)-1
      count1[i,j] = length(find(bins[j] .<= INT .< bins[j+1]))
    end
  end
  return count1
end
end

function fano_train(Fcount, error_code)
  if Fcount == -5
    return -5, -5, -5
  else
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
end

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
    if x1 == x2
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

function moments(x)
  n = length(x)
  m = mean(x)
  skew = zeros(n)
  kurt = zeros(n)
  vary = zeros(n)
  for i =1:n
    d = x[i] - m
    skew[i] = d^3
    kurt[i] = d^4
    vary[i] = d^2
  end
  ms = mean(skew)
  mk = mean(kurt)
  mv = mean(vary)
  skewt = ms/(mv^1.5)
  KURT = mk/(mv^2)
  t = (m-1.)/sqrt(mv)
  return m, mv, t, skewt, KURT
end
