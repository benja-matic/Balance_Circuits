###Various and sundry functions I may have thrown out
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

function adapt_switch(adapt_top, adapt_bot, th, tl)
  aa = adapt_bot .+ adapt_top
  at = adapt_top ./ aa
  at = at[2:end]
  flags, times = WLD_A(at, th, tl)
  return flags, times
end


# function low_pass_DD(times, mint)
#   #mint is the threshold for reversions; if shorter than mint, we will assume dominance throughout
#   xtime = []
#   i = 1
#   while i < (length(times) -1)
#     c = i
#     while (((times[c+1][1] - times[c][2]) < mint) & (c <= length(times) -1))
#       println(times[c+1][1] - times[c][2], " and ", (times[c+1][1] - times[c][2]) < mint, " and ", c, " and ", i)
#       c +=1
#     end
#     println("just exited while loop")
#     push!(xtime, [times[i][1], times[c][2]])
#     i += c
#   end
#   return xtime
# end

  # xtime = [times[i+1][1] - times[i][2] for i = 1:length(times)-1] #times between dominances
  # revs = [] #index where there's < mint time between your end and next start
  # for i in eachindex(xtime)
  #   if xtime[i] < mint
  #     push!(revs, i)
  #   end
  # end
  # dx = diff(revs)
  # x2 = []
  # for i in eachindex(dx)
  #   if dx[i] == 1
  #     push!(x2, i)
  #   end
  # end
  #
  # ftime = []
  # dtime = []
  # for i = 1:length(times) - 1
  #   if i in revs
  #     push!(ftime, (times[i][1], times[i+1][2]))
  #     push!(dtime, i+1)
  #   elseif
  #
  #
  #
  # for i = 1:length(times)-1
  #   if (times[i+1][1] - times[i][2]) < mint
  #     push!(xtime, (times[i][1], times[i+1][2])
  #   end
  # end



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

# function low_pass_DD(flags, times, mint)
#   xf = []
#   xt = []
#   l = find(flags.=="lose")
#   w = find(flags.=="win")
#   d = find(flags.=="draw")
#   lm = minimum(l)
#   ld = minimum(w)
#   m = minimum([lm, ld])
#   fx = flags[m:end]
#   tx = times[m:end]
#   d = diff(tx)
#   dx = find(d .>= mint)
#   return fx[dx], tx[dx]
# end

function low_pass_DD(flags, times, mint)
  wflag = true
  tn, fn, wflag = diff_eater(times, flags, mint)
  while wflag == true
    tn, fn, wflag = diff_eater(tn, fn, mint)
  end
  return tn, fn
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


#given count trains for a set of neurons in two pools, get pearson correlation for all N^2 pairs
#I abandoned this in favor of the random sampling technique
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

#Given two input time series, get some statistics
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

function analyze_ligated_inputs(input)
  m = mean(input)
  s = std(input)
  t = (m-1.)/s
  return m, s, t
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

#Instead of convolving two signals, circularly shift one
#Get pearson correlation at each lag
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

#forget if this actually does anything useful
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

#take two arrays of the same length
#shift one by a time lag tau, padding with zeros
#get the pearson correlation of the non-zero padded parts
#do this for time lags -empezar:empezar

function tau_cor(s1, s2, empezar)
  ea = abs(empezar)
  n = ea*2
  cc = zeros(n+1) #correlate with 0 lag, and +/- ea
  cc[ea+1] = cor(s1, s2) #julia can't index at zero so have to do this manually
  for i = 2:ea+1
    k1 = cor(s1[i:end][:], s2[1:end-(i-1)])
    k2 = cor(s2[i:end][:], s1[1:end-(i-1)])
    cc[ea+2-i] = k1
    cc[ea+i] = k2
  end
  return cc
end

function tau_cov(s1, s2, empezar)
  ea = abs(empezar)
  n = ea*2
  cc = zeros(n+1) #correlate with 0 lag, and +/- ea
  cc[ea+1] = dot(s1, s2)/length(s1) #julia can't index at zero so have to do this manually
  for i = 2:ea+1
    k1 = dot(s1[i:end][:], s2[1:end-(i-1)])/(length(s1)-(i-1))
    k2 = dot(s2[i:end][:], s1[1:end-(i-1)])/(length(s1)-(i-1))
    cc[ea+2-i] = k1
    cc[ea+i] = k2
  end
  return cc
end

#Forget if this does anything
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

#Given a signal s, measure mean over an interval 'window', and increment by 'shift'
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
