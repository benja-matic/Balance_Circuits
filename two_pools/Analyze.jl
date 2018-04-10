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

function nt_diff_angle(t, r, ntotal, P2s, P2e, netd_binsize)

  netd_bins = collect(1:netd_binsize:ntotal)
  ntd = zeros(length(netd_bins)-1)
  nts = zeros(length(netd_bins)-1)

  for j = 2:length(netd_bins)
    tf = find(netd_bins[j-1] .<= t .< netd_bins[j])
    T = sum(P2s .< r[tf] .< P2e)
    B = length(tf) - T
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
function splice_flags(flags, times, netd_binsize)
  l = find(flags.=="lose")
  w = find(flags.=="win")
  d = find(flags.=="draw")
  bot = [[((times[i]-1) .* netd_binsize)+1, ((times[i+1]-1) .* netd_binsize)] for i in l]
  top = [[((times[i]-1) .* netd_binsize)+1, ((times[i+1]-1) .* netd_binsize)] for i in w]
  nmz = [[((times[i]-1) .* netd_binsize)+1, ((times[i+1]-1) .* netd_binsize)] for i in d]
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
  if Neurons == -5
    return -5
  end
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
  if Fcount == -5
    return -5
  end
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
  return a
end
#If the simulation measured inputs coming into neurons and the network is rivaling, partition the inputs according to dominance times
function ligate_inputs(introns, input1, input2)
  i1, i2 = Float64[0], Float64[0]
  for i in introns
    s1 = input1[i[1]:i[2]]
    s2 = input2[i[1]:i[2]]
    i1 = cat(1, i1, s1)
    i2 = cat(1, i2, s2)
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
  # mcv = emptiness(a, mean, -5)
  # medcv = emptiness(a, median, -5)
  # stdcv = emptiness(a, std, -5)
  return a
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
  if Neurons == -5
    return -5
  end
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
  return x
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

function WTAN_Analysis(t, r, Input, end_trans, Ne, Ni, half_e, min_e_neurons, min_i_neurons, fbinsize, cbinsize, FPe)

  P1s = round(Int,Ne/4-FPe)
  P1e = round(Int,Ne/4+FPe)
  P2s = round(Int,3*Ne/4-FPe)
  P2e = round(Int,3*Ne/4+FPe)

  P1 = P1s:P1e
  P2 = P2s:P2e

  rt = ((ntotal - end_trans)/1000.)*h

  ex = find(r .<= Ne)
  te = t[ex]
  re = r[ex]

  ix = find(r .> Ne)
  ti = t[ix]
  ri = r[ix]

  E_Input = Input[1:Ne, :]
  I_Input = Input[Ne+1:end, :]

  E_bot = E_Input[P1s:P1e, :]
  E_top = E_Input[P2s:P2e, :]

  e_top_pt = E_top[:, end_trans:end]
  e_bot_pt = E_bot[:, end_trans:end]
  i_pop_pt = I_Input[:, end_trans:end]

  tem = find(te.> end_trans)
  te_pt = te[tem]
  re_pt = re[tem]
  tim = find(ti.> end_trans)
  ti_pt = ti[tim]
  ri_pt = ri[tim]

  wta_ness, bias = score_analysis(re_pt, Ne)
  top_e_neurons, bot_e_neurons = Neurons_tb_ns(re_pt, half_e, 10, min_e_neurons)
  if ((top_e_neurons == -5) & (bot_e_neurons == -5))
    return "too few neurons spiking for analysis"
  else
  I_Neurons = Neuron_finder(ri_pt, 10, min_i_neurons)
  #rates
  E_Rates = [length(find(re_pt .== i))/rt for i=1:Ne]
  I_Rates = [length(find(ri_pt .== i))/rt for i=1:Ni]

  #counts
  E_count_top = count_train_intron(fbinsize, te_pt, re_pt, top_e_neurons, length(top_e_neurons), false)
  E_count_bot = count_train_intron(fbinsize, te_pt, re_pt, bot_e_neurons, length(bot_e_neurons), false)
  I_count_all = count_train_intron(fbinsize, ti_pt, ri_pt, I_Neurons, length(I_Neurons), false)
  #FF
  FF_TOP = fano_train(E_count_top, -5)
  FF_BOT = fano_train(E_count_bot, -5)
  FF_INH = fano_train(I_count_all, -5)
  #cv
  CV_TOP = CV_ISI_ALLTIME(top_e_neurons, te_pt, re_pt)
  CV_BOT = CV_ISI_ALLTIME(bot_e_neurons, te_pt, re_pt)
  CV_INH = CV_ISI_ALLTIME(I_Neurons, ti_pt, ri_pt)
  #synchrony
  RSC_TOP = rand_pair_cor(cbinsize, te_pt, re_pt, top_e_neurons, 1000)
  RSC_BOT = rand_pair_cor(cbinsize, te_pt, re_pt, bot_e_neurons, 1000)
  RSC_INH = rand_pair_cor(cbinsize, ti_pt, ri_pt, I_Neurons, 500)
  #inputs
  ETI = zeros(length(P2))
  EBI = zeros(length(P1))
  III = zeros(Ni)
  for i = 1:length(P1)
    ETI[i] = mean(e_top_pt[i,:][:])
    EBI[i] = mean(e_bot_pt[i,:][:])
  end
  for i = 1:Ni
    III[i] = mean(i_pop_pt[i,:][:])
  end

  return wta_ness, bias, E_Rates, I_Rates, CV_TOP, CV_BOT, CV_INH, RSC_TOP, RSC_BOT, RSC_INH, ETI, EBI, III, FF_TOP, FF_BOT, FF_INH, te_pt, re_pt, e_top_pt, e_bot_pt, i_pop_pt
end
end

function WTAN_no_input(t, r, end_trans, Ne, Ni, half_e, min_e_neurons, min_i_neurons, fbinsize, cbinsize, FPe)
  P1s = round(Int,Ne/4-FPe)
  P1e = round(Int,Ne/4+FPe)
  P2s = round(Int,3*Ne/4-FPe)
  P2e = round(Int,3*Ne/4+FPe)

  P1 = P1s:P1e
  P2 = P2s:P2e

  rt = ((ntotal - end_trans)/1000.)*h

  ex = find(r .<= Ne)
  te = t[ex]
  re = r[ex]

  ix = find(r .> Ne)
  ti = t[ix]
  ri = r[ix]
tem = find(te.> end_trans)
  te_pt = te[tem]
  re_pt = re[tem]
  tim = find(ti.> end_trans)
  ti_pt = ti[tim]
  ri_pt = ri[tim]

  wta_ness, bias = score_analysis(re_pt, Ne)
  top_e_neurons, bot_e_neurons = Neurons_tb_ns(re_pt, half_e, 10, min_e_neurons)
  if ((top_e_neurons == -5) & (bot_e_neurons == -5))
    return "too few neurons spiking for analysis"
  else
  I_Neurons = Neuron_finder(ri_pt, 10, min_i_neurons)
  #rates
  E_Rates = [length(find(re_pt .== i))/rt for i=1:Ne]
  I_Rates = [length(find(ri_pt .== i))/rt for i=1:Ni]

  #counts
  E_count_top = count_train_intron(fbinsize, te_pt, re_pt, top_e_neurons, length(top_e_neurons), false)
  E_count_bot = count_train_intron(fbinsize, te_pt, re_pt, bot_e_neurons, length(bot_e_neurons), false)
  I_count_all = count_train_intron(fbinsize, ti_pt, ri_pt, I_Neurons, length(I_Neurons), false)
  #FF
  FF_TOP = fano_train(E_count_top, -5)
  FF_BOT = fano_train(E_count_bot, -5)
  FF_INH = fano_train(I_count_all, -5)
  #cv
  CV_TOP = CV_ISI_ALLTIME(top_e_neurons, te_pt, re_pt)
  CV_BOT = CV_ISI_ALLTIME(bot_e_neurons, te_pt, re_pt)
  CV_INH = CV_ISI_ALLTIME(I_Neurons, ti_pt, ri_pt)
  #synchrony
  RSC_TOP = rand_pair_cor(cbinsize, te_pt, re_pt, top_e_neurons, 1000)
  RSC_BOT = rand_pair_cor(cbinsize, te_pt, re_pt, bot_e_neurons, 1000)
  RSC_INH = rand_pair_cor(cbinsize, ti_pt, ri_pt, I_Neurons, 500)
  #inputs

  return wta_ness, bias, E_Rates, I_Rates, CV_TOP, CV_BOT, CV_INH, RSC_TOP, RSC_BOT, RSC_INH, FF_TOP, FF_BOT, FF_INH, te_pt, re_pt
end
end



function Rivalry_Analysis(t, r, Input, adapt, end_trans, Ne, Ni, half_e, min_e_neurons, min_i_neurons, fbinsize, cbinsize, netd_binsize, FPe)

  P1s = round(Int,Ne/4-FPe)
  P1e = round(Int,Ne/4+FPe)
  P2s = round(Int,3*Ne/4-FPe)
  P2e = round(Int,3*Ne/4+FPe)

  P1 = P1s:P1e
  P2 = P2s:P2e

  rt = ((ntotal - end_trans)/1000.)*h

  ex = find(r .<= Ne)
  te = t[ex]
  re = r[ex]

  ix = find(r .> Ne)
  ti = t[ix]
  ri = r[ix]

  aeb = zeros(length(adapt[1,:][:]))
  aet = zeros(length(aeb))
  for i = 1:half_e
    aeb .+= adapt[i,:][:]
    aet .+= adapt[i+half_e,:][:]
  end

  a_diff = aeb .- aet
  a_sum = aeb .+ aet
  a_frac = a_diff ./ a_sum
  a_frac = a_frac[2:end] #first entry is sometimes NaN
  a_flags, a_times = WLD_01(a_frac, -.333, .333)

  ntd, nts = nt_diff_H(te, re, ntotal, half_e, netd_binsize)
  s = ntd ./ nts #signal for dominances
  flags, times = WLD_01(s, -.333, .333)

  d = convert(Array{Float64}, diff(netd_binsize/(1000./h) .* times))
  cvd = cv(d) ###Raw estimate of CVD, likely to include very rapid switches which should really be smoothed out

  LP = .3

  dx = []
  for i in d
      if i > LP
          push!(dx, i)
      end
  end
  dx = convert(Array{Float64}, dx)
  cvdlp = cv(dx) ###Low-pass filter measure of CVD

  ###Use this code for spiking statistics, which includes only dominance and suppression times; mixed percept time is absorbed into either of those
  t2, f2 = splice_reversions(flags, times) ###Another way to get rid of rapid switches that aren't really there
  fw = find(f2 .== "win")
  fl = find(f2 .== "lose")
  tx = t2 .* netd_binsize
  d = convert(Array{Float64}, diff(tx))
  cvd2 = cv(d) ### Another estimate of CVD

  TN, BN = Neurons_tb_ns(re, half_e, 10, 100) #neurons in either pool who fired at least 10 spkes in simulation
  top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times, netd_binsize) #find win, lose, and draw times
  tbf, rbf = ligase(bot, bdom, te, re, BN) #bottom pool up states
  ttf, rtf = ligase(top, tdom, te, re, TN) #top pool up states
  tbdf, rbdf = ligase(top, tdom, te, re, BN) #bottom pool down states
  ttdf, rtdf = ligase(bot, bdom, te, re, TN) #top pool down states
  countFT = count_train_intron(fbinsize, ttf, rtf, TN, length(TN), false)
  countFB = count_train_intron(fbinsize, tbf, rbf, BN, length(BN), false)
  countFBD = count_train_intron(fbinsize, tbdf, rbdf, BN, length(BN), false)
  countFTD = count_train_intron(fbinsize, ttdf, rtdf, TN, length(TN), false)
  FF_TOP = fano_train(countFT, -5)
  FF_BOT = fano_train(countFB, -5)
  FF_TOPD = fano_train(countFTD, -5)
  FF_BOTD = fano_train(countFBD, -5)
  #correlations
  cwTu = rand_pair_cor(cbinsize, ttf, rtf, TN, 1000)
  cwBu = rand_pair_cor(cbinsize, tbf, rbf, BN, 1000)
  cwBd = rand_pair_cor(cbinsize, ttdf, rtdf, TN, 1000)
  cwTd = rand_pair_cor(cbinsize, tbdf, rbdf, BN, 1000)

  CV_TU = CV_ISI(top, TN, te, re)
  CV_BU = CV_ISI(bot, BN, te, re)
  CV_BD = CV_ISI(top, BN, tbdf, rbdf)
  CV_TD = CV_ISI(bot, TN, ttdf, rtdf)

  E_Input = Input[1:Ne, :]
  I_Input = Input[Ne+1:end, :]

  E_bot = E_Input[P1s:P1e, :]
  E_top = E_Input[P2s:P2e, :]

  etu = zeros(length(P2), Int64(tdom) + length(top));
  ebd = zeros(length(P1), Int64(tdom) + length(top));
  ebu = zeros(length(P2), Int64(bdom) + length(bot));
  etd = zeros(length(P1), Int64(bdom) + length(bot));

  println("Ligating Input Time Series for All Neurons\nThis part takes forever")
  println("If you aren't interested in measuring inputs, you can subsample in time or neurons, or not measure inputs at all")
  println("To do that, you'll have to mess with Euler_W.jl")
  println("Alternatively, you can comment this part out")
  for i = 1:length(P1)
    println("On the $(i)th iteration")
    etu[i,1:end], ebd[i,1:end] = ligate_inputs(top, Input[P2[i],:][:], Input[P1[i],:][:])
    ebu[i,1:end], etd[i,1:end] = ligate_inputs(bot, Input[P1[i],:][:], Input[P2[i],:][:])
  end

  etwm = zeros(length(P2))
  ebwm = zeros(length(P2))
  etlm = zeros(length(P2))
  eblm = zeros(length(P2))
  for i in eachindex(P2)
    etwm[i] = mean(etu[i,:])
    etlm[i] = mean(etd[i,:])
    ebwm[i] = mean(ebu[i,:])
    eblm[i] = mean(ebd[i,:])
  end

  return te_pt, re_pt, TN, BN, d, cvd, flags, times, cvdlp, t2, f2, cvd2, top, tdom, bot, bdom, nmz, tnmz, FF_TOP, FF_BOT, FF_TOPD, FF_BOTD, cwTu, cwTd, cwBu, cwBd, CV_TU, CV_BU, CV_TD, CV_BD, etwm, ebwm, etlm, eblm
end

function write_array(filename, a)
  newfile = open(filename, "w")
  for i in eachindex(a)
    text = "$(a[i])\n"
    write(newfile, text)
  end
  close(newfile)
end
#
function write_raster(filename, t, r)
  newfile = open(filename, "w")
  for i in eachindex(t)
    text = "$(t[i]), $(r[i])\n"
    write(newfile, text)
  end
  close(newfile)
end

function ligate_inputs(introns, input1, input2)
  i1, i2 = Float64[0], Float64[0]
  for i in introns
    s1 = input1[i[1]:i[2]]
    s2 = input2[i[1]:i[2]]
    i1 = cat(1, i1, s1)
    i2 = cat(1, i2, s2)
  end
  s1, s2 = i1[2:end], i2[2:end]
  return s1, s2
end

function zscore(x)
    m = mean(x)
    s = std(x)
    x1 = x .- m
    x2 = x1 ./s
    return x2
end


# for i = 1:runtime
#   EXW[i] = sum(Exc[P1s:P1e, i][:])
#   IXW[i] = sum(Inh[P1s:P1e, i][:])
#   EXL[i] = sum(Exc[P2s:P2e, i][:])
#   IXL[i] = sum(Inh[P2s:P2e, i][:])
# end


# for i in range(N):
#   for j in range(N):
#     for k in range(N):
#       text = "julia ./run_simplest.jl {0} {1} {2}\n".format(Aee[i], Aie[j], Aie_NL[k])
#       newfile.write(text)
function sim_2_theory(SEE, SEI, SIE, SIEL, SII, fe, fi, cth, re1, re2, ri1, ri2, n)

    FE = fe - cth
    FI = fi - cth

    n0 = div(n, 2)

    see1 = zeros(n0)
    see2 = zeros(n0)
    sie1 = zeros(n0)
    sie2 = zeros(n0)
    sieL1 = zeros(n0)
    sieL2 = zeros(n0)
    sii1 = zeros(n0)
    sii2 = zeros(n0)
    sei1 = zeros(n0)
    sei2 = zeros(n0)

    for i = 1:n0
        see1[i] = mean(SEE[i,:][:])
        see2[i] = mean(SEE[i+n0,:][:])
        sie1[i] = mean(SIE[i,:][:])
        sie2[i] = mean(SIE[i+n0,:][:])
        sieL1[i] = mean(SIEL[i,:][:])
        sieL2[i] = mean(SIEL[i+n0,:][:])
        sei1[i] = mean(SEI[i,:][:])
        sei2[i] = mean(SEI[i+n0,:][:])
        sii1[i] = mean(SII[i,:][:])
        sii2[i] = mean(SII[i+n0,:][:])
    end

    WEE_1 = mean(see1)/re1
    WEE_2 = mean(see2)/re2
    WIE_1 = mean(sie1)/re1
    WIE_2 = mean(sie2)/re2
    WIEL_1 = mean(sieL1)/re2
    WIEL_2 = mean(sieL2)/re1
    WEI_1 = mean(sei1)/ri1
    WEI_2 = mean(sei2)/ri2
    WII_1 = mean(sii1)/ri1
    WII_2 = mean(sii2)/ri2

    return WEE_1, WEE_2, WIE_1, WIE_2, WIEL_1, WIEL_2, WEI_1, WEI_2, WII_1, WII_2, FE, FI
end

function estimate_I(SEE, SEI, SIE, SIEL, SII, s_strength, n)

  n0 = div(n, 2)
  sE1 = zeros(n0)
  sE2 = zeros(n0)
  sI1 = zeros(n0)
  sI2 = zeros(n0)

  for i = 1:n0
      sE1[i] = mean(SEE[i,:][:] .+ SEI[i,:][:] .+ s_strength)
      sE2[i] = mean(SEE[i+n0,:][:] .+ SEI[i+n0,:][:] .+ s_strength)
      sI1[i] = mean(SIE[i,:][:] .+ SIEL[i,:][:] .+ SII[i,:][:])
      sI2[i] = mean(SIE[i+n0,:][:] .+ SIEL[i+n0,:][:] .+ SII[i+n0,:][:])
  end

  sE1m = mean(sE1)
  sE2m = mean(sE2)
  sI1m = mean(sI1)
  sI2m = mean(sI2)

  return sE1m, sE2m, sI1m, sI2m
end

function estimated_gtiX(I)
  return (I-.6)*.05
end

function estimated_gti(I)
  return I*.05
end

function theory_rates(WEE_1, WEE_2, WIE_1, WIE_2, WIEL_1, WIEL_2, WEI_1, WEI_2, WII_1, WII_2, FE, FI)

  RE1 = ((WII_1*FE) - (WEI_1*FI))/(((WEI_1*(WIE_1+WIEL_1))-(WII_1*WEE_1)))
  RI1 = ((FE*(WIE_1+WIEL_1))-(WEE_1*FI))/(((WEI_1*(WIE_1+WIEL_1))-(WII_1*WEE_1)))

  RE2 = ((WII_2*FE) - (WEI_2*FI))/(((WEI_2*(WIE_2+WIEL_2))-(WII_2*WEE_2)))
  RI2 = ((FE*(WIE_2+WIEL_2))-(WEE_2*FI))/(((WEI_2*(WIE_2+WIEL_2))-(WII_2*WEE_2)))

  return RE1, RI1, RE2, RI2
end

function sim_2_theory_2x2(SEE, SEI, SIE, SII, fe, fi, cth, re1, ri1, n)

    FE = fe - cth
    FI = fi - cth

    see1 = zeros(n)
    sie1 = zeros(n)
    sii1 = zeros(n)
    sei1 = zeros(n)

    for i = 1:n
        see1[i] = mean(SEE[i,:][:])
        sie1[i] = mean(SIE[i,:][:])
        sei1[i] = mean(SEI[i,:][:])
        sii1[i] = mean(SII[i,:][:])
    end

    WEE_1 = mean(see1)/re1
    WIE_1 = mean(sie1)/re1
    WEI_1 = mean(sei1)/ri1
    WII_1 = mean(sii1)/ri1

    return WEE_1, WIE_1, WEI_1, WII_1, FE, FI
end

function theory_rates_2x2(WEE_1, WIE_1, WEI_1, WII_1, FE, FI)

  RE1 = ((WII_1*FE) - (WEI_1*FI))/((WEI_1*WIE_1)-(WII_1*WEE_1))
  RI1 = ((WIE_1*FE) - (WEE_1*FI))/((WEI_1*WIE_1)-(WII_1*WEE_1))

  return RE1, RI1
end

function theory_pos_eigV(b, c, d)
  return b + sqrt(c+d)
end

#renormalized gain function values, which are W coefficients, set to 1 for now

function get_b(wee, wii)
  p = (wee .- wii .- 2)
  return p./2.
end

function get_c(wee, wei, wie, wii)
  p1 = (wii .+ wee) .^ 2.
  p2 = wei .* wie .* 4.
  return p1 .- p2
end

function get_d(wei, wieL)
  return 4. .* wei .* wieL
end

#renormalized gain function values, which are W coefficients, are free parameters

function get_b_g(wee, wii, ge, gi)
  p = ((ge*wee) .- (gi*wii) .- 2)
  return p./2.
end

function get_c_g(wee, wei, wie, wii, ge, gi)
  p1 = ((gi*wii) .+ (ge*wee)) .^ 2.
  p2 = wei .* wie .* 4. .* ge .* gi
  return p1 .- p2
end

function get_d_g(wei, wieL, ge, gi)
  return 4. .* wei .* wieL .* ge .* gi
end

function FI_G(S, vth, tau)
  denom = tau*log((S-(vth/tau))/S)
  return -1/denom
end

function firing_period(tau, vth, S)
  return -tau*log(1-(vth/(tau*S)))
end

function single_neuron(tau, vth, S, runtime, h)
  v0 = randn()*vth
  ntotal = Int64(div(runtime, h))
  v = v0
  m_leak = h/tau
  spikes = []
  for i = 1:ntotal
    v += (S*h) - (v * m_leak)
    if v > vth
      push!(spikes, i)
      v -= vth
    end
  end
  return spikes
end

function single_neuron_noise(tau, vth, S, runtime, h)
  v0 = randn()*vth
  ntotal = Int64(div(runtime, h))
  v = v0
  m_leak = h/tau
  spikes = []
  sh = sqrt(h)
  sigma = sqrt(S)
  ssh = sh*sigma
  eta = (randn(ntotal)*ssh) + S*h
  for i = 1:ntotal
    v = v - (v*m_leak) + eta[i]
    if v > vth
      push!(spikes, i)
      v -= vth
    end
  end
  return spikes
end
