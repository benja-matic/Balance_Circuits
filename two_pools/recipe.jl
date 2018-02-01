include("simplest_analysis.jl")
# kee = parse(Float64, ARGS[1])
# kei = parse(Float64, ARGS[2])
# kie_L = parse(Float64, ARGS[3])
# kii = parse(Float64, ARGS[4])
#
# Aee = parse(Float64, ARGS[5])
# Aei = parse(Float64, ARGS[6])
# Aie_L = parse(Float64, ARGS[7])
# Aii= parse(Float64, ARGS[8])
# s_strength = parse(Float64, ARGS[9])
# p = parse(Float64, ARGS[10])
# Ne = parse(Float64, ATGS[11])

function sd_2_k(sd)
  return 1./(2*pi*((sd)^2))
end

function recipe_part_1(c_max)

  c = 0

  s_strength = 2.
  p = .2

  h = .1
  Ne = 1600
  Ne = Int64(Ne)
  Ni = div(Ne, 4)

  min_e_neurons = 200
  min_i_neurons = 50
  half_e = div(Ne, 2)
  half_i = div(Ni, 2)

  runtime = 20000 #ms
  ntotal = round(runtime/h) #time points
  fbinsize = 400/h
  cbinsize = 100/h
  vth = 20.

  end_trans = 1000
  rt = ((ntotal - end_trans)/1000.)*h
  success = false


  while c < c_max
    kee=sd_2_k(rand(15:100)/100.)
    kei=sd_2_k(rand(20:100)/100.)
    kie_L=sd_2_k(rand(20:100)/100.)
    kii=sd_2_k(rand(15:100)/100.)

    Aee=rand(1:100)
    Aei=rand(1:1000)
    Aie_L=rand(1:1000)
    Aii=rand(1:1000)

    Aee /= p
    Aei /= p
    Aie_L /= p
    Aii /= p

    wee,wei,wie,wii = weights(Ne,Ni,kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, p,p,p,p);
    te,re,ti,ri,kill_flag = lif(h,runtime,Ne,Ni,wee,wei,wie,wii,s_strength,s_strength);

    tem = find(te.> end_trans)
    te_pt = te[tem]
    re_pt = re[tem]
    tim = find(ti.> end_trans)
    ti_pt = ti[tim]
    ri_pt = ri[tim]

    wta_ness, bias = score_analysis(re_pt, Ne)
    top_e_neurons, bot_e_neurons = Neurons_tb_ns(re_pt, half_e, 10, min_e_neurons)#22

    if ((top_e_neurons == -5) & (bot_e_neurons == -5))
      println("##RESULT -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
      #Set stuff equal to error code
      #println("##RESULT stuff")
    else
      I_Neurons = Neuron_finder(ri_pt, 10, min_i_neurons)
      #rates
      E_rate = (length(re_pt)/rt)/Ne
      I_rate = (length(ri_pt)/rt)/Ni

      #counts
      E_count_top = count_train(fbinsize, te_pt, re_pt, top_e_neurons, length(top_e_neurons))
      E_count_bot = count_train(fbinsize, te_pt, re_pt, bot_e_neurons, length(bot_e_neurons))
      I_count_all = count_train(fbinsize, ti_pt, ri_pt, I_Neurons, length(I_Neurons))
      #FF
      ntf = network_fano(te_pt, re_pt, cbinsize, Int64(ntotal))
      E_FANO_mean_top, E_FANO_median_top, E_FANO_std_top = fano_train(E_count_top, -5)
      E_FANO_mean_bot, E_FANO_median_bot, E_FANO_std_bot = fano_train(E_count_bot, -5)
      I_FANO_mean, I_FANO_median, I_FANO_std = fano_train(I_count_all, -5)
      #cv
      E_CV_mean_top, E_CV_median_top, E_CV_STD_top = CV_ISI_ALLTIME(top_e_neurons, te_pt, re_pt)
      E_CV_mean_bot, E_CV_median_bot, E_CV_STD_bot = CV_ISI_ALLTIME(bot_e_neurons, te_pt, re_pt)
      I_CV_mean, I_CV_median, I_CV_STD = CV_ISI_ALLTIME(I_Neurons, ti_pt, ri_pt)
      #synchrony
      E_spike_correlation_top = rand_pair_cor(cbinsize, te_pt, re_pt, top_e_neurons, 1000)
      E_spike_correlation_bot = rand_pair_cor(cbinsize, te_pt, re_pt, bot_e_neurons, 1000)
      I_spike_correlation = rand_pair_cor(cbinsize, ti_pt, ri_pt, I_Neurons, 500)

      c+=1
      if (2. > E_FANO_mean_top > .8) & (2. > E_FANO_mean_bot > .8) & (2. > I_FANO_mean > .8) & (E_spike_correlation_top < .1) & (E_spike_correlation_bot < .1) & (I_spike_correlation < .1)
        success = true
        println("##RESULT success_1, $(c), $(E_rate), $(I_rate), $(wta_ness), $(bias), $(ntf), $(E_FANO_mean_top), $(E_FANO_mean_bot), $(I_FANO_mean), $(E_CV_mean_top), $(E_CV_mean_bot), $(I_CV_mean), $(E_spike_correlation_top), $(E_spike_correlation_bot), $(I_spike_correlation), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
        return c, kee, kei, kie_L, kii, Aee*p, Aei*p, Aie_L*p, Aii*p, wta_ness
      else
        println("##RESULT $(success), $(c), $(E_rate), $(I_rate), $(wta_ness), $(bias), $(ntf), $(E_FANO_mean_top), $(E_FANO_mean_bot), $(I_FANO_mean), $(E_CV_mean_top), $(E_CV_mean_bot), $(I_CV_mean), $(E_spike_correlation_top), $(E_spike_correlation_bot), $(I_spike_correlation), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
      end
    end
    if c >= c_max -1
      return c, -5, -5, -5, -5, -5, -5, -5, -5, wta_ness
    end
  end
end

function recipe_part_2(kee, kei, kie_L, kiim, Aee, Aei, Aie_L, Aii, c_max)

  c = 0
  s_strength = 2.
  p = .2

  h = .1
  Ne = 1600
  Ne = Int64(Ne)
  Ni = div(Ne, 4)

  min_e_neurons = 200
  min_i_neurons = 50
  half_e = div(Ne, 2)
  half_i = div(Ni, 2)

  runtime = 20000 #ms
  ntotal = round(runtime/h) #time points
  fbinsize = 400/h
  cbinsize = 100/h
  vth = 20.

  end_trans = 10000
  rt = ((ntotal - end_trans)/1000.)*h
  success = false

  Aei /= p
  Aie_L /= p

  while c < c_max

    Aee = rand(50:200)
    Aii = rand(1:1000)
    kee=sd_2_k(rand(15:40)/100.)
    kii=sd_2_k(rand(15:40)/100.)

    Aee /= p
    Aii /= p

    #println("preparing to simulate")

    wee,wei,wie,wii = weights(Ne,Ni,kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, p,p,p,p)
    te,re,ti,ri,kill_flag = lif(h,runtime,Ne,Ni,wee,wei,wie,wii,s_strength,s_strength)

    #println("simulation complete, preparing to analyze")

    tem = find(te.> end_trans)
    te_pt = te[tem]
    re_pt = re[tem]
    tim = find(ti.> end_trans)
    ti_pt = ti[tim]
    ri_pt = ri[tim]

    wta_ness, bias = score_analysis(re_pt, Ne)
    top_e_neurons, bot_e_neurons = Neurons_tb_ns(re_pt, half_e, 10, min_e_neurons)#22

    if ((top_e_neurons == -5) & (bot_e_neurons == -5))
      if kill_flag == true
        println("killed due to super fast spiking")
      end
      println("##RESULT -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
      #Set stuff equal to error code
      #println("##RESULT stuff")
    else
      I_Neurons = Neuron_finder(ri_pt, 10, min_i_neurons)
      #rates
      E_rate = (length(re_pt)/rt)/Ne
      I_rate = (length(ri_pt)/rt)/Ni

      #counts
      E_count_top = count_train(fbinsize, te_pt, re_pt, top_e_neurons, length(top_e_neurons))
      E_count_bot = count_train(fbinsize, te_pt, re_pt, bot_e_neurons, length(bot_e_neurons))
      I_count_all = count_train(fbinsize, ti_pt, ri_pt, I_Neurons, length(I_Neurons))
      #FF
      ntf = network_fano(te_pt, re_pt, cbinsize, Int64(ntotal))
      E_FANO_mean_top, E_FANO_median_top, E_FANO_std_top = fano_train(E_count_top, -5)
      E_FANO_mean_bot, E_FANO_median_bot, E_FANO_std_bot = fano_train(E_count_bot, -5)
      I_FANO_mean, I_FANO_median, I_FANO_std = fano_train(I_count_all, -5)
      #cv
      E_CV_mean_top, E_CV_median_top, E_CV_STD_top = CV_ISI_ALLTIME(top_e_neurons, te_pt, re_pt)
      E_CV_mean_bot, E_CV_median_bot, E_CV_STD_bot = CV_ISI_ALLTIME(bot_e_neurons, te_pt, re_pt)
      I_CV_mean, I_CV_median, I_CV_STD = CV_ISI_ALLTIME(I_Neurons, ti_pt, ri_pt)
      #synchrony
      E_spike_correlation_top = rand_pair_cor(cbinsize, te_pt, re_pt, top_e_neurons, 1000)
      E_spike_correlation_bot = rand_pair_cor(cbinsize, te_pt, re_pt, bot_e_neurons, 1000)
      I_spike_correlation = rand_pair_cor(cbinsize, ti_pt, ri_pt, I_Neurons, 500)

      c+=1
      println("##RESULT $(success), $(c), $(E_rate), $(I_rate), $(wta_ness), $(bias), $(ntf), $(E_FANO_mean_top), $(E_FANO_mean_bot), $(I_FANO_mean), $(E_CV_mean_top), $(E_CV_mean_bot), $(I_CV_mean), $(E_spike_correlation_top), $(E_spike_correlation_bot), $(I_spike_correlation), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
      if bias == 1
        if (2. > E_FANO_mean_top > .8) & (2. > I_FANO_mean > .8) & (E_spike_correlation_top < .1) & (wta_ness > .8)
          success = true
          #plot(te, re, "g.")
          println("##RESULT success_2, $(c), $(E_rate), $(I_rate), $(wta_ness), $(bias), $(ntf), $(E_FANO_mean_top), $(E_FANO_mean_bot), $(I_FANO_mean), $(E_CV_mean_top), $(E_CV_mean_bot), $(I_CV_mean), $(E_spike_correlation_top), $(E_spike_correlation_bot), $(I_spike_correlation), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
          return c, kee, kei, kie_L, kii, Aee*p, Aei*p, Aie_L*p, Aii*p, wta_ness
        end
      else bias == 2
          if (2. > E_FANO_mean_bot > .8) & (2. > I_FANO_mean > .8) & (E_spike_correlation_bot < .1) & (wta_ness > .8)
            success = true
            plot(te, re, "g.")
            println("##RESULT success_2, $(c), $(E_rate), $(I_rate), $(wta_ness), $(bias), $(ntf), $(E_FANO_mean_top), $(E_FANO_mean_bot), $(I_FANO_mean), $(E_CV_mean_top), $(E_CV_mean_bot), $(I_CV_mean), $(E_spike_correlation_top), $(E_spike_correlation_bot), $(I_spike_correlation), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(p), $(Ne)")
            return c, kee, kei, kie_L, kii, Aee*p, Aei*p, Aie_L*p, Aii*p, wta_ness
          end
      end
    end
  end
end

function recipe_part_3(kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, c_max)

  c = 0
  s_strength = 2.
  p = .2

  h = .1
  Ne = 1600
  Ne = Int64(Ne)
  Ni = div(Ne, 4)

  min_e_neurons = 200
  min_i_neurons = 50
  half_e = div(Ne, 2)
  half = half_e
  half_i = div(Ni, 2)

  runtime = 20000 #ms
  ntotal = round(runtime/h) #time points
  fbinsize = 400/h
  cbinsize = 100/h
  vth = 20.

  end_trans = 1000
  rt = ((ntotal - end_trans)/1000.)*h
  success = false
  tauadapt = 1000.

  Aee /= p
  Aei /= p
  Aie_L /= p
  Aii /= p
  srand(12345678)
  g = [.00001*(1.09^i) for i=1:100]
  for i in g
    c+=1

    #gamma = exp(-rand(299000000:1381000000)/100000000.)
    gamma = i

    wee,wei,wie,wii = weights(Ne,Ni,kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, p,p,p,p);
    #te,re,ti,ri = euler_lif(h,runtime,Ne,Ni,wee,wei,wie,wii,s_strength,s_strength,gamma,tauadapt);
    te,re,ti,ri,kill_flag, e_top, e_bot, aet, aeb, ev = lif(h,runtime,Ne,Ni,wee,wei,wie,wii,s_strength,s_strength,i,tauadapt);


    ntd, nts = nt_diff(te, re, ntotal, half, 25/h)
    s = ntd./nts #signal for dominances
    flags, times = WLD_01(s)
    top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times)

    if (tdom > 20000) & (bdom > 20000)
      #rivalry
      println("RIVALING")
      TN, BN = Neurons_tb_ns(re, half, 10, 100)
      if ((TN == -5) | (BN == -5))
        println("##RESULT -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, $(success), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(c), $(gamma)")
        cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10
      else
      t2, f2 = splice_reversions(flags, times)
      MD = mean(diff(250*t2))
      MDT, MDB = (tdom*1./length(bot)), (bdom*1./length(top))
      tbf, rbf = ligase(bot, bdom, te, re, BN)
      ttf, rtf = ligase(top, tdom, te, re, TN)
      tbdf, rbdf = ligase(top, tdom, te, re, BN)
      ttdf, rtdf = ligase(bot, bdom, te, re, TN)
      #countCT = count_train_intron(cbinsize, ttf, rtf, TN, 50, true)
      #countCB = count_train_intron(cbinsize, tbf, rbf, BN, 50, true)
      countFT = count_train_intron(fbinsize, ttf, rtf, TN, length(TN), false)
      countFB = count_train_intron(fbinsize, tbf, rbf, BN, length(BN), false)
      #cwT, cwB = correlate_within(countCT, -5), correlate_within(countCB, -5)
      cwT = rand_pair_cor(cbinsize, ttf, rtf, TN, 1000)
      cwB= rand_pair_cor(cbinsize, tbf, rbf, BN, 1000)
      MFT, medFT, stdFT = fano_train(countFT, -5)
      MFB, medFB, stdFB = fano_train(countFB, -5)
      mcvt, medcvt, stdcvt = CV_ISI(top, TN, te, re)
      mcvb, medcvb, stdcvb = CV_ISI(bot, BN, te, re)
      if (MFT > .8) & (MFT < 2.5) & (MFB > .8) & (MFB < 2.5) & (cwT < .1) & (cwB < .1) & (MD > 5000)
        println("party time")
        success = true
        println("##RESULT 1, $(MDT), $(MDB), $(MD), $(cwT), $(cwB), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(success), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(c), $(gamma)")
        return c, kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, gamma, tauadapt, MDT
      end
      println("##RESULT 1, $(MDT), $(MDB), $(MD), $(cwT), $(cwB), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(success), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(c), $(gamma)")
    end

    elseif tnmz > 160000
      println("NORMALIZATION")
      #normalization
      TN, BN = Neurons_tb_ns(re, half, 10, 200)
      if ((TN == -5) | (BN == -5))
        println("##RESULT -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, $(success), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(c), $(gamma)")
        cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8
      else
      tbf, rbf = ligase(nmz, tnmz, te, re, BN)
      ttf, rtf = ligase(nmz, tnmz, te, re, TN)
      #countCT = count_train_intron(cbinsize, ttf, rtf, TN, 50, true)
      #countCB = count_train_intron(cbinsize, tbf, rbf, BN, 50, true)
      countFT = count_train_intron(fbinsize, ttf, rtf, TN, length(TN), false)
      countFB = count_train_intron(fbinsize, tbf, rbf, BN, length(BN), false)
      cwT = rand_pair_cor(cbinsize, ttf, rtf, TN, 1000)
      cwB = rand_pair_cor(cbinsize, tbf, rbf, BN, 1000)
      #cwT, cwB = correlate_within(countCT, -5), correlate_within(countCB, -5)
      MFT, medFT, stdFT = fano_train(countFT, -5)
      MFB, medFB, stdFB = fano_train(countFB, -5)
      mcvt, medcvt, stdcvt = CV_ISI(nmz, TN, te, re)
      mcvb, medcvb, stdcvb = CV_ISI(nmz, BN, te, re)
      println("##RESULT 0, -7, -7, $(tnmz), $(cwT), $(cwB), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(success), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(c), $(gamma)")
    end

      #mtd, vtd, mbu, vbu = input_analysis(bot, e_top, e_bot)
    elseif (tdom > 180000) | (bdom > 180000)
      println("WINNER TAKE ALL")
      TN, BN = Neurons_tb_ns(re, half, 10, 200)
      n1, n2 = length(TN), length(BN)
      println("here's $(n1), and $(n2)")
      if ((TN == -5) & (BN == -5))
        println("##RESULT -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, $(success), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(c), $(gamma)")
        cwT, cwB, MFT, MFB, mcvt, mcvb, vtu, mtu, vtd, mtd, vbu, mbu, vbd, mbd = -8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8
      else
      if n1 > n2
        Neurons = TN
        s1 = e_top
        s2 = e_bot
        exon = top
        dom = tdom
      else
        Neurons = BN
        s1 = e_bot
        s2 = e_top
        exon = bot
        dom = bdom
      end
      tbf, rbf = ligase(exon, dom, te, re, Neurons)
      #countCB = count_train_intron(cbinsize, tbf, rbf, Neurons, 10, true)
      countFB = count_train_intron(fbinsize, tbf, rbf, Neurons, length(Neurons), false)
      cwB = rand_pair_cor(cbinsize, tbf, rbf, Neurons, 1000)
      #cwB = correlate_within(countCB, -5)
      MFB, medFB, stdFB = fano_train(countFB, -5)
      mcvb, medcvb, stdcvb = CV_ISI(exon, Neurons, te, re)
      mcvt, medcvt, stdcvt, MFT, medFT, stdFT = -7, -7, -7, -7, -7, -7
      cwT = -7
      println("##RESULT 2, -7, -7, $(tnmz), -7, $(cwB), $(MFT), $(medFT), $(stdFT), $(MFB), $(medFB), $(stdFB), $(mcvt), $(medcvt), $(stdcvt), $(mcvb), $(medcvb), $(stdcvb), $(success), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(c), $(gamma)")
    end
    else
      #some kind of error
      println("##RESULT -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9,-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, $(success), $(kee), $(kei), $(kie_L), $(kii), $(Aee*p), $(Aei*p), $(Aie_L*p), $(Aii*p), $(s_strength), $(c), $(gamma)")
    # end
    end #end if statement

  end # end while loop
end #end function






c_irg, kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, wta_ness = recipe_part_1(1000)

if kee == -5
  exit()
else
println("Step 1 Complete")

c_wta, kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, wta_ness = recipe_part_2(kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, 1000)
end
if kee == -5
  exit()
else
include("r5.jl")

c_riv, kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, gamma, tauadapt, MDT = recipe_part_3(kee, kei, kie_L, kii, Aee, Aei, Aie_L, Aii, 1000)

end
end



#defualt paramters
# s_strength = 2.
# p = .2
# Ne = 1600
# Ni = 400
#
# Aee = 100
# Aei = 500
# Aie_L = 500
# Aii = 500
#
# kee = sd_2_k(.25)
# kei = sd_2_k(.8)
# kie_L = sd_2_k(.8)
# kii = sd_2_k(.5)
#
Aee /= p
Aei /= p
Aie_L /= p
Aii /= p
#
# h = .1
# runtime = 20000.
#
#
# function tally(zd)
#     ret = zeros(Int64, K)
#     for k in zd
#         ret[k] += 1
#     end
#     return ret
# end




#
