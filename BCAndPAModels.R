
# Guillaume Deffuant 2025
# Models of information consultation
# See paper entitled: "How opinions get more extreme in an era of abundant information"
# Authors: Guillaume Deffuant, Marijn Keijzer, Sven Banisch

#================
# examples of BC model runs
#=================

# running one trajectory and visu
# tr = trajBC(150, nContent = 100, visu = TRUE)
# running 10000 trajectories of 150 steps and visu distribution
# sm = distribBC(steps = 150, nTraj = 10000, nContent = 100)


#================
# examples of PA model runs
#=================

# running one trajectory and voisu
# tr = trajPAT(150, N = 100, beta = 0.15, visu = TRUE)
# running 10000 trajectoies of 150 steps and visu
# sm = PATDistribRun(steps = 150, nTraj = 10000, N = 100, beta = 0.15, visu = TRUE)

# running Markov model
# sm = runPATDens(steps = 150, N = 100, beta = 0.15)


#================================================
# Bounded Confidence model (agent version only)
#================================================

#-----------------
# Function of shrinking confidence bound with extremity of attitude (x)
boundc = function(x, bmax, beta, gamma = 2){
  
  #gamma = 2
  return(bmax * exp(-beta * abs(x)^gamma))
}


#------------------
#One BC step taking content in the whole set of contents
# returns new value of x
stepBC = function(x, content, bmax, beta, mu, gamma){
  
  b = boundc(x, bmax, beta, gamma)
  y = sample(content, 1)
  if (abs(x - y) <= b) x = x + mu *(y-x)  
  return(x)
}

#------------------
#One BC step taking content in agent's confidence segment
# returns new value of x
stepBCIn = function(x, b, content, mu){
  
  sp = content[content >= (x - b) & content <= (x + b)]
  if (length(sp) > 0){ 
    x = x + mu *(sample(sp, 1)-x)  
  }
  return(x)
}


#------------------
#One BC step taking content in the whole set of contents
# returns new value of x
stepBC = function(x, b, content, mu){
  
  c = sample(content, 1)
  if (abs(x - c) <= b) x = x + mu *(c-x)  
  return(x)
}

#---------------------
# if sig = 0, uniform proba of content, otherwise Gaussian of mean sig, 
# Content a set of values chosen from this proba with a number = nContent
# single trajectory
trajBC = function(steps, sig = 0, nContent, bmax = 0.4, beta = 2, mu = 0.5, prejud = 0, gamma = 2, visu = FALSE){
  
  
  if (sig == 0) content = runif(nContent, -1, 1)
  else {
    content = rnorm(nContent, 0, sig)
    content = content[abs(content) <= 1]
  }
  if (visu){
    ksteps = 5
    dec = 4
    pchIn = 19
    pchOut = 4
    par(mar = c(4, 4, 0.1, 0.1), cex = 1.2)
    plot(-10,-1, xlim = c(-1,1), ylim = c(-dec, steps), xlab = "Attitude", ylab = "Item consultations", xaxt = "n", yaxt = "n", axes = FALSE)
    axis(1, at = c(-1, 0, 1), col = "white", col.ticks = "black")
    axis(2, at = (1:ksteps)* trunc(steps / ksteps), col = "white", col.ticks = "black")
    for(i in 1:nContent) lines(c(content[i], content[i]), c(-dec,0), col = colExtr((abs(content[i]))), lwd = 2)
  }
  tr = rep(0, steps)
  x = prejud
  tr[1] = x
  for (j in (1:steps)){
    b = boundc(x, bmax, beta, gamma)
    if (visu) lines(c(x-b, x+b), c(j, j), col = "lightgrey", lwd = 2)
    inSeg = FALSE
    c = sample(content, 1)
    if (abs(x - c) <= b){ 
      x = x + mu *(c-x)
      inSeg = TRUE
    }
    tr[j+1] = x
    if (visu){
      points(c, j, col = colExtr(abs(c)), pch = ifelse(inSeg, pchIn, pchOut))
      lines(c(tr[j], x), c(j, j+1), col = "black", lwd = 2)
    }
  }
  return(tr)
}

#-------------
# to compute histogram in PATDens()
inLims = function(x, lims){
  
  for(i in 2:length(lims)){
    if (x <= (lims[i]+lims[i-1])/2) return(i-1)
  }
  return(length(lims))
}


#---------------------
# if sig = 0, uniform proba of content, otherwise Gaussian of mean sig, 
# Content a set of values chosen from this proba with a number = nContent
# watches is the number of contents watched in one day
# nTraj trajectories
# Returns matrix of densities
distribBC = function(steps, sig = 0, nTraj, nContent, bmax = 0.4, beta = 2, mu = 0.5, prejud = 0, ng = 50, gamma = 2, inSeg = FALSE, visu = TRUE){
  
  ngt = 2*ng+1
  svm = matrix(rep(0, ngt*(steps)), nrow = ngt, ncol = steps)
  colnames(svm) = 1:steps
  lims = -1 + (0:(2*ng))/ng
  rownames(svm) = -1 + (0:(2*ng))/ng
  bm = bmax * exp(- beta)
  for(i in 1:nTraj){
    cat("\r   ",i)
    if (sig == 0) content = runif(nContent, -1, 1)
    # (sig == 0) content = runif(nContent, -1-bm, 1+bm)
    else {
      content = rnorm(nContent, 0, sig)
      content = content[abs(content) <= 1]
    }
    x = prejud
    for (j in (1:steps)){
      b = boundc(x, bmax, beta, gamma)
      if (inSeg) x = stepBCIn(x, b, content, mu)
      else x = stepBC(x, b, content, mu)
      k = inLims(x, lims)
      svm[k, j] = svm[k, j]+1
    }
  }
  svm = svm/nTraj
  if (visu) visuMatKSteps(sv = svm, ksteps = 5, lwd = 6)
  return(svm)
}

#================================================
# Argument model in agent version
#================================================

#---------------
# prejud is the prejudice
# beta is a parameter equal to beta*alpha in the paper. 
# if beta < 0, beta should be in (0, 0.5). proba = 0.5-beta to believe congruent item (0.5+beta) to believe incongruent items)
# if beta >= 0  proba = 1 / (1 + exp(- beta * (b+prejud)) to believe positive items
bel = function(b, prejud, beta){
  
  if (beta >= 0) return(1 / (1 + exp(- beta * (b + prejud))))
  if (b + prejud > 0) return(0.5- beta)
  if (b+prejud == 0) return(0.5)
  return(0.5+beta)
}


#-------------------------------------------
# One step PAT with one arg chosen at random
# changes beliefs that can switch from 0 to 1 or from 1 to 0
stepPAT = function(beliefsp, beliefsn, prejud, beta){
  
  b = sum(beliefsp) - sum(beliefsn)
  pp = bel(b, prejud, beta)
  nta = length(beliefsp)
  i = sample(1:(2*nta), 1)
  if (i <= nta){
    if (runif(1) < pp) beliefsp[i] = 1 else beliefsp[i] = 0
    bel = beliefsp[i]
    }
  else {
    i = i - nta
    if (runif(1) < (1-pp)) beliefsn[i] = 1 else beliefsn[i] = 0
    bel = beliefsn[i]
    i = -i
  }
  return(list('bp'= beliefsp, 'bn' = beliefsn, 'arg' = i, 'bel' = bel))
}


#-----------------------------------------------------
# Returns one trajectory of attitudes of steps iterations
# qteps = nb steps of the trajectory
# N = total number of argument N should be even N = 2*nta
trajPAT = function(steps, N, beta, prejud = 0, visu = FALSE){
  
  aa = exp(- beta * (prejud))
  # vv = exp(-  (prejud + sum(bp) - sum(bn)))
  aa = (1-aa)/(1+aa) 
  traja = aa
  trajb = prejud
  if (N %% 2 != 0) N = N+1 # N should be even
  nta = N / 2
  bp = rep(0, nta)
  bn = rep(0, nta)
  if (visu){
    ksteps = 5
    dec = 4
    pch = c(4, 19)
    cl = c("white", "green")
    par(mar = c(4, 4, 0.1, 0.1), cex = 1.2)
    plot(-10,-1, xlim = c(-1,1), ylim = c(0, steps), xlab = "Attitude", ylab = "Item consultations", xaxt = "n", yaxt = "n", axes = FALSE)
    axis(1, at = c(-1, 0, 1), col = "white", col.ticks = "black")
    axis(2, at = (1:ksteps)* trunc(steps / ksteps), col = "white", col.ticks = "black")
    lines(c(0,0), c(0, steps), col = "black")
  }  
  for (i in 1:steps){
    ll = stepPAT(bp, bn, prejud, beta)
    if (visu){
      for(ii in 1:nta){
        xx = ii/nta
        #if(bn[ii] == ll$bn[ii]) lines(c(-xx, -xx), c(i-1, i), col = cl[bn[ii]+1], lwd = 2)
        if(bn[ii] == ll$bn[ii] & bn[ii] == 1) points(-xx, i, col = cl[bn[ii]+1], pch = 16, cex = 0.2)
        #if(bp[ii] == ll$bp[ii]) lines(c(xx, xx), c(i-1, i), col = cl[bp[ii]+1], lwd = 2)
        if(bp[ii] == ll$bp[ii] & bp[ii] == 1) points(xx, i, col = cl[bp[ii]+1], pch = 16, cex = 0.2)
      }
      points(ll$arg/nta, i, pch = pch[ll$bel+1], col = 'darkgreen', lwd = 2)
    }
    bp = ll$bp
    bn = ll$bn
    bb = sum(bp) - sum(bn)
    aa = exp(- beta * (prejud + sum(bp) - sum(bn)))
    # vv = exp(-  (prejud + sum(bp) - sum(bn)))
    aa = (1-aa)/(1+aa)
    traja = c(traja, aa)
    trajb = c(trajb, bb)
    if(visu) lines(c(traja[i], aa), c(i-1, i), col = "black", lwd = 2)
  }
  return(list('attitude' = traja, "belief" = trajb))
}


#-------------------------
# returns distribution of attitudes from trajectories
# nTraj = nb of trajectories
PATDistribRun = function(steps, nTraj, N, beta = 0.15, prejud = 0, visu = FALSE){
  
  nbs = rep(0, N+1)
  na = N /2
  svm = matrix(rep(0, (N+1)*(steps)), nrow = N+1, ncol = steps)
  colnames(svm) = 1:steps
  avals = exp(- beta*(-(na):na))
  avals = (1-avals)/(1+avals)
  rownames(svm) = avals
  for(i in 1:nTraj){
    cat("\r   ",i)
    trajs = trajPAT(steps,N, beta)
    for(j in 1:steps){
      k = trajs$belief[j] + na+1 
      svm[k, j] = svm[k, j]+1
    }
  }
  svm = svm / nTraj
  if (visu) visuMatKSteps(sv = svm, ksteps = 5, decr = TRUE)
  return(svm)
}

#=================
# Argument model in Markov chain version
#=================
runPATDens = function(steps, N, beta = 0.15, visu = TRUE){
  
  if (N %% 2 != 0) N = N+1 # N should be even
  na = N/2
  sm = matrix(data = rep(0, (na+1)^2), nrow = na+1, ncol = na+1)
  sm[1,1] = 1
  #  print(sm[51,51])
  st = matrix(data = rep(0, steps*(N+1)), nrow = N+1, ncol = steps)
  colnames(st) = 1:steps
  avals = exp(- beta*(-(na):na))
  avals = (1-avals)/(1+avals)
  rownames(st) = avals
  for (tt in (1:steps)){
    cat("\r   ",tt)
    delta = matrix(data = rep(0, (na+1)^2), nrow = na+1, ncol = na+1)
    for (p in (0:na)){
      for (n in (0:na)){
        prb = 1 / (1 + exp(- beta * (p -n)))
        # print(paste("p: ", p, "  n: ", n))
        dd = sm[p+1, n+1]
        if (dd > 0){
          if (p < na){ 
            dp = dd * prb * (1 - p/na)/2
            delta[p+2, n+1] = delta[p+2, n+1] + dp
            delta[p+1, n+1] = delta[p+1, n+1] -dp
          }
          if (p > 0) { 
            dp = dd * (1 - prb) * p/ (2*na)
            delta[p, n+1] = delta[p, n+1] + dp
            delta[p+1, n+1] = delta[p+1, n+1] -dp
          }
          if (n < na){ 
            dp = dd * (1-prb) * (1 - n/na)/2
            delta[p+1, n+2] =  delta[p+1, n+2] + dp
            delta[p+1, n+1] = delta[p+1, n+1] -dp
          } 
          if (n > 0){ 
            dp = dd * prb * n/ na /2
            delta[p+1, n] = delta[p+1, n] + dp
            delta[p+1, n+1] = delta[p+1, n+1] -dp
          } 
        }
      }
    }
    sm = sm + delta
    states = rep(0, N+1)
    for (p in (0:na)){
      for (n in (0:na)){
        ii = p - n + na+1
        states[ii] = states[ii]+sm[p+1, n+1]
      }
    }
    st[,tt] = states
  } 
  if (visu) visuMatKSteps(sv = st, ksteps = 5, decr = TRUE)
  return(st)
}

#=================
# Visualisation of distribution evolution (for both models)
#=================

#Associates colour to extremity
#--------------------
colExtr = function(extr){
  
  rd = c(1, 0, 0)
  gr = c(0, 1, 0)
  or = c(1, 165/255, 0)
  if (extr <= 0.5) {
    coef = extr*2
    vals = gr * (1 - coef) + or * coef
  } 
  else{
    coef = (extr - 0.5)*2
    vals = or*(1-coef)+rd *coef
  }
  return(rgb(vals[1], vals[2], vals[3]))
}


#Visualises the evolution of distribution from matrix sm
# It represents ksteps distributions at regular intervals
# decr TRUE the bars width decreases when extremity increases
#----------------
visuMatKSteps = function(sv, ksteps, lwd = 2, top = FALSE, decr = FALSE){
  
  sm = ncol(sv)
  intSteps = trunc(sm / ksteps)
  #print(intSteps)
  vals = as.numeric(rownames(sv))
  par(mar = c(4, 4, 0.1, 0.1), cex = 1.2)
  if (top){
    #par(mar = c(0.1, 4, 4, 0.1), cex = 1.2)
    plot(1,-100, xlim = c(-1,1), ylim = c(sm, 0), xlab = "Attitude", ylab = "Item consultations", xaxt = "n", yaxt = "n", axes = FALSE)
  }
  else {
    plot(1,-1, xlim = c(-1,1), ylim = c(intSteps, sm+intSteps), xlab = "Attitude", ylab = "Item consultations", xaxt = "n", yaxt = "n", axes = FALSE)}
  axis(1, at = c(-1, 0, 1), col = "white", col.ticks = "black")
  axis(2, at = (1:ksteps)*intSteps, col = "white", col.ticks = "black")
  mx = max(sv[,(1:ksteps)*intSteps])*0.9
  if (top) mult = -1 else mult = 1
  if (decr) {
    fac = 0.9
    dd = 1
  }
  else {
    fac = 1
    dd = 0
  }
  for (ss in 1:ksteps){
    yy = ss*intSteps
    mx = max(sv[,yy])*1.4
    dist = sv[,ss*intSteps]*mult
    for(i in 1:length(dist)) 
      if (dist[i] > 0) lines(c(vals[i], vals[i]), c(yy,yy + dist[i]*intSteps/mx), col = colExtr(abs(vals[i])), lwd = lwd*(1-abs(vals[i]*dd)*fac))
    lines(c(-1,1), c(yy, yy), col = "grey", lwd = 2)
  }
}


