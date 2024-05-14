library(combinat)
cycle = function(x)
{
  x = c(x[2:length(x)],x[1])
  return(x)
}

one = function(t) matrix(1, nrow = t)

initial.design = function(n,m)
{
  if(m < 6) 
  {
    temp = unlist(permn(m))
    d = matrix(temp, ncol = m, byrow = TRUE)
    d = d[sample(factorial(m),size = n),]
    # rm(temp)
  } else d = t(replicate(n, sample(1:m, m)))
  return(d)
}

pwo = function(x)
{
  m = length(x) #x is a numeric vector with integer numbers
  O = vector(mode = "numeric", length = m*(m-1)/2)
  cnt = 0
  for(i in 1:(m-1))
  {
    for(j in (i+1):m)
    {
      cnt = cnt + 1
      p1 = which(x == i)
      p2 = which(x == j)
      if(p1 < p2) O[cnt] = 1 else O[cnt] = -1
    }
  }
  return(O)
}

swap = function(x, i)
{
  temp = x[i]
  if(i < length(x)) 
  {
    x[i] = x [i+1] 
    x[i+1] = temp
  }  else {
            x[i] = x[1]
            x[1] = temp
          }  
  return(x)
}

doofa.pwo = function(n,m)
{
  p = m*(m-1)/2 + 1
  detM = 0
  try.init = 0
  while(detM < 10^(-10) & try.init <= 10000)
  {
    try.init = try.init + 1
    d0 = initial.design(n,m)
    P0 = t(apply(d0, 1, pwo)) 
    x = cbind(one(n), P0)
    xpx = t(x)%*%x
    M = xpx/n
    detM = det(M)
    if(try.init == 10000 & detM < 10^(-10)) stop("try.init maximum limit reached. try increasing try.init value")
  }
  deff0 = (detM/(((m+1)^(m-1))/3^(p-1)))^(1/p)
  repeat
  {
    deff = deff0
    for(i in 1:n)
    {
      di = d0[i,]
      pi = P0[i,]
      if(i > 1) 
      {
          p1 = P0[1:(i-1),] 
          dim(p1) = c(i-1,ncol(P0))
      }
      if(i < n) 
      {
        p2 = P0[(i+1):n,] 
        dim(p2) = c(n-i,ncol(P0))
      }
      if(i > 1 & i < n) 
      {
           a = t(p1)%*%one(nrow(p1)) + pi + t(p2)%*%one(nrow(p2)) 
           A = t(p1)%*%p1 + t(p2)%*%p2 + pi%*%t(pi)
      } 
      if(i == 1) 
      {
        a =  pi + t(p2)%*%one(nrow(p2)) 
        A =  t(p2)%*%p2 + pi%*%t(pi)
      }  
      if(i == n)
      {
        a = t(p1)%*%one(nrow(p1)) + pi
        A = t(p1)%*%p1 + pi%*%t(pi)
      }      
      for(j in 1:m)
      {
        di1 = swap(di,j)
        pistar = pwo(di1)
        astar = a - pi + pistar
        Astar = A - pi%*%t(pi) + pistar %*% t(pistar)
        xpxstar = rbind(cbind(n,t(astar)),cbind(astar,Astar))
        Mstar = xpxstar/n
        detMstar = det(Mstar)
        if(detMstar < 10^(-8)) detMstar = 0
        deffstar = (detMstar/(((m+1)^(m-1))/3^(p-1)))^(1/p)
        if(deffstar > deff) 
        {
          di = di1
          a = astar
          A = Astar
          pi = pistar
          deff = deffstar
        }
      }  
      d0[i,] = di
      P0[i,] = pi
    }
    if(deff <= deff0) break
    deff0 = deff
  }  
  result = list(design = d0, D.eff = deff)    
  return(result)
}  

shuffle = function(x)  
{
  chosen.row = sample(1:nrow(x),1)
  x[chosen.row,] = sample(x[chosen.row,])
  return(x)
}
  
gen.design2= function(n, m, num.repeat = 10)
{
  out = replicate(num.repeat, doofa.pwo(n,m), simplify=FALSE)
  d.efficiencies.out = sapply(out, "[[", 2)
  s = which.max(d.efficiencies.out)
  return(s)
}

bin = function(x,m)
{
  temp = rep(0,m)
  temp[x] = 1
  return(temp)
}

revbin = function(x)
{
  return(which(x == 1))
}

vbin = function(x) as.vector(sapply(x, bin, m = length(x)))

vrevbin = function(x, m) apply(matrix(x, ncol = m), 2, revbin)

library(lpSolve)

getBjk = function(A, Aj, j, k, n, m, lambda, TLMAT)
{
  obj = rep(0, n)	
  constr1 = rep(1, n)  
  constr2 = t(A) 
  if(k > 1) constr3 = t(Aj) else constr3 = NULL
  if(is.null(TLMAT)) constr4 = NULL else constr4 = t(TLMAT)
  constr = rbind(constr1, constr2, constr3, constr4)
  dir1 = c("=")	
  dir2 = rep("=", times = ncol(A))
  if(k > 1) dir3 = rep("=", times = (k-1)) else dir3 = NULL
  if(is.null(TLMAT)) dir4 = NULL else dir4 = rep("<", times = ncol(TLMAT))
  dir = c(dir1, dir2, dir3, dir4)
  rhs1 = lambda*(m-1)
  rhs2 = rep(lambda, ncol(A)) 
  for(col in 1:ncol(A))
  {
    if((col-k)%% m == 0) rhs2[col] = 0
  }
  if(k > 1) rhs3 = rep(0, k-1) else rhs3 = NULL
  if(is.null(TLMAT)) rhs4 = NULL else rhs4 = rep(lambda*(m-1), times = ncol(TLMAT))
  rhs = c(rhs1, rhs2, rhs3, rhs4)
  types = rep("B", times = n)	
  sol = lp (direction = "min", obj, constr, dir, rhs,transpose.constraints = TRUE, all.bin=TRUE, use.rw=TRUE)
  if (sol$status == 0) out = sol$solution else out = 0
  return(out) 
}

onev = function(t) matrix(1,t,1)
getBj = function(j, AM, n, m, lambda, TLM)
{
  AJ = NULL
  for(k in 1:(m-1)) {
    soln = getBjk(AM, AJ, j, k, n, m, lambda, TLM)    
    if(length(soln) > 1) {
      Ajk = soln
      if(k == 1) AJ = Ajk else AJ = cbind(AJ, Ajk) 
    } else {     
      out = 0
      return(out)
    }     
  }
  AJ = cbind(AJ, onev(n)-rowSums(AJ[,1:(m-1)]))  
  out = AJ
  return(out)
}  

alternate.B = function(j1, AMAT, n, m, lambda, TMAT)
{
  success = 0
  try = 0
  ntry = 50*m
  while(success == 0 & try <= ntry)
  {
    try = try + 1  
    r = sample((1:j1),1)
    TMATt = cbind(TMAT, AMAT[,(((r-1)*m+1):(r*m))])
    Bt = AMAT[,-c((((r-1)*m)+1):(r*m))]
    soln1 = getBj(j1, AM = Bt, n, m, lambda, TLM = TMATt)
    if(length(soln1) > 1) {
      TMAT = TMATt
      AMAT = cbind(Bt, soln1)
      success = 1
    }
  }
  if(success == 0) AMAT = 0 
  out = list(AMAT = AMAT, TMAT = TMAT)
  return(out)
}

coa = function(m, lambda, ntrial)
{
  trial = 0
  win = 0
  TL = NULL
  n = lambda*m*(m-1)
  while(trial <= ntrial & win == 0)
  {
    trial = trial + 1    
    d1 = matrix(rep(1:m,lambda*(m-1)), ncol = 1)
    d1 = matrix(sample(d1,n), ncol = 1)
    B = t(apply(d1,1,bin,m))
    j = 2
    while(j <= m)
    {
      soln = getBj(j, AM = B, n, m, lambda, TLM = TL)
      if(length(soln) > 1) {
        B = cbind(B, soln) 
    }  else {
       soln.alt = alternate.B(j1 = (j-1), AMAT = B, n, m, lambda, TMAT = TL)
        if(length(soln.alt$AMAT) > 1) {
          TL = soln.alt$TMAT
          B =  soln.alt$AMAT
          j = j - 1
        } else break
      } 
      j = j + 1
      if(!is.null(TL)) { 
        if (ncol(TL) > 50*m) break
      }
    }
    if(ncol(B) == m*m) {
      win = 1
      D = t(apply(B, 1, vrevbin, m))
      return(D)
    } else TL = B
  }  
  if(win == 0) return("No solution found. You may try again.")
}