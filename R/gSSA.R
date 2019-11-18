#' stochastic simulation function
#' @include RcppExports.R
#' @include getBeta.R
#' @export gSSA

gSSA <- function(N, R0, gamma, sigma, rho, epsilon, omega, c, m, omegam, mu, mj, phi, s, 
                 jspan, span, k, preg, Ea0, Ej0, En0, Ra0, Rj0, Rn0, Na, Nj, Nn,
                 whichBeta, whichRho, whichGamma, tmax=10, inc=52, maxtau=1/365, nevents=56) {

  #calculate beta
  beta <- max(0, getBeta(R0, rho, epsilon, sigma, m, gamma, N, whichBeta))
  
  #get approximate equilibrium age class sizes
  #Na <- round(N/(1+m/mu+b/(omegam+mj)))
  #Nj <- round(Na*m/mu)
  #Nn <- N - Na - Nj
  
  #order parameters and specifications for simulation (note that rho and p have been renamed)
  par <- as.numeric(c(beta=beta, gamma=gamma, rho=sigma, p=rho, epsilon=epsilon, omega=omega, c=c, m=m, omegam=omegam, 
           mu=mu, mj=mj, phi=phi, s=s, k=k, preg=preg))
  spec <- as.numeric(c(whichBeta,whichRho,whichGamma))
  
  #set initial population sizes
  Sa0 <- max(0, Na - 10 - Ea0 - Ra0)
  Sj0 <- max(0, Nj - Ej0 - Rj0)
  Sn0 <- max(0, Nn - En0 - Rn0)
  
  init <- as.integer(c(ceiling(Sa0/2),ceiling(Ea0/2),5,ceiling(Ra0/2),max(0, Nn-En0),max(0,Nj-Ej0),
                       Rn0,Ej0,0,Rj0,En0,0,0,floor(Sa0/2),floor(Ea0/2),5,floor(Ra0/2)))
  
  #simulate and store
  sim <- simRcpp(init, par, spec, nevents, tmax, inc, maxtau)
  SF <- sim[1:(tmax*inc)]
  EF <- sim[(tmax*inc+1):(2*tmax*inc)]
  IF <- sim[(2*tmax*inc+1):(3*tmax*inc)]
  RF <- sim[(3*tmax*inc+1):(4*tmax*inc)]
  SJ0 <- sim[(4*tmax*inc+1):(5*tmax*inc)]
  SJ <- sim[(5*tmax*inc+1):(6*tmax*inc)]
  MJ0 <- sim[(6*tmax*inc+1):(7*tmax*inc)]
  EJ <- sim[(7*tmax*inc+1):(8*tmax*inc)]
  IJ <- sim[(8*tmax*inc+1):(9*tmax*inc)]
  RJ <- sim[(9*tmax*inc+1):(10*tmax*inc)]
  EJ0 <- sim[(10*tmax*inc+1):(11*tmax*inc)]
  IJ0 <- sim[(11*tmax*inc+1):(12*tmax*inc)]
  RJ0 <- sim[(12*tmax*inc+1):(13*tmax*inc)]
  SM <- sim[(13*tmax*inc+1):(14*tmax*inc)]
  EM <- sim[(14*tmax*inc+1):(15*tmax*inc)]
  IM <- sim[(15*tmax*inc+1):(16*tmax*inc)]
  RM <- sim[(16*tmax*inc+1):(17*tmax*inc)]  
  S <- SF + SM
  E <- EF + EM
  I <- IF + IM
  R <- RF + RM
  
  Sall <- S + SJ0 + SJ
  Eall <- E + EJ0 + EJ
  Iall <- I + IJ0 + IJ
  Rall <- R + RJ0 + RJ + MJ0
  
  #any <- which(Eall > 0 | Iall > 0)
  #S <- S[any]
  #E <- E[any]
  #I <- I[any]
  #R <- R[any]
  #SF <- SF[any]
  #EF <- EF[any]
  #IF <- IF[any]
  #RF <- RF[any]
  #SM <- SM[any]
  #EM <- EM[any]
  #IM <- IM[any]
  #RM <- RM[any]
  #SJ <- SJ[any]
  #EJ <- EJ[any]
  #IJ <- IJ[any]
  #RJ <- RJ[any]    
  #SJ0 <- SJ0[any]
  #EJ0 <- EJ0[any]
  #IJ0 <- IJ0[any]
  #RJ0 <- RJ0[any]
  #MJ0 <- MJ0[any]
  #Sall <- Sall[any]
  #Eall <- Eall[any]
  #Iall <- Iall[any]
  #Rall <- Rall[any]
  N <- Sall + Eall + Iall + Rall
  Na <- S + E + I + R
  Nj <- SJ + EJ + IJ + RJ
  Nn <- SJ0 + EJ0 + IJ0 + RJ0 + MJ0
  
  return(list(#S=S,E=E,I=I,R=R,SJ=SJ,EJ=EJ,IJ=IJ,RJ=RJ,
              #SJ0=SJ0,EJ0=EJ0,IJ0=IJ0,RJ0=RJ0,MJ0=MJ0,
              #Sall=Sall,Eall=Eall,Iall=Iall,Rall=Rall,N=N,
              #Na=Na,Nj=Nj,Nn=Nn,SF=SF,EF=EF,IF=IF,RF=RF,
              #SM=SM,EM=EM,IM=IM,RM=RM,
              PosA=IF+IM+EF+EM+RF+RM,
              NegA=SF+SM,
              PosJ=IJ+IJ0+EJ0+EJ+RJ0+RJ+MJ0,
              NegJ=SJ+SJ0,
              I=IJ+IJ0+IM+IF,
              EA=EF+EM,
              EJ=EJ,
              EN=EJ0))
}