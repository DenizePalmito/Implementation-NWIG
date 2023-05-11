#_________________________________________________________________________________________________________
#                          Upload the following packages
#_________________________________________________________________________________________________________
require(gamlss)
require(numDeriv)
#
#_________________________________________________________________________________________________________
#               Next New Weibull inverse Gaussian (NWIG) implementation in gamlss 
#_________________________________________________________________________________________________________
N_Weibull_IG<- function (mu.link = "log", sigma.link="log", nu.link = "identity", tau.link = "identity")
{
  mstats <- checklink("mu.link", "New Weibull inverse Gaussian", substitute(mu.link), 
                      c("inverse", "log", "identity", "logit"))
  dstats <- checklink("sigma.link", "New Weibull inverse Gaussian", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink(   "nu.link", "New Weibull inverse Gaussian", substitute(nu.link),    
                         c("1/nu^2", "log", "identity", "own"))
  tstats <- checklink(  "tau.link", "New Weibull inverse Gaussian", substitute(tau.link),   
                        c("1/tau^2", "log", "identity", "own")) 
  structure(
    list(family = c("N_Weibull_IG", "New Weibull inverse Gaussian"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
         nopar = 4, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         nu.link = as.character(substitute(nu.link)), 
         tau.link = as.character(substitute(tau.link)), 
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         nu.linkfun = vstats$linkfun,
         tau.linkfun = tstats$linkfun,  
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         tau.linkinv = tstats$linkinv, 
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         nu.dr = vstats$mu.eta,
         tau.dr = tstats$mu.eta, 
#_________________________________________________________________________________________________________
#                          First, Second and Cross Numerical Derivatives 
#_________________________________________________________________________________________________________
         dldm = function(y,mu,sigma,nu,tau){  
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiN_Weibull_IG(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           dldm
         },
         d2ldm2 = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiN_Weibull_IG(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           d2ldm2 <- -dldm * dldm
           d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
           d2ldm2 
         },     
         dldd = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,x,nu,tau){log(dauxiN_Weibull_IG(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
           dldd
         } ,
         d2ldd2 = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,x,nu,tau){log(dauxiN_Weibull_IG(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
           d2ldd2 <- -dldd*dldd
           d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
           d2ldd2
         },   
         dldv = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiN_Weibull_IG(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau,method='simple')
           dldv
         },
         d2ldv2 = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiN_Weibull_IG(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau,method='simple')
           d2ldv2<- -dldv * dldv
           d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                   
           d2ldv2
         },
         dldt = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiN_Weibull_IG(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           dldt
         } ,
         d2ldt2 = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiN_Weibull_IG(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau,method='simple')
           d2ldt2<- -dldt * dldt
           d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
           d2ldt2
         },
         d2ldmdd = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiN_Weibull_IG(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,x,nu,tau){log(dauxiN_Weibull_IG(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau)
           d2ldmdd = -(dldm * dldd)
           d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
           d2ldmdd                 
         },
         d2ldmdv = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiN_Weibull_IG(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiN_Weibull_IG(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
           d2ldmdv = -(dldm * dldv)
           d2ldmdv				
         },
         
         d2ldmdt = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiN_Weibull_IG(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiN_Weibull_IG(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           d2ldmdt <- -(dldm*dldt)
           d2ldmdt
         },
         
         d2ldddv = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,x,nu,tau){log(dauxiN_Weibull_IG(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiN_Weibull_IG(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
           d2ldddv = -(dldd * dldv)
           d2ldddv	
         },
         d2ldddt = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,x,nu,tau){log(dauxiN_Weibull_IG(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau)
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiN_Weibull_IG(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           d2ldddt <- -(dldd*dldt) 
           d2ldddt 
         },
         d2ldvdt = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiN_Weibull_IG(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiN_Weibull_IG(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           d2ldvdt <- -(dldv*dldt) 
           d2ldvdt 
         },
         
         G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
         { 
           -2*dN_Weibull_IG(y,mu,sigma,nu,tau,log=TRUE)
         } ,                     
         rqres = expression(   
           rqres(pfun="pN_Weibull_IG", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
         
         mu.initial = expression( mu <- rep(0.5, length(y))),
         sigma.initial = expression( sigma <-rep(0.5, length(y))),
         nu.initial = expression(nu <- rep(0.5, length(y))), 
         tau.initial = expression(tau <-rep(0.5, length(y))), 
         
         mu.valid = function(mu) all(mu> 0), 
         sigma.valid = function(sigma) all(sigma > 0),
         nu.valid = function(nu) all(nu > 0), 
         tau.valid = function(tau) all(tau > 0), 
         y.valid = function(y) all(y>0)),
    class = c("gamlss.family","family"))
}
#_________________________________________________________________________________________________________
#                          NWIG probability density function
#_________________________________________________________________________________________________________
dN_Weibull_IG <- function(x, mu = 0.5, sigma = 0.5, nu = 0.5, tau = 0.5, log = FALSE){
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))
  
  pdfsn <- dIG(x,mu=mu,sigma=sigma)
  cdfsn <- pIG(x,mu=mu,sigma=sigma)
  
  fy1 <- ((nu*tau*pdfsn)/(cdfsn))*((-log(cdfsn))^(tau-1))*(exp(-nu*(-log(cdfsn))^(tau)))
  
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}    
#_________________________________________________________________________________________________________
#                          NWIG cumulative density function 
#_________________________________________________________________________________________________________
pN_Weibull_IG <- function(q, mu = 0.5, sigma = 0.5, nu = 0.5, tau = 0.5, lower.tail = TRUE, log.p = FALSE){  
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))
  if (any(q < 0))  stop(paste("q must be positive", "\n", ""))          
  
  cdfsn<-pIG(q,mu=mu,sigma=sigma)
  
  cdf1 <- exp(-nu*(-log(cdfsn))^(tau))
  
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#_________________________________________________________________________________________________________
#                          NWIG quantile function 
#_________________________________________________________________________________________________________
qN_Weibull_IG <-  function(p, mu = 0.5, sigma = 0.5, nu = 0.5, tau = 0.5, lower.tail = TRUE, log.p = FALSE){   
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))
  
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  
  u <- exp(-(log(p)/(-nu))^(1/tau))
  
  q <- qIG(u, mu=mu , sigma=sigma)
  q
}
#_________________________________________________________________________________________________________
#                          NWIG random generating function
#_________________________________________________________________________________________________________
rN_Weibull_IG <- function(n, mu = 0.5, sigma = 0.5, nu = 0.5, tau = 0.5){
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
  if (any(n < 0))  stop(paste("n must be positive", "\n", "")) 
  
  uni<- runif(n = n,0,1)
  
  r <- qN_Weibull_IG(uni,mu =mu, sigma =sigma, nu=nu, tau=tau)
  r
}
#_________________________________________________________________________________________________________
#                          NWIG auxiliar function for derivates 
#_________________________________________________________________________________________________________
dauxiN_Weibull_IG <- function(t,mu,sigma,nu,tau){ 
  pdfsn <- dIG(t,mu=mu,sigma=sigma)
  cdfsn <- pIG(t,mu=mu,sigma=sigma)
  
  fy1 <- ((nu*tau*pdfsn)/(cdfsn))*((-log(cdfsn))^(tau-1))*(exp(-nu*(-log(cdfsn))^(tau)))
  fy1}
#_________________________________________________________________________________________________________
#                          End of NWIG implemented in gamlss
#_________________________________________________________________________________________________________