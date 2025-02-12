################################################################################
t0 <- Sys.time()
require(runjags)
require(dplyr)
################################################################################
# data
load(file="int.data_240826.RData")
################################################################################
# JAGS Model Settings
# Parameters monitored
parameters <- c("theta", "PSR", 
								"gamma0", "gamma1", 
                "beta0", "beta1", "beta2", "beta3",
                "alpha0", "alpha1", #"alpha2", "alpha3",
                "p2","p3","eps2", "eps3",
								"p.f.z1 ", "p.m.z1", # proportion spawners by sex
								"PSI.male", "PSI.female", # mean trans prob
								"P.male.ret", "P.female.ret", # average interval probabilities
								"barI.male", "barI.female", # average intervals
                "deviance",
								"T1zero", "T1zero.rep", "T1zero.sim",
								#"T1one", "T1one.rep", "T1one.sim", # cannot observe y=1 in y_1
								"T2zero", "T2zero.rep", "T2zero.sim",
								"T2one", "T2one.rep", "T2one.sim",
								"T3zero", "T3zero.rep", "T3zero.sim",
								"T3one", "T3one.rep", "T3one.sim")	
# MCMC parameters
nchains <- 3 # number of MCMC chains
nadapt <- 1000 # adaption phase
nburn <- 30000 # discard these draws
niter <- 20000 # length of MCMC chains
nsamp <- 10000 # number of samples to take from each MCMC chains
thin_ <- round(niter/nsamp)  
################################################################################
# Functions for JAGS model evaluation
## Initialize unknown z states function (altered slightly from K&S p276)
	# sample with equal probability from among possible states for all cases where 
  # state is unknown
	ms.init.z <- function(ch, f, stoch){
		for (i in 1:nrow(ch)){ch[i,1:f[i]] <- NA} # for each i, first capture is NA
		states <- max(ch, na.rm = TRUE) # largest state, equals number of states
		known.states <- 1:(states-1) # vector of known states (w/o unobserved/max)
		v <- which(ch==states|ch==stoch) # occs w/ unobserved or stochastic state
		ch[-v] <- NA #when state observed, equals NA (coming in as data)
		ch[v] <- sample(known.states, length(v), replace = TRUE)  # draw rand states
		return(ch)
	}
## Initialize uninformed G.prior with G.prior=0.5
  init.G.prior <- function(G.prior, nind.g.na, inds.g.na){ 
    G.prior <- 
  for (i in 1:nind.g.na){ G.prior[inds.g.na[i]] <- 0.5} 
  }
## Initialize JAGS model	
	inits <- function(chains=1, dat=int.data){ 
		attach(int.data)
		inits.out <- list()
		for (c in 1:chains){
			inits.i <- list()
			inits.i$theta = rbeta(1, 3, 3)
			inits.i$gamma0 = rnorm(1, 0, sqrt(2))
			inits.i$gamma1 = rnorm(1, 0, sqrt(2))
			inits.i$beta0 = rnorm(2, 0, sqrt(2))
			inits.i$beta1 = rnorm(2, 0, sqrt(2))
			inits.i$beta2 = rnorm(2, 0, sqrt(2))
			inits.i$beta3 = rnorm(2, 0, sqrt(2))
			inits.i$alpha0 = rnorm(2, 0, sqrt(2))
			inits.i$alpha1 = rnorm(2, 0, sqrt(2))
      inits.i$eps2 = runif(1, 0, 1)	
      inits.i$eps3 = runif(1, 0, 1)
			inits.i$p2 = array(runif(4, 0, 1), c(2,2))
			inits.i$p3 = runif(2, 0, 1)
			inits.i$z = ms.init.z(dat$min_y_it, f, 2)
			inits.i$.RNG.name = "base::Wichmann-Hill"
			inits.i$.RNG.seed = runif(1, 1, 1e3)
			inits.out[[c]] <- inits.i
		}
		detach(int.data)
		return(inits.out)
	}		
################################################################################
# Run models
runfun <- function(model, keepdir){
  run.jags(data=int.data, model=model, 
           monitor=parameters, n.chains=nchains, thin=thin_, adapt=nadapt,
           burnin=nburn, sample=nsamp, inits=inits(chains=nchains,int.data), 
           summarise=F, keep.jags.files=keepdir, method='parallel')
}
	
# null temperature model
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3_null.txt", keepdir="runjagsfiles9.0.3_null"))

# X.DOY5
# additive
int.data$X <- int.data$X.DOY5%>% scale() %>% as.vector()
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3_additive.txt", keepdir="runjagsfiles9.0.3_DOY5_additive"))
# interactive
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3.txt", keepdir="runjagsfiles9.0.3_DOY5_interaction"))

# X.T100
# additive
int.data$X <- int.data$X.T100%>% scale() %>% as.vector()
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3_additive.txt", keepdir="runjagsfiles9.0.3_T100_additive"))
# interactive
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3.txt", keepdir="runjagsfiles9.0.3_T100_interaction"))

# X.icex
# additive
int.data$X <- int.data$X.icex%>% scale() %>% as.vector()
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3_additive.txt", keepdir="runjagsfiles9.0.3_icex_additive"))
# interactive
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3.txt", keepdir="runjagsfiles9.0.3_icex_interaction"))

# X.Dgt5lag1
# additive
int.data$X <- int.data$X.Dgt5lag1 %>% scale() %>% as.vector()
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3_additive.txt", keepdir="runjagsfiles9.0.3_Dgt5lag1_additive"))
# interactive
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3.txt", keepdir="runjagsfiles9.0.3_Dgt5lag1_interaction"))

# X.DOY5lag1
# additive
int.data$X <- int.data$X.DOY5lag1%>% scale() %>% as.vector()
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3_additive.txt", keepdir="runjagsfiles9.0.3_DOY5lag1_additive"))
# interactive
set.seed(999)
system.time(ms0 <- runfun(model="IMSCJS_m9.0.3.txt", keepdir="runjagsfiles9.0.3_DOY5lag1_interaction"))

