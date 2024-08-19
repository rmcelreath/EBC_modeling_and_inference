# sythnetic simulation of Pan brain/vocal data

library(rethinking)

source("simulation_functions.r")

############################################
# combine brain and vocal models

# scaffolded model formulas

formulaB <- alist(
    # B model
     brain_obs ~ dgampois(mu,tau),
     mu <- B[ID]^theta * brain_base ,
     save> vector[n_individuals]:B <- inv_logit(beta[S]*(A-alpha) + brain_z[ID]),
     brain_z[ID] ~ normal( 0 , 1 ),
    # priors
     delta ~ exponential(1),
     vector[2]:gamma ~ normal(0,1),
     alpha ~ normal(5,1),
     vector[2]:beta ~ normal(0,1),
     tau ~ exponential(1),
     theta ~ exponential(1),
     bXV ~ normal(0,0.5),
     bXB ~ normal(0,0.5)
)

formulaV <- alist(
    # Y > 0
    Y|Y>0 ~ custom( log(p[BID]*( 1-exp(-bA*A[IDV[i]])*exp(-bB*B[IDV[i]]) )) + poisson_log_lpmf(Y|log_lambda1) ),
    # Y = 0
    Y|Y==0 ~ custom( log_sum_exp( 
        log(p[BID]*( 1-exp(-bA*A[IDV[i]])*exp(-bB*B[IDV[i]]) )) + poisson_log_lpmf(0|log_lambda1) , 
        log1m(p[BID]*( 1-exp(-bA*A[IDV[i]])*exp(-bB*B[IDV[i]]) ))
        ) ),
    #c(bA,bB) ~ normal(0,0.5),
    bA ~ half_normal(0,0.5),
    bB ~ half_normal(0,0.5),
    # matrix[n_individuals,M]:Q <- p[i] * ( 1-exp(-bA*A[BID])*exp(-bB*B) ) ,
    # log-linear model for rate of behavior
    log_lambda1 <- log(L[BID]) + log(D[J]),
    # priors
    vector[M]:L ~ exponential(1),
    vector[M]:p ~ beta(2,2),
    vector[n_individuals]:B ~ beta(4,4)
)

formulaV2 <- alist(
    # Y > 0
    Y|Y>0 ~ custom( log(p[BID]) + poisson_log_lpmf(Y|log_lambda1) ),
    # Y = 0
    Y|Y==0 ~ custom( log_sum_exp( 
        log(p[BID]) + poisson_log_lpmf(0|log_lambda1) , 
        log1m(p[BID])
        ) ),
    #c(bA,bB) ~ normal(0,0.5),
    bA ~ half_normal(0,0.5),
    bB ~ half_normal(0,0.5),
    # matrix[n_individuals,M]:Q <- p[i] * ( 1-exp(-bA*A[BID])*exp(-bB*B) ) ,
    # log-linear model for rate of behavior
    log_lambda1 <- log(L[BID]) + log(D[J]),
    # priors
    vector[M]:L ~ exponential(1),
    vector[M]:p ~ beta(2,2)
)

formulaV3 <- alist(
    # Y > 0
    Y|Y>0 ~ custom( log(p[BID]*( 1-exp(-bA*A[IDV[i]]) )) + poisson_log_lpmf(Y|log_lambda1) ),
    # Y = 0
    Y|Y==0 ~ custom( log_sum_exp( 
        log(p[BID]*( 1-exp(-bA*A[IDV[i]]) )) + poisson_log_lpmf(0|log_lambda1) , 
        log1m(p[BID]*( 1-exp(-bA*A[IDV[i]]) ))
        ) ),
    bA ~ half_normal(0,0.5),
    # matrix[n_individuals,M]:Q <- p[i] * ( 1-exp(-bA*A[BID])*exp(-bB*B) ) ,
    # log-linear model for rate of behavior
    log_lambda1 <- log(L[BID]) + log(D[J]),
    # priors
    vector[M]:L ~ exponential(1),
    vector[M]:p ~ beta(2,2)
)

formulaVB <- alist(
    # Y > 0
    Y|Y>0 ~ custom( log_pm[i] + poisson_log_lpmf(Y|log_lambda1) ),
    # Y|Y>0 ~ custom( log(p[BID]*( 1-exp(-bA*A[IDV[i]])*exp(-bB*B[IDV[i]]) )) + poisson_log_lpmf(Y|log_lambda1) ),
    # Y = 0
    Y|Y==0 ~ custom( log_sum_exp( 
        log_pm[i] + poisson_log_lpmf(0|log_lambda1) , 
        log1m_exp(log_pm[i])
        ) ),
    #Y|Y==0 ~ custom( log_sum_exp( 
    #    log(p[BID]*( 1-exp(-bA*A[IDV[i]])*exp(-bB*B[IDV[i]]) )) + poisson_log_lpmf(0|log_lambda1) , 
    #    log1m(p[BID]*( 1-exp(-bA*A[IDV[i]])*exp(-bB*B[IDV[i]]) ))
    #    ) ),
    #c(bA,bB) ~ normal(0,0.5),
    bA ~ half_normal(0,0.5),
    bB ~ half_normal(0,0.5),
    # matrix[n_individuals,M]:Q <- p[i] * ( 1-exp(-bA*A[BID])*exp(-bB*B) ) ,
    log_pm <- log( p[BID]*( 1-exp(-bA*A[IDV[i]])*exp(-bB*B[IDV[i]]) ) ),
    # log-linear model for rate of behavior
    log_lambda1 <- log(L[BID]) + log(D[J]),
    # priors
    vector[M]:L ~ exponential(1),
    vector[M]:p ~ beta(2,2),

    # B model
     brain_obs ~ dgampois(mu,tau),
     mu <- exp( theta*log(B[ID]) + log(brain_base) ),
     save> vector[n_individuals]:B <- inv_logit(beta[S]*(A-alpha) + brain_z[ID]),
     brain_z[ID] ~ normal( 0 , 1 ),
    # priors
     delta ~ exponential(1),
     vector[2]:gamma ~ normal(0,1),
     alpha ~ normal(5,1),
     vector[2]:beta ~ normal(0,1),
     tau ~ exponential(1),
     theta ~ exponential(1),
     bXV ~ normal(0,0.5),
     bXB ~ normal(0,0.5)
)

dat <- sim_brains(n_individuals=100)

datV <- sim_repertoire(N=dat$n_individuals,M=20,NR=10,age=dat$A,brain=dat$brain_true,beta=c(0.1,1))
datV$a <- NULL

# ulam version with post-sampling computaiton of repertoires
# convert Y matrix to long form with m column
YY <- rep(NA,length(datV$Y))
MM <- YY
ID <- YY
JJ <- YY
k <- 1
for ( j in 1:nrow(datV$Y) ) for ( m in 1:ncol(datV$Y) ) {
    YY[k] <- datV$Y[j,m]
    MM[k] <- m
    ID[k] <- datV$id[j]
    JJ[k] <- j
    k <- k + 1
}

dat_merge <- dat
dat_merge$Y <- YY
dat_merge$IDV <- ID
dat_merge$BID <- MM
dat_merge$J <- JJ
dat_merge$D <- datV$d
dat_merge$M <- max(MM)

m3 <- ulam(
    formulaVB,
    data=dat_merge, 
    constraints=list(delta="lower=0",gamma="lower=0",alpha="lower=0"),
    chains=3 , cores=3 , sample=TRUE )

precis(m3)

post <- extract.samples(m3)

# vocalization p inference
plot( datV$p , apply(post$p,2,mean) , xlab="p (true)" , ylab="p (posterior mean)" , xlim=c(0,1) , ylim=c(0,1) )
ci <- apply((post$p),2,PI,prob=0.5)
for ( i in 1:datV$M ) {
    lines( rep(datV$p[i],2) , ci[,i] , col=grau(0.5) )
}
abline(a=0,b=1,lty=2)

# brain inference
plot( dat$brain_true , apply(post$B,2,mean) , xlab="brain (true)" , ylab="brain (posterior mean)" )
ci <- apply((post$B),2,PI,prob=0.5)
for ( i in 1:dat$n_individuals ) {
    lines( rep(dat$brain_true[i],2) , ci[,i] , col=grau(0.5) )
}
abline(a=0,b=1,lty=2)


############################################
# now with missing data

# Brain missing as function of age
RB <- rbern(n_individuals,prob=inv_logit(age*0.1))
Bo <- brain_obs/base
Bo[RB==1] <- NA # missing

# Vocal missing as random, 10% missing
RV <- rbern(n_individuals,prob=ifelse(RB==0,0.8,0.1))
Vo <- vocal
Vo[RV==1] <- NA # missing

# cbind( Bo , Vo )
table( RB , RV )

dat2 <- list(
    n_individuals = n_individuals,
    A = age,
    S = sex,
    brain_obs = Bo,
    V = ifelse( is.na(Vo) , -10 , Vo ),
    X = X
)

m2 <- ulam(
    formula,
    data=dat2, 
    constraints=list(delta="lower=0",gamma="lower=0",alpha="lower=0"),
    chains=3 , cores=3 , sample=TRUE )

#blank(w=3,h=0.8)
yoff <- (-0.3)

par(mfrow=c(2,1))

plot(precis(m1,2,pars=c("gamma","beta","delta")))
mtext("no missing data")
points( delta , 1+yoff , col=2 , pch=3 , lwd=3 )
points( b_rate[2] , 2+yoff , col=2 , pch=3 , lwd=3)
points( b_rate[1] , 3+yoff , col=2 , pch=3 , lwd=3)
points( v_rate[2] , 4+yoff , col=2 , pch=3 , lwd=3)
points( v_rate[1] , 5+yoff , col=2 , pch=3 , lwd=3)

plot(precis(m2,2,pars=c("gamma","beta","delta")),xlim=c(0,2))
mtext("missing data")
points( delta , 1+yoff , col=2 , pch=3 , lwd=3 )
points( b_rate[2] , 2+yoff , col=2 , pch=3 , lwd=3)
points( b_rate[1] , 3+yoff , col=2 , pch=3 , lwd=3)
points( v_rate[2] , 4+yoff , col=2 , pch=3 , lwd=3)
points( v_rate[1] , 5+yoff , col=2 , pch=3 , lwd=3)
