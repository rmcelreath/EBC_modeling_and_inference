# sythnetic simulation of Pan brain/vocal data

sim_brains <- function(
    n_individuals = 100 
) {

    # age and sex have no observed influences
    age <- runif(n_individuals,1,20)
    sex <- sample(1:2,size=n_individuals,replace=TRUE)

    # simulate a contextual variable that can influence both B and V
    X <- rbern(n_individuals)
    bXB <- 0
    bXV <- 0

    # brain structure influenced by age and sex
    # assume logistic relationship with age for each sex
    b_rate <- c(0.9,1.1)
    alpha <- 7
    # curve( inv_logit( b_rate[1]*(x-5) ) , from=0 , to=20 , ylim=c(0,1) , xlab="age" , ylab="brain" , lwd=5 )
    brain_offset <- rnorm(n_individuals,0,0.5) # offset for each individual trajectory
    brain_plat <- rnorm(n_individuals,1,0.01)
    brain <- 1 * inv_logit( b_rate[sex]*(age-alpha) + bXB*X + brain_offset )

    if ( FALSE ) {
    # plot trajectories
    plot( NULL , xlim=c(0,20) , ylim=c(0,max(brain_plat)) , xlab="age" , ylab="brain" )
    for ( i in 1:n_individuals )
        curve( brain_plat[i] * inv_logit( b_rate[2]*(x-alpha) + brain_offset[i] ) , from=0 , to=20 , add=TRUE )

    }

    # now observe brain structure - this adds observation error
    base <- 1000
    brain_obs <- rgampois( n_individuals , brain * base , 10 )
    # plot( brain , brain_obs )

    # vocal behavior influenced by age,sex,brain
    # simulate count of unique utterences
    # assume for dramatic effect that brain has NO influence
    v_rate <- c(1,1)
    delta <- 0.5
    lambda <- brain^v_rate[sex] * age^(delta + bXV*X)
    vocal <- rpois(n_individuals, lambda )

    # curve( ( inv_logit( b_rate[1]*(x-5) ) )^v_rate[2] * x^(delta) , from=0 , to=20 , xlab="age" , ylab="vocal" , lwd=5 )

    # plot pairwise associations
    #pairs( ~ brain + vocal + age + sex , col=2 , lwd=2 , cex=1.5 )

    # estimate

    dat <- list(
        n_individuals = n_individuals,
        ID = 1:n_individuals,
        A = age,
        S = sex,
        brain_obs = brain_obs,
        brain_n = brain_obs/base,
        V = vocal,
        X = X,
        brain_true = brain,
        brain_base = base
    )

    return(dat)

}#sim_brains

sim_repertoire <- function(
    N = 50, # number of individuals
    M = 50, # number of unique vocalizations
    NR = 20, # average recordings per individual
    age_model = "DEFAULT",
    beta = c(0,2),
    age,brain,
    p,L,
    X
) {

    a <- NA
    
    if ( missing(p) )
        p <- sort(rbeta(M,5,2)) # prob of each vocalization, ordered rare to common, among individuals
    if ( missing(L) )
        L <- rlnorm(M,0.25,0.5)/10 # rate of each vocalization, when individual uses it

    # simulate individual repertoires
    # individuals in rows, tokens in columns
    x <- matrix(NA,nrow=N,ncol=M)
    if ( age_model=="NULL" )
        x <- t( replicate( N , rbern(M,prob=p) ) )

    # age as individual covariate
    if ( age_model=="DEFAULT" ) {
        # continuous age
        if ( missing(age) )
            a <- runif(N)
        else
            a <- age
        alpha <- p # asymptotic probabilities of each 
        # curve( 0.6*(1-exp(-beta*x)) , from=0 , to=1 , ylim=c(0,1) )
        for ( i in 1:N ) {
            x[i,] <- rbern( M , prob=alpha*( 1 - exp(-beta[1]*a[i])*exp(-beta[2]*brain[i]) ) )
        }#i
    }

    # general categorical covariates, table X supplied for individual covariates, p matrix supplied as argument in which X values already used to generate p values for each individual (set of covariates)
    if ( age_model=="X" ) {
        for ( i in 1:N ) {
            x[i,] <- rbern( M , prob=p[i,] )
        }#i
    }

    # simulate vocal samples with varying sample sizes for each individual
    max_tokens <- 20 # max number of tokens per recording
    n_recordings <- rpois(N,NR) + 1
    id <- rep(1:N,times=n_recordings) # individual making each recording
    y <- matrix(NA,sum(n_recordings),max_tokens)
    d <- 0.2 * rpois(sum(n_recordings),5) + 1 # durations

    if ( FALSE ) {
    # use Gillespie algorithm to simulate vocalizations
    # this gives intervals between tokens, but we don't need that right now
    for ( i in 1:sum(n_recordings) ) {    
        t0 <- 0 # time in recording
        j <- 1
        while( (t0 < d[i]) & (j <= max_tokens) ) {
            r <- runif(2)
            R <- sum(L*x[id[i],]) # total rate of vocalizations
            t1 <- (1/R)*log(1/r[1]) # time to next vocalization
            # which vocalization?
            k <- L*x[id[i],] # rates for this individual
            y[i,j] <- sample(1:M,size=1,prob=k) # draw 1 proportional to relative rates
            # update
            t0 <- t1
            j <- j + 1
            # check for max_tokens breech and truncate recording at this point
            if ( j > max_tokens ) d[i] <- t1
        }#t0 < d
    }#i
    }

    # count number of each token in each recording
    XX <- matrix(NA,sum(n_recordings),M)
    for ( i in 1:sum(n_recordings) ) {
        for ( j in 1:M )
            XX[i,] <- rpois( M , L*x[id[i],]*d[i] )
    }#i

    dat <- list(
        N=N, # number of individuals
        M=M, # number of unique vocalizations
        J=sum(n_recordings), # number of recordings (all individuals)
        d=d, # durations of each recording
        id=id, # individual IDs, length J
        Y=XX, # matrix with recordings on rows and counts of vocalizations in columns
        a=a, # age covariate
        L=L, # rates
        p=p, # base rates
        Rtrue=x ) # true repertoires 

    return(dat)
}#sim_repertoire

# curve( 1 - exp(-2*x) , from=0 , to=1 )

