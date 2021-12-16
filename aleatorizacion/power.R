## Clase sobre poder estad�stico
## Antonella Bandiera
## antonella.bandiera@itam.mx
##
## Realizado para la Red de Aprendizaje del AccLabPY
## 2021-12-13


# Análisis de poder para un estudio con 80 observaciones y 
# efecto de tamaño 0.25
#install.packages("pwr")
library(pwr)
pwr.t.test(n = 40, d = 0.25, sig.level = 0.05,
           power = NULL)



# Análisis de poder para un estudio con 250 observaciones y 
# efecto de tamaño 0.25
library(pwr)
pwr.t.test(n = 150, d = 0.25, sig.level = 0.05,
           power = NULL)


# Etudio de tamaño muestral dado el poder y efecto
# de tamaño 0.25
library(pwr)
pwr.t.test(n = NULL, d = 0.25, sig.level = 0.05,
           power = 0.8)


# install.packages("randomizr")
# install.packages("estimatr")
library(randomizr)
library(estimatr)

## Y0 va fijo 
## Solo lo generamos una vez:
make_Y0 <- function(N){  rnorm(n = N) }
repeat_experiment_and_test <- function(N, Y0, tau){
    Y1 <- Y0 + tau
    Z <- complete_ra(N = N)
    Yobs <- Z * Y1 + (1 - Z) * Y0
    estimator <- lm_robust(Yobs ~ Z)
    pval <- estimator$p.value[2]
    return(pval)
}


power_sim <- function(N,tau,sims){
    Y0 <- make_Y0(N)
    pvals <- replicate(n=sims,
      repeat_experiment_and_test(N=N,Y0=Y0,tau=tau))
    pow <- sum(pvals < .05) / sims
    return(pow)
}

set.seed(12345)
power_sim(N=80,tau=.25,sims=100)
power_sim(N=80,tau=.25,sims=100)

# install.packages("DeclareDesign")
# install.packages("tidyverse")
library(DeclareDesign)
library(tidyverse)


my_design <- declare_model(N = 100,
U = rnorm(N),
Y_Z_0 = U,
Y_Z_1 = U + rnorm(N, mean = 2, sd = 2)) +
declare_assignment(Z = complete_ra(N, prob = 0.5)) +
declare_inquiry(ATE = mean(Y_Z_1 - Y_Z_0)) +
declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
declare_estimator(Y ~ Z, inquiry = "ATE")
diagnose_design(my_design)



my_design2 <- declare_model(N = 100,
U = rnorm(N),
Y_Z_0 = U,
Y_Z_1 = U + rnorm(N, mean = 1, sd = 3)) +
declare_assignment(Z = complete_ra(N, prob = 0.5)) +
declare_inquiry(ATE = mean(Y_Z_1 - Y_Z_0)) +
declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
declare_estimator(Y ~ Z, inquiry = "ATE")
diagnose_design(my_design2)



my_design3 <- declare_model(N = 100,
U = rnorm(N),
Y_Z_0 = U,
Y_Z_1 = U + rnorm(N, mean = 1, sd = 3)) +
declare_assignment(Z = complete_ra(N, prob = 0.5)) +
declare_inquiry(ATE = mean(Y_Z_1 - Y_Z_0)) +
declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
declare_estimator(Y ~ Z, model = lm_robust, inquiry = "ATE")
diagnose_design(my_design3)


## Y0 va fijo
make_Y0_cov <- function(N){
    u0 <- rnorm(n = N)
    x <- rpois(n=N,lambda=2)
    Y0 <- .5 * sd(u0) * x + u0
    return(data.frame(Y0=Y0,x=x))
 }
##  X predice Y0.
test_dat <- make_Y0_cov(100)
test_lm  <- lm_robust(Y0~x,data=test_dat)
summary(test_lm)

## la simulación
repeat_experiment_and_test_cov <- function(N, tau, Y0, x){
    Y1 <- Y0 + tau
    Z <- complete_ra(N = N)
    Yobs <- Z * Y1 + (1 - Z) * Y0
    estimator <- lm_robust(Yobs ~ Z+x,data=data.frame(Y0,Z,x))
    pval <- estimator$p.value[2]
    return(pval)
}
## crear los datos una vez, asignar de manera aleatoria el tratamiento
## reportar la proporción que devuelve p-valor < 0.05
power_sim_cov <- function(N,tau,sims){
    dat <- make_Y0_cov(N)
    pvals <- replicate(n=sims, repeat_experiment_and_test_cov(N=N,
                          tau=tau,Y0=dat$Y0,x=dat$x))
    pow <- sum(pvals < .05) / sims
    return(pow)
}


set.seed(12345)
power_sim_cov(N=80,tau=.25,sims=100)




N = 100

model <- declare_model(
  N = N, 
  X = sort(rnorm(N)), 
  u = rnorm(N), 
  potential_outcomes(Y ~ Z - X * Z + u/10)
)
inquiry     <- declare_inquiry(ATE = mean(Y_Z_1 - Y_Z_0))
blocking    <- declare_step(fabricate, couples = rep(1:(N/2), each = 2))
assignment  <- declare_assignment(Z = block_ra(blocks = couples))
measurement <- declare_measurement(Y = reveal_outcomes(Y ~ Z))
estimator   <- declare_estimator(Y ~ Z)

design_likes <-  model + inquiry +
                 blocking + assignment + measurement + estimator
diagnose_design(design_likes)


## Y0 va fijo
make_Y0_clus <- function(n_indivs,n_clus){
    # n_indivs numero de individuos en un cluster
    # n_clus numero de clusters
    clus_id <- gl(n_clus,n_indivs)
    N <- n_clus * n_indivs
    u0 <- fabricatr::draw_normal_icc(N=N,clusters=clus_id,ICC=.1)
    Y0 <- u0
    return(data.frame(Y0=Y0,clus_id=clus_id))
 }

test_dat <- make_Y0_clus(n_indivs=10,n_clus=100)
# confirmar que esto produce data con 10 en cada uno de los 100 clusters 
table(test_dat$clus_id)
# confirmar el icc
ICC::ICCbare(y=Y0,x=clus_id,data=test_dat)


repeat_experiment_and_test_clus <- function(N, tau, Y0, clus_id){
    Y1 <- Y0 + tau
    # randomizamos z al nivel del cluster
    Z <- cluster_ra(clusters=clus_id)
    Yobs <- Z * Y1 + (1 - Z) * Y0
    estimator <- lm_robust(Yobs ~ Z, clusters = clus_id,
                    data=data.frame(Y0,Z,clus_id), se_type = "CR2")
    pval <- estimator$p.value[2]
    return(pval)
  }
power_sim_clus <- function(n_indivs,n_clus,tau,sims){
    dat <- make_Y0_clus(n_indivs,n_clus)
    N <- n_indivs * n_clus
    # aleatorizar el tratamiento sims veces
    pvals <- replicate(n=sims,
                repeat_experiment_and_test_clus(N=N,tau=tau,
                                Y0=dat$Y0,clus_id=dat$clus_id))
    pow <- sum(pvals < .05) / sims
    return(pow)
}


set.seed(12345)
power_sim_clus(n_indivs=8,n_clus=100,tau=.25,sims=100)
power_sim_clus(n_indivs=8,n_clus=100,tau=.25,sims=100)



some_ns <- seq(10,800,by=10)
pow_by_n <- sapply(some_ns, function(then){
    pwr.t.test(n = then, d = 0.25, sig.level = 0.05)$power
            })
plot(some_ns,pow_by_n,
    xlab="Tamaño muestral",
    ylab="Poder")
abline(h=.8)
## See https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html
## for fancier plots
## ptest <-  pwr.t.test(n = NULL, d = 0.25, sig.level = 0.05, power = .8)
## plot(ptest)


some_taus <- seq(0,1,by=.05)
pow_by_tau <- sapply(some_taus, function(thetau){
    pwr.t.test(n = 200, d = thetau, sig.level = 0.05)$power
            })
plot(some_taus,pow_by_tau,
    xlab="Efecto causal promedio (estandarizado)",
    ylab="Poder")
abline(h=.8)

