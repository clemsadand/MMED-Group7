library(deSolve)# Load libary to be used for numerical integration
library(ggplot2)

#Model parameters --------------------------------------------------------------

seirv_params <- list(
  a = 10,
  c = 0,
  beta_AA = 2,
  beta_AT = 2,
  beta_TA = 2,
  beta_TT = 2,
  gamma = 1/4,#incubation rate 
  tau = 1/365,#loss immunity rate
  delta = 1/10,#recovery rate
  nu = 1/(0.4*365),#vaccination rate
  b = 1/(70*365),#birth rate
  epsilon = 0.1#disease death rate 
)

# Initial state ----------------------------------------------------------------
N0_A = 1000#
N0_T = 1000#
flu_zero <- c(S_A = N0_A-1, E_A = 0, I_A = 1, R_A = 0, V_A = 0,
              S_T = N0_T-1, E_T = 0, I_T = 1, R_T = 0, V_T = 0)


# Model definition -------------------------------------------------------------

seirv <- function(t, y, parms){
  with(c(as.list(y), parms), {
    N_A <- S_A + E_A + I_A + R_A + V_A
    N_T <- S_T + E_T + I_T + R_T + V_T
    N <- N_A + N_T
    I <- I_A + I_T
    mu = b#background death rate
    #Transmission coefficients
    lambda_AA <- beta_AA * exp(-a*I/N)* I_A/N_A; 
    lambda_TA <- beta_TA * exp(-a*I/N)*I_T/N_T;
    lambda_TT <- beta_TT*exp(-c*I/N)*I_T/N_T ; 
    lambda_AT <- beta_AT*exp(-c*I/N)*I_A/N_A;
    #
    # FOI_AA <- lambda_AA * I/N; FOI_TA <- lambda_TA * I_T/N;
    # FOI_TT <- lambda_TT * I/N; FOI_AT <- lambda_AT * I_A/N;
    
    # risk averse
    dS_Adt <- (b*N_A + tau * V_A) - (mu + nu + lambda_AA + lambda_TA)*S_A 
    dE_Adt <- (lambda_AA + lambda_TA) * S_A - (gamma + mu) * E_A
    # dS_Adt <- (b*N_A + tau * V_A) - (mu + nu + FOI_AA + FOI_TA)*S_A
    # dE_Adt <- (FOI_AA + FOI_TA) * S_A - (gamma + mu) * E_A
    dI_Adt <- gamma * E_A - (delta + epsilon + mu) * I_A
    dR_Adt <- delta * I_A - mu * R_A
    dV_Adt <- nu * S_A - (tau + mu) * V_A
    
    #risk tolerant
    dS_Tdt <- (b*N_T + tau * V_T) - (mu + nu + lambda_TT + lambda_AT)*S_T
    dE_Tdt <- (lambda_TT + lambda_AT) * S_T - (gamma + mu) * E_T
    # dS_Tdt <- (b*N_T + tau * V_T) - (mu + nu + FOI_TT + FOI_AT)*S_T
    # dE_Tdt <- (FOI_TT + FOI_AT) * S_T - (gamma + mu) * E_T
    dI_Tdt <- gamma * E_T - (delta + epsilon + mu) * I_T
    dR_Tdt <- delta * I_T - mu * R_T
    dV_Tdt <- nu * S_T - (tau + mu) * V_T
    # return a list 
    return(list(c(dS_Adt, dE_Adt, dI_Adt, dR_Adt, dV_Adt, dS_Tdt, dE_Tdt, dI_Tdt, dR_Tdt, dV_Tdt)))
  })
}


#### SEIRv with inytervention -------------------------------------------------------------------------

fn.control <- function(t, cMax, cRate, cHalf) {
  min(1, 1 - cMax / (1 + exp(-cRate * (t - cHalf))))
}

seirv_intervention <- function(t, y, parms){
  with(c(as.list(y), parms), {
    N_A <- S_A + E_A + I_A + R_A + V_A
    N_T <- S_T + E_T + I_T + R_T + V_T
    N <- N_A + N_T
    I <- I_A + I_T
    mu = b#background death rate
    #Transmission coefficients
    cMax = 1.#
    cRate = 0.2#
    cHalf = 10#
    
    interv <- fn.control(t, cMax, cRate, cHalf)##intervention
    lambda_AA <- beta_AA * exp(-a*I/N) * interv * I_A/N_A; 
    lambda_TA <- beta_TA * exp(-a*I/N) * interv * I_T/N_T;
    lambda_TT <- beta_TT * I_T/N_T * interv ;
    lambda_AT <- beta_AT * I_A/N_A * interv ;
    
    # risk averse
    dS_Adt <- (b*N_A + tau * V_A) - (mu + nu + lambda_AA + lambda_TA)*S_A 
    dE_Adt <- (lambda_AA + lambda_TA) * S_A - (gamma + mu) * E_A
    # dS_Adt <- (b*N_A + tau * V_A) - (mu + nu + FOI_AA + FOI_TA)*S_A
    # dE_Adt <- (FOI_AA + FOI_TA) * S_A - (gamma + mu) * E_A
    dI_Adt <- gamma * E_A - (delta + epsilon + mu) * I_A
    dR_Adt <- delta * I_A - mu * R_A
    dV_Adt <- nu * S_A - (tau + mu) * V_A
    
    #risk tolerant
    dS_Tdt <- (b*N_T + tau * V_T) - (mu + nu + lambda_TT + lambda_AT)*S_T
    dE_Tdt <- (lambda_TT + lambda_AT) * S_T - (gamma + mu) * E_T
    # dS_Tdt <- (b*N_T + tau * V_T) - (mu + nu + FOI_TT + FOI_AT)*S_T
    # dE_Tdt <- (FOI_TT + FOI_AT) * S_T - (gamma + mu) * E_T
    dI_Tdt <- gamma * E_T - (delta + epsilon + mu) * I_T
    dR_Tdt <- delta * I_T - mu * R_T
    dV_Tdt <- nu * S_T - (tau + mu) * V_T
    # return a list 
    return(list(c(dS_Adt, dE_Adt, dI_Adt, dR_Adt, dV_Adt, dS_Tdt, dE_Tdt, dI_Tdt, dR_Tdt, dV_Tdt)))
  })
}

# Test the function seirv() ---------------------------------------------------
# seirv(0.1, flu_zero, seirv_params)


# Integration of the ODEs ------------------------------------------------------

## This function compute the ODEs' solution
simEpidemic <- function(initial_state, timestep, model, params){
  ## Integration with
  ts <- data.frame(lsoda(
    y = initial_state,# Initial conditions for population
    times = timestep, # Timepoints for evaluation
    func = model,# Function to evaluate
    parms = params# Vector of parameters
  ))
  N = 
  ts$P <- with(ts, (I_A+I_T)/(N0_A+N0_T))
  # ts$P_T <- with(ts, I_T/N0_T)
  return(ts)
}


# Data simulation --------------------------------------------------------------

## Function to 'sample' the population:
## From a simulated epidemic, measure prevalence at several time points by drawing
## cross-sectional samples of individuals at each time, testing them, and then calculating sample
## prevalence and associated binomial confidence intervals

# sampleEpidemic <- function(simDat # Simulated "data" which we treat as real 
#                            , sampleDates = seq(0, 60, by = 10) # Sample every 3 days
#                            , numSamp = rep(80, length(sampleDates)) # Number of individuals sampled at each time point
# ){
#   # Risk averse
#   prev_at_sample_times <- simDat[simDat$time %in% sampleDates, 'P_A']
#   numPos_A <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
#   lci_A <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos_A, n = numSamp)
#   uci_A <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos_A, n = numSamp)    
#   # Risk tolerant
#   prev_at_sample_times <- simDat[simDat$time %in% sampleDates, 'P_T']
#   numPos_T <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
#   lci_T <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos_T, n = numSamp)
#   uci_T <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos_T, n = numSamp)    
#   #
#   return(data.frame(time = sampleDates, numPos_A, numPos_T, numSamp,
#                     sampPrev_A =  numPos_A/numSamp,
#                     sampPrev_T =  numPos_T/numSamp,
#                     lci_A = lci_A, uci_A = uci_A,
#                     lci_T = lci_T, uci_T = uci_T
#                     ))
# }

sampleEpidemic <- function(simDat # Simulated "data" which we treat as real 
                           , sampleDates = seq(0,100, by = 3) # Sample every 3 days
                           , numSamp = rep(80, length(sampleDates)) # Number of individuals sampled at each time point
){
  prev_at_sample_times <- simDat[simDat$time %in% sampleDates, 'P']
  numPos <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
  lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos, n = numSamp)
  uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos, n = numSamp)    
  return(data.frame(time = sampleDates, numPos, numSamp, sampPrev =  numPos/numSamp,
                    lci = lci, uci = uci))
}

## Set the time
num_years <- 0.5
time.out <- seq(0, num_years*365, by=1)

## Integration
####
# Taks: I want to plot for differenr value of the t_interv
######
par(mfrow = c(2,2), oma = c(0,0,2,0))
for(sdVal in c(.05, .1, .5, 2)) {
  ts_seirv <- simEpidemic(flu_zero, time.out, seirv_intervention, seirv_params, t_intervention)
  plot(posteriorSample, xlab='iteration', ylab = 'prevalence',
       main = bquote(sigma==.(sdVal)),
       ylim = c(0,1),
       type = 'l', las = 1)
}
mtext('posterior sample by proposer sd', side = 3, line = 0, outer = T)


# Plotting ---------------------------------------------------------------------

#### Infected-----------------------------------------------
max_y <- max(c(max(ts_seirv$I_T), max(ts_seirv$I_A)))
plot(ts_seirv$time,               # Time on the x axis
     ts_seirv$I_T,                  # Number infected (I) on the y axis
     xlab = "Time (days)",     # Label the x axis
     ylab = "Number of infected",  # Label the y axis
     main = "Influenza",    # Plot title
     xlim = c(0, max(ts_seirv$time)),           #
     ylim = c(0, max_y * 1.1),
     col="green",
     type = "l",                # Use a line plot
     bty = "n")                 # Remove the box around the plot
lines(ts_seirv$time, ts_seirv$I_A, col="red")
legend("topright", legend = c("Risk tolerant", "Risk averse"), col = c("green", "red"), lwd = 2)

#### Prevalence ---------------------
# max_y <- max(c(max(ts_seirv$I_T), max(ts_seirv$P)))
plot(ts_seirv$time,               # Time on the x axis
     ts_seirv$P,                  # Number infected (I) on the y axis
     xlab = "Time (days)",     # Label the x axis
     ylab = "Prevalence",  # Label the y axis
     main = "Influenza",    # Plot title
     xlim = c(0, max(ts_seirv$time)),           #
     ylim = c(0, max(ts_seirv$P) * 1.1),
     col="green",
     type = "l",                # Use a line plot
     bty = "n")                 # Remove the box around the plot

####Susceptibles -------------------------------------------------
# max_y <- max(c(max(ts_seirv$S_T), max(ts_seirv$S_A)))
# plot(ts_seirv$time,               # Time on the x axis
#      ts_seirv$S_T,                  # Number infected (I) on the y axis
#      xlab = "Time (days)",     # Label the x axis
#      ylab = "Number of susceptibles",  # Label the y axis
#      main = "Influenza",    # Plot title
#      xlim = c(0, max(ts_seirv$time)),           #
#      ylim = c(0, max_y * 1.1),
#      col="green",
#      type = "l",                # Use a line plot
#      bty = "n")                 # Remove the box around the plot
# lines(ts_seirv$time, ts_seirv$S_A, col="red")
# legend("topright", legend = c("Risk tolerant", "Risk averse"), col = c("green", "red"), lwd = 2)

## Simulated data ----------------------------------------------------

# with(ts_seirv, plot(time, P, xlab = '', ylab = 'Prevalence', type = 'l', ylim = c(0,max(P) * 1.1), col='green', las = 1))
# ## Take cross-sectional sample of individuals to get prevalence estimates at multiple time points:
# # set.seed(0) # Initiate the random number generator (so everyone's simulation looks the same)
# myDat <- sampleEpidemic(ts_seirv) # Simulate data from the sampling process (function defined above)
# points(myDat$time, myDat$sampPrev, col = 'green', pch = 20, cex = 2) # Plot sample prevalence at each time point
# arrows(myDat$time, myDat$uci, myDat$time, myDat$lci, col = 'green', len = .025, angle = 90, code = 3) # Plot 95% CIs around the sample prevalences
##
##
##
# with(ts_seirv, plot(time, P_A, xlab = '', ylab = 'Prevalence in A', type = 'l', ylim = c(0,.4), col='red', las = 1))
# ## Take cross-sectional sample of individuals to get prevalence estimates at multiple time points:
# set.seed(7) # Initiate the random number generator (so everyone's simulation looks the same)
# myDat <- sampleEpidemic(ts_seirv) # Simulate data from the sampling process (function defined above)
# points(myDat$time, myDat$sampPrev_A, col = 'red', pch = 20, cex = 2) # Plot sample prevalence at each time point
# arrows(myDat$time, myDat$uci_A, myDat$time, myDat$lci_A, col = 'red', len = .025, angle = 90, code = 3) # Plot 95% CIs around the sample prevalences



# # Plotting the transmission coefficient
# ts_seirv$prob_inf_A <- with(ts_seirv, exp(-seirv_params$a * P))
# with(ts_seirv, plot(P, prob_inf_A, type="l", xlab="Prevalence", ylab="Transmission coeff in Risk averse", col="red"))


##
## Plot the results
# out_long <- tidyr::pivot_longer(ts_seirv, cols = -time, names_to = "compartment", values_to = "value")
# 
# ggplot(out_long, aes(x = time, y = value, color = compartment)) +
#   geom_line() +
#   theme_bw()+
#   labs(title = "SEIR Model with Vaccination for Two Populations",
#        x = "Time (days)", y = "Number of individuals") +
#   theme_minimal()
# 
# head(ts_seirv, 20)


## Negative log likelihood

# ## Return the negative log of likelihood by convention
# nllikelihood <- function(params, initial_state, obsDat=myDat) {#simEpidemic <- function(initial_state, timestep, model, params)
#   simDat <- simEpidemic(initial_state, obsDat$time, seirv, params)
#   ## What are the rows from our simulation at which we have observed data?
#   matchedTimes <- simDat$time %in% obsDat$time
#   nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P[matchedTimes], log = T)
#   return(sum(nlls))
# }
# 
# nllikelihood(seirv_params, flu_zero, myDat)
# 
# # Initial guess for the parameters
# initial_guess <- list(
#   a = 8,
#   c = 1,
#   beta_AA = 1.2,
#   beta_AT = 2.6,
#   beta_TA = 1.1,
#   beta_TT = 1.0,
#   gamma = 0.1,
#   tau = 0.023,
#   delta = 1/10,
#   nu = 0.02,
#   b = 0.001,
#   epsilon = 0.2
# )

# Observed data >>>> myDat

# # Perform the optimization
# optim_result <- optim(
#   par = seirv_params,
#   fn = nllikelihood,
#   initial_state = flu_zero,
#   obsDat = myDat,
#   control = list(trace = TRUE, maxit = 800, reltol = 10^-7),
#   method = "Nelder-Mead",
#   hessian = TRUE
# )
# 
# optim_result
# MLEfits <- optim_result$par
# exp(unname(MLEfits))
