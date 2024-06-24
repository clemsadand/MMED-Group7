library(deSolve)# Load libary to be used for numerical integration
library(ggplot2)

#Model parameters --------------------------------------------------------------

seirv_params <- list(
  a = 50,
  c = 30,
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


# fn.control <- function(t, cMax, cRate, cHalf){
#   min(1,1-cMax/(1+exp(-cRate*(t-cHalf))))
# }

# Model definition -------------------------------------------------------------

seirv <- function(t, y, parms){
  with(c(as.list(y), parms), {
    N_A <- S_A + E_A + I_A + R_A + V_A
    N_T <- S_T + E_T + I_T + R_T + V_T
    N <- N_A + N_T
    I <- I_A + I_T
    mu = b#background death rate
<<<<<<< HEAD
    #Transmission coefficients
=======
    # lambda_AA <- beta_AA*I_A/N_A; lambda_TA <- beta_TA*I_T/N_T;
    # lambda_TT <- beta_TT*I_T/N_T; lambda_AT <- beta_AT*I_A/N_A;
    # lambda_AA <- beta_AA * exp(-a*I_A/N_A); lambda_TA <- beta_TA * exp(-a*I_T/N_T);
    # lambda_TT <- beta_TT * exp(-a*I_T/N_T); lambda_AT <- beta_AT * exp(-a*I_A/N_A);
>>>>>>> f81ceca01a60225043573a95eac3b83bf79b9570
    lambda_AA <- beta_AA * exp(-a*I/N)* I_A/N_A; 
    lambda_TA <- beta_TA * exp(-c*I/N)*I_T/N_T;
    lambda_TT <- beta_TT*I_T/N_T ; 
    lambda_AT <- beta_AT *exp(-c*I/N)*I_A/N_A;
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
  ts$P_A <- with(ts, I_A/N0_A)
  ts$P_T <- with(ts, I_T/N0_T)
  return(ts)
}


# Data simulation --------------------------------------------------------------

## Function to 'sample' the population:
## From a simulated epidemic, measure prevalence at several time points by drawing
## cross-sectional samples of individuals at each time, testing them, and then calculating sample
## prevalence and associated binomial confidence intervals

sampleEpidemic <- function(simDat # Simulated "data" which we treat as real 
                           , sampleDates = seq(0, 60, by = 10) # Sample every 3 days
                           , numSamp = rep(80, length(sampleDates)) # Number of individuals sampled at each time point
){
  # Risk averse
  prev_at_sample_times <- simDat[simDat$time %in% sampleDates, 'P_A']
  numPos_A <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
  lci_A <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos_A, n = numSamp)
  uci_A <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos_A, n = numSamp)    
  # Risk tolerant
  prev_at_sample_times <- simDat[simDat$time %in% sampleDates, 'P_T']
  numPos_T <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
  lci_T <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos_T, n = numSamp)
  uci_T <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos_T, n = numSamp)    
  #
  return(data.frame(time = sampleDates, numPos_A, numPos_T, numSamp,
                    sampPrev_A =  numPos_A/numSamp,
                    sampPrev_T =  numPos_T/numSamp,
                    lci_A = lci_A, uci_A = uci_A,
                    lci_T = lci_T, uci_T = uci_T
                    ))
}

## Set the time
num_years <- 0.5
time.out <- seq(0, num_years*365, by=1)

## Integration
ts_seirv <- simEpidemic(flu_zero, time.out, seirv, seirv_params)

## Simulated data

with(ts_seirv, plot(time, P_T, xlab = '', ylab = 'Prevalence in T', type = 'l', ylim = c(0,.4), col='green', las = 1))
## Take cross-sectional sample of individuals to get prevalence estimates at multiple time points:
set.seed(5) # Initiate the random number generator (so everyone's simulation looks the same)
myDat <- sampleEpidemic(ts_seirv) # Simulate data from the sampling process (function defined above)
points(myDat$time, myDat$sampPrev_T, col = 'green', pch = 20, cex = 2) # Plot sample prevalence at each time point
arrows(myDat$time, myDat$uci_T, myDat$time, myDat$lci_T, col = 'green', len = .025, angle = 90, code = 3) # Plot 95% CIs around the sample prevalences
##
##
##
with(ts_seirv, plot(time, P_A, xlab = '', ylab = 'Prevalence in A', type = 'l', ylim = c(0,.4), col='red', las = 1))
## Take cross-sectional sample of individuals to get prevalence estimates at multiple time points:
set.seed(7) # Initiate the random number generator (so everyone's simulation looks the same)
myDat <- sampleEpidemic(ts_seirv) # Simulate data from the sampling process (function defined above)
points(myDat$time, myDat$sampPrev_A, col = 'red', pch = 20, cex = 2) # Plot sample prevalence at each time point
arrows(myDat$time, myDat$uci_A, myDat$time, myDat$lci_A, col = 'red', len = .025, angle = 90, code = 3) # Plot 95% CIs around the sample prevalences


# Plotting ---------------------------------------------------------------------

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

max_y <- max(c(max(ts_seirv$S_T), max(ts_seirv$S_A)))

plot(ts_seirv$time,               # Time on the x axis
     ts_seirv$S_T,                  # Number infected (I) on the y axis
     xlab = "Time (days)",     # Label the x axis
     ylab = "Number of susceptibles",  # Label the y axis
     main = "Influenza",    # Plot title
     xlim = c(0, max(ts_seirv$time)),           #
     ylim = c(0, max_y * 1.1),
     col="green",
     type = "l",                # Use a line plot
     bty = "n")                 # Remove the box around the plot


lines(ts_seirv$time, ts_seirv$S_A, col="red")

legend("topright", legend = c("Risk tolerant", "Risk averse"), col = c("green", "red"), lwd = 2)

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