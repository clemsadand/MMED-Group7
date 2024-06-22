library(deSolve)# Load libary to be used for numerical integration

seirv_params <- list(
  a = 0.5,
  beta_AA = 0.1,
  beta_AT = 0.1,
  beta_TA = 0.51,
  beta_TT = 0.51,
  N_A = 100000,
  N_T = 100000,
  gamma = 1/4,#incubation rate 
  tau = 1/365,#loss immunity rate
  delta = 1/10,#recovery rate
  nu = 1/(0.4*365),#vaccination rate
  b = 1/(70*365),#birth rate
  epsilon = 0.1
)

#Initial state
flu_zero <- c(S_A = seirv_params$N_A - 1, E_A = 0, I_A = 1, R_A = 0, V_A = 0,
              S_T = seirv_params$N_T - 1, E_T = 0, I_T = 1, R_T = 0, V_T = 0)


fn.control <- function(t, cMax, cRate, cHalf){
  min(1,1-cMax/(1+exp(-cRate*(t-cHalf))))
}

seirv <- function(t, y, parms){
  with(c(as.list(y), parms), {
    mu = b#background death rate
    # lambda_AA <- beta_AA * exp(-a*I_A/N_A); lambda_TA <- beta_TA * exp(-a*I_T/N_T); 
    # lambda_TT <- beta_TT * exp(-a*I_T/N_T); lambda_AT <- beta_AT * exp(-a*I_T/N_T);
    lambda_AA <- beta_AA; lambda_TA <- beta_TA; 
    lambda_TT <- beta_TT; lambda_AT <- beta_AT;
    #
    # FOI_AA <- lambda_AA * fn.control(t,1.0, 0.2,0.5) * (I_A+I_T)/(N_A+N_T); FOI_TA <- lambda_TA * fn.control(t,1.0, 0.2,0.5) * (I_A+I_T)/(N_A+N_T); 
    # FOI_TT <- lambda_TT * fn.control(t,1.0, 0.2,0.5) * (I_A+I_T)/(N_A+N_T); FOI_AT <- lambda_AT * fn.control(t,1.0, 0.2,0.5) * (I_A+I_T)/(N_A+N_T);
    #risk averse
    dS_Adt <- (b*N_A + tau * V_A) - (mu + nu + lambda_AA + lambda_TA)*S_A
    dE_Adt <- (lambda_AA + lambda_TA) * S_A - (gamma + mu) * E_A
    # dS_Adt <- (b*N_A + tau * V_A) - (mu + nu + FOI_AA + FOI_TA)*S_A
    # dE_Adt <- (FOI_AA + FOI_TA) * S_A - (gamma + mu) * E_A
    dI_Adt <- gamma * E_A - (delta + beta_AA + beta_AT + epsilon + mu) * I_A
    dR_Adt <- delta * I_A - mu * R_A
    dV_Adt <- nu * S_A - (tau + mu) * V_A
    #risk tolerant
    dS_Tdt <- (b*N_T + tau * V_T) - (mu + nu + lambda_TT + lambda_AT)*S_T
    dE_Tdt <- (lambda_TT + lambda_AT) * S_T - (gamma + mu) * E_T
    # dS_Tdt <- (b*N_T + tau * V_T) - (mu + nu + FOI_TT + FOI_AT)*S_T
    # dE_Tdt <- (FOI_TT + FOI_AT) * S_T - (gamma + mu) * E_T
    dI_Tdt <- gamma * E_T - (delta + beta_TT + beta_TA + epsilon + mu) * I_T
    dR_Tdt <- delta * I_T - mu * R_T
    dV_Tdt <- nu * S_T - (tau + mu) * V_T
    return(list(c(dS_Adt, dE_Adt, dI_Adt, dR_Adt, dV_Adt, dS_Tdt, dE_Tdt, dI_Tdt, dR_Tdt, dV_Tdt)))
  })
}

seirv(0.1, flu_zero, seirv_params)

#Set the time
num_years <- 2
time.out <- seq(0, num_years*365, by=0.1)

parameters <- c(b = 0.01, nu = 0.02, tau = 0.03, gamma = 0.1, delta = 0.1, epsilon = 0.01,
                beta_AA = 0.5, beta_TA = 0.3, beta_TT = 0.4, beta_AT = 0.2, a = 0.05,
                N_A = 2000, N_T = 2000)

ts_seirv <- data.frame(lsoda(
  y = flu_zero,               # Initial conditions for population
  times = time.out,             # Timepoints for evaluation
  func = seirv,                   # Function to evaluate
  parms = seirv_params          # Vector of parameters
))

head(ts_seirv)

# max_y <- max(c(max(ts_seirv$I_T), max(ts_seirv$I_A)))


plot(ts_seirv$time,               # Time on the x axis
     ts_seirv$I_T,                  # Number infected (I) on the y axis
     xlab = "Tim0e in weeks",     # Label the x axis
     ylab = "Number of infected",  # Label the y axis
     main = "Influenza",    # Plot title
     xlim = c(0, num_years*365),           #
     # ylim = c(0, max_y * 1.1),
     col="blue",
     type = "l",                # Use a line plot
     bty = "n")                 # Remove the box around the plot


lines(ts_seirv$time, ts_seirv$I_A, col="red")

legend("topright", legend = c("Risk tolerant", "Risk averse"), col = c("blue", "red"), lwd = 2)
