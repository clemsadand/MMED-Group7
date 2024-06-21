library(deSolve)                # Load libary to be used for numerical integration

#Initialization of parameters
seivr_params <- list(
  beta.AA = 0.5/3,
  beta.AT = 0.5/3,
  beta.TA = 0.38/3,
  beta.TT = 0.38/3,
  N.A = 100000,
  N.T = 100000,
  ave_lifespan = 70*365,#years
  ave_influenza_duration = 5,#days, from 3 to 7
  ave_incubation_period = 2,#days, from 1 to 4 
  ave_recovery_rate = 10,#days
  ave_vaccination_period = 8,#????September to April
  ave_vaccination_duration = 6*31#months
)

seiv_pre_vax <- within(seivr_params, {
  ave_vaccination_period <- Inf
})

#Initial state
flu_zero <- c(S.A = seivr_params$N.A - 10, E.A = 0, I.A = 10, V.A = 1, R.A = 0,
              S.T = seivr_params$N.T - 21, E.T = 0, I.T = 1, V.T = 20, R.T = 0)

seivr <- function(t, y, params){
  with(c(as.list(y), params),{
    death <- y * (1/ave_lifespan)
    disease_death <- I.A * (1/ave_influenza_duration)
    birth <- sum(death) + disease_death
    
    ##risk averse
    #processes relative to the group
    infection.A <- beta.AA * S.A * I.A / N.A
    contact_susc_to_infectious.A <- beta.AA*I.A
    contact_suscT_to_suscA <- beta.TA*I.T
    incubation.A <- 1 / (ave_incubation_period) * E.A
    recovery.A <- 1 / (ave_recovery_rate) * I.A
    vaccination.A <- 1 / (ave_vaccination_period) * S.A
    immunity_loss.A <- 1 / (ave_vaccination_duration) * V.A
    #change in time
    dS.Adt <- (birth+immunity_loss.A)-(infection.A+vaccination.A+contact_susc_to_infectious.A)#+contact_suscT_to_suscA)
    dE.Adt <- infection.A-incubation.A
    dI.Adt <- incubation.A-recovery.A-disease_death
    dR.Adt <- recovery.A
    dV.Adt <- vaccination.A-immunity_loss.A
    
    ##risk tolerant
    infection.T <- beta.TT * S.T * I.T / N.T
    contact_susc_to_infectious.T <- beta.TT*I.T
    contact_suscA_to_suscT <- beta.AT*I.A
    incubation.T <- 1 / (ave_incubation_period) * E.T
    recovery.T <- 1 / (ave_recovery_rate) * I.T
    vaccination.T <- 1 / (ave_vaccination_period) * S.T
    immunity_loss.T <- 1 / (ave_vaccination_duration) * V.T
    #change in time
    dS.Tdt <- (birth+immunity_loss.T)-(infection.T+vaccination.T+contact_susc_to_infectious.T)#+contact_suscA_to_suscT)
    dE.Tdt <- infection.T-incubation.T
    dI.Tdt <- incubation.T-recovery.T-disease_death
    dR.Tdt <- recovery.T
    dV.Tdt <- vaccination.T-immunity_loss.T
    return(list(c(dS.Adt,dE.Adt,dI.Adt,dR.Adt,dV.Adt,dS.Tdt,dE.Tdt,dI.Tdt,dR.Tdt,dV.Tdt)-death))
  })
}


time.out <- seq(0, 2*365, 0.1)
# 
# ts_seivr <- tail(data.frame(lsoda(
#   y = flu_zero,               # Initial conditions for population
#   times = time.out,             # Timepoints for evaluation
#   func = seivr,                   # Function to evaluate
#   parms = seivr_params                # Vector of parameters
# )), 100)

ts_seivr <- tail(data.frame(lsoda(
  y = flu_zero,               # Initial conditions for population
  times = time.out,             # Timepoints for evaluation
  func = seivr,                   # Function to evaluate
  parms = seiv_pre_vax                # Vector of parameters
)), 100)

ts_seivr$I.T


plot(ts_seivr$time,               # Time on the x axis
     ts_seivr$I.A,                  # Number infected (I) on the y axis
     xlab = "Time in days",     # Label the x axis
     ylab = "Number infected",  # Label the y axis
     main = "Influenza",    # Plot title
     xlim = c(0, 2*365),           #
     type = "l",                # Use a line plot
     bty = "n")                 # Remove the box around the plot

