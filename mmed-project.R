library(deSolve)# Load libary to be used for numerical integration

#Initialization of parameters
seirv_params <- list(
  beta.AA = 0.1,
  beta.AT = 0.01,
  beta.TA = 0.001,
  beta.TT = 0.1,
  N.A = 100000,
  N.T = 100000,
  ave_lifespan = 70*365,#years
  ave_influenza_duration = 4,#days, from 3 to 7
  ave_incubation_period = 4,#days, from 1 to 4 
  ave_recovery_rate = 10,#days
  ave_vaccination_period = 365*0.8,#????September to April
  ave_vaccination_duration = 365#months
)

seirv_pre_vax <- within(seirv_params, {
  ave_vaccination_period <- Inf
})

#Initial state
flu_zero <- c(S.A = seirv_params$N.A - 14, E.A = 5, I.A = 8, R.A = 0, V.A = 1,
              S.T = seirv_params$N.T - 2, E.T = 0, I.T = 1, R.T = 0, V.T = 1)



seirv <- function(t, y, params){
    with(c(as.list(y), params),{
        death <- y * (1/ave_lifespan)
        disease_death <- I.A * (1/ave_influenza_duration)
        birth <- sum(death) + disease_death
        
        ##risk averse
        ##processes relative to the group
        infection.A <- beta.AA*I.A*S.A/N.A
        contact_suscT_to_suscA <- beta.TA*I.T*S.A/N.T
        incubation.A <- 1 / (ave_incubation_period) * E.A
        recovery.A <- 1 / (ave_recovery_rate) * I.A
        vaccination.A <- 1 / (ave_vaccination_period) * S.A
        immunity_loss.A <- 1 / (ave_vaccination_duration) * V.A
        #change in time
        dS.Adt <- (birth + immunity_loss.A) - (infection.A + contact_suscT_to_suscA + vaccination.A)
        dE.Adt <- infection.A - incubation.A
        dI.Adt <- incubation.A - (recovery.A + disease_death)
        dR.Adt <- - (recovery.A + disease_death)
        dV.Adt <- vaccination.A - immunity_loss.A

        #Risk tolerant
        ##processes relative to the group
        infection.T <- beta.TT*I.T*S.T/N.T
        contact_suscA_to_suscT <- beta.AT*I.A*S.T/N.A
        incubation.T <- 1 / (ave_incubation_period) * E.T
        recovery.T <- 1 / (ave_recovery_rate) * I.T
        vaccination.T <- 1 / (ave_vaccination_period) * S.T
        immunity_loss.T <- 1 / (ave_vaccination_duration) * V.T
        #change in time
        dS.Tdt <- (birth + immunity_loss.T) - (infection.T + contact_suscA_to_suscT + vaccination.T)
        dE.Tdt <- infection.T - incubation.T
        dI.Tdt <- incubation.T - (recovery.T + disease_death)
        dR.Tdt <- - (recovery.T + disease_death)
        dV.Tdt <- vaccination.T - immunity_loss.T

        # Return the rate of change for each compartment adjusted for natural deaths
        return(list(c(dS.Adt,dE.Adt,dI.Adt,dR.Adt,dV.Adt,dS.Tdt,dE.Tdt,dI.Tdt,dR.Tdt,dV.Tdt)-death))
    })
}

#Set the time
num_years <- 4
time.out <- seq(0, num_years*365, 7)

ts_seirv <- tail(data.frame(lsoda(
  y = flu_zero,               # Initial conditions for population
  times = time.out,             # Timepoints for evaluation
  func = seirv,                   # Function to evaluate
  parms = seirv_params              # Vector of parameters
)), 100)

head(ts_seirv,20)


plot(ts_seirv$time,               # Time on the x axis
     ts_seirv$I.T,                  # Number infected (I) on the y axis
     xlab = "Tim0e in weeks",     # Label the x axis
     ylab = "Number of infected",  # Label the y axis
     main = "Influenza",    # Plot title
     xlim = c(0, num_years*365),           #
     col="blue",
     type = "l",                # Use a line plot
     bty = "n")                 # Remove the box around the plot


lines(ts_seirv$time, ts_seirv$I.A, col="red")

legend("topright", legend = c("Risk tolerant", "Risk averse"), col = c("blue", "red"), lwd = 2)