library(deSolve)# Load libary to be used for numerical integration

#Initialization of parameters
seirv_params <- list(
  beta.AA = 0.1,
  beta.AT = 0.1,
  beta.TA = 0.1,
  beta.TT = 0.1,
  N.A = 100000,
  N.T = 100000,
  ave_lifespan = 70*365,#years
  ave_influenza_duration = 5,#days, from 3 to 7
  ave_incubation_period = 2,#days, from 1 to 4 
  ave_recovery_rate = 10,#days
  ave_vaccination_period = 8*31,#????September to April
  ave_vaccination_duration = 6*31#months
)

seirv_pre_vax <- within(seirv_params, {
  ave_vaccination_period <- Inf
})

#Initial state
flu_zero <- c(S.A = seirv_params$N.A - 15, E.A = 5, I.A = 8, R.A = 1, V.A = 1,
              S.T = seirv_params$N.T - 2, E.T = 0, I.T = 1, R.T = 1, V.T = 1)



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
        dR.Adt <- - recovery.A
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
        dR.Tdt <- - recovery.T
        dV.Tdt <- vaccination.T - immunity_loss.T

        # Return the rate of change for each compartment adjusted for natural deaths
        return(list(c(dS.Adt,dE.Adt,dI.Adt,dR.Adt,dV.Adt,dS.Tdt,dE.Tdt,dI.Tdt,dR.Tdt,dV.Tdt)-death))
    })
}

#Set the time
num_years <- 3
time.out <- seq(0, num_years*365, 7)

#test if it works 
out <- seirv(time.out[2], flu_zero, seirv_params)
names(out)
# out <- out + seirv(time.out[2], flu_zero, seirv_params)#*0.1

out

f.xbytime <- function(x0, parms, deltat, maxt) {
    # initialise x: state variables
    
    x <- x0
    S.A <- x['S.A']
    E.A <- x['E.A']
    I.A <- x['I.A']
    R.A <- x['R.A']
    V.A <- x['V.A']
    S.T <- x['S.T']
    E.T <- x['E.T']
    I.T <- x['I.T']
    R.T <- x['R.T']
    V.T <- x['V.T']
    # N <- sum(x) # S + I + R
    
    # initialise y: columns for time and state variables 
    
    y <- data.frame(time = 0, S.A = S.A, E.A = E.A, I.A = I.A, R.A = R.A, V.A = V.A,
        S.T = S.T, E.T = E.T, I.T = I.T, R.T = R.T, V.T = V.T
    )
    #
    for (currenttime in seq(from = deltat, to = maxt, by = deltat)) {
        xx <- names(x[[1]])
        x <- xx + seirv(currenttime, x, parms)*deltat
        y <- rbind(y
                , c(time = currenttime, S.A = x[['S.A']], E.A = x[['E.A']], I.A = x[['I.A']], R.A = x[['R.A']], V.A = x[['V.A']],
                S.T = x[['S.T']], E.T = x[['E.T']], I.T = x[['I.T']], R.T = x[['R.T']], V.T = x[['V.T']]
                ))
        
    }
    
    return(y)
}

f.xbytime(flu_zero, seirv_params, deltat=0.1, 1)

ts_seirv <- tail(data.frame(lsoda(
  y = flu_zero,               # Initial conditions for population
  times = time.out,             # Timepoints for evaluation
  func = seirv,                   # Function to evaluate
  parms = seirv_params#seirv_pre_vax                # Vector of parameters
)), 100)

head(ts_seirv)


plot(ts_seirv$time,               # Time on the x axis
     ts_seirv$S.T,                  # Number infected (I) on the y axis
     xlab = "Time in days",     # Label the x axis
     ylab = "Number infected",  # Label the y axis
     main = "Influenza",    # Plot title
     xlim = c(0, num_years*365),           #
     type = "l",                # Use a line plot
     bty = "n")                 # Remove the box around the plot


