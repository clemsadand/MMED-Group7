# Load necessary library
library(deSolve)
library(ggplot2)

# Define the SEIRV model with behavioral change
SEIRV_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Total populations for risk-averse (SA + EA + IA + RA + VA) and risk-tolerant (ST + ET + IT + RT + VT)
    N_A <- SA + EA + IA + RA + VA
    N_T <- ST + ET + IT + RT + VT
    
    # Total infected individuals
    total_infected <- IA + IT
    
    # Adjust beta_AA based on the total infections (IA + IT)
    beta_AA_adj <- beta_AA / (1 + c * total_infected)
    
    # Proportion of population shifting to risk-tolerant (stopping wearing masks) based on infection levels
    shifting_rate <- k * total_infected
    
    # Rate of change for risk-averse group
    dSA <- b * (SA + EA + IA + RA) - beta_AA_adj * SA * IA / N_A - mu * SA - tau * SA - shifting_rate * SA
    dEA <- beta_AA_adj * SA * IA / N_A - gamma * EA - mu * EA
    dIA <- gamma * EA - delta * IA - epsilon * IA - mu * IA
    dRA <- delta * IA - mu * RA
    dVA <- tau * SA - mu * VA
    
    # Rate of change for risk-tolerant group
    dST <- b * (ST + ET + IT + RT) - beta_TA * ST * IA / N_T - beta_TT * ST * IT / N_T - mu * ST - tau * ST + shifting_rate * SA
    dET <- beta_TA * ST * IA / N_T + beta_TT * ST * IT / N_T - gamma * ET - mu * ET
    dIT <- gamma * ET - delta * IT - epsilon * IT - mu * IT
    dRT <- delta * IT - mu * RT
    dVT <- tau * ST - mu * VT
    
    # Return the rate of change
    list(c(dSA, dEA, dIA, dRA, dVA, dST, dET, dIT, dRT, dVT))
  })
}

# Parameters
parameters <- c(
  b = 0.01,      # Birth rate
  beta_AA = 0.5, # Transmission rate within risk-averse
  beta_TA = 0.3, # Transmission rate from risk-tolerant to risk-averse
  beta_TT = 0.4, # Transmission rate within risk-tolerant
  gamma = 0.1,   # Incubation rate
  delta = 0.05,  # Recovery rate
  epsilon = 0.01, # Disease-induced death rate
  mu = 0.02,     # Natural death rate
  tau = 0.01,    # Vaccination rate
  c = 0.001,     # Behavioral change parameter for beta_AA
  k = 0.002     # Rate at which people stop wearing masks as infections decrease
)

# Initial state
initial_state <- c(
  SA = 1000, EA = 10, IA = 5, RA = 0, VA = 0, 
  ST = 1000, ET = 10, IT = 5, RT = 0, VT = 0
)

# Time vector
times <- seq(0, 365, by = 1)  # Simulate for 1 year

# Solve the system of differential equations
output <- ode(y = initial_state, times = times, func = SEIRV_model, parms = parameters)

# Convert output to data frame for easier plotting
output_df <- as.data.frame(output)

# Plot the results
# Plot total infected over time for both groups
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = SA, color = "Risk-averse Susceptible")) +
  geom_line(aes(y = EA, color = "Risk-averse Exposed")) +
  geom_line(aes(y = IA, color = "Risk-averse Infected")) +
  geom_line(aes(y = RA, color = "Risk-averse Recovered")) +
  geom_line(aes(y = ST, color = "Risk-tolerant Susceptible")) +
  geom_line(aes(y = ET, color = "Risk-tolerant Exposed")) +
  geom_line(aes(y = IT, color = "Risk-tolerant Infected")) +
  geom_line(aes(y = RT, color = "Risk-tolerant Recovered")) +
  labs(title = "SEIRV Model with Behavioral Change", x = "Time (days)", y = "Population") 
theme_minimal()

plot(output,times)

