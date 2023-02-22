library(deSolve)
library(reshape2)
library(ggplot2)

initial_susceptible <- 1000000
initial_exposed <- 1100
initial_infected <- 500
initial_quarantined <- 500
initial_recovered <- 10000
initial_dead <- 100
inital_vaccinated <- 0

initial_state_values <- c(S=initial_susceptible, E=initial_exposed, 
                          I=initial_infected, Q=initial_quarantined, 
                          R=initial_recovered, D=initial_dead, 
                          V=inital_vaccinated)

parameters <- c(LAMBDAA=2300, beta=0.000000858, alpha=0.005, mu=0.00003, 
                gamma=0.1818, sigma=0.05, 
                delta=0.2631578947,kappa=0.014, 
                lambda=0.1,rho=0.06666667)

times <-seq(from = 0, to = 100, by = 1)

seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    N <- S+E+I+Q+R+D+V
    
    dS <- LAMBDAA - beta*S*I - alpha*S - mu*S
    dE <- beta*S*I - gamma*E + sigma*beta*V*I - mu*E
    dI <- gamma*E - delta*I - mu*I
    dQ <- delta*I - (1-kappa)*lambda*Q - kappa*rho*Q - mu*Q
    dR <- (1-kappa)*lambda*Q - mu*R
    dD <- kappa*rho*Q
    dV <- alpha*S - sigma*beta*V*I - mu*V
    
    return(list(c(dS,dE,dI,dQ,dR,dD,dV)))
  })
}

output <- as.data.frame(ode(y = initial_state_values,
                            times = times,
                            func = seir_model,
                            parms = parameters))

output_long <- melt(as.data.frame(output), id = "time")

ggplot(data = output_long, 
       aes(x = time, y = value, color = variable, group = variable)) + 
  geom_line() + 
  xlab("Time(days)") + 
  ylab("Number of people") + 
  labs(title = paste("SEIQRDV Model"), colour = "Compartment")