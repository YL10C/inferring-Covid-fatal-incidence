## Xiaotian Xing s2661110;

# nolint start

# Set the working directory to the path of the local folder
setwd("D:/OneDrive/文档/Yilin/Edinburgh/Statistical Programming/Assignment2")
# setwd("D:/PG-CAM/Statistical Programming")
# Read the data file, which contains the column names in the header
data <- read.table("engcov.txt", header = TRUE)

# Daily deaths and julian for the previous 150 days
deaths <- data$nhs[1:150]
julian_days <- data$julian[1:150]


# Used to generate a distribution of duration from infection to death
# meanlog and sdlog are the parameters of a lognormal distribution
# max_duration: the maximum number of days in the distribution
generate_infection_to_death_dist <- function(meanlog, sdlog, max_duration) {
  probabilities <- dlnorm(1:max_duration, meanlog, sdlog)
  probabilities / sum(probabilities)
}


# Used to generate the initial infection time for each deceased case
# deaths: the vector of the number of deaths per day
# julian_days: the corresponding calendar day
# infection_to_death_probs: the duration distribution from infection to death
generate_initial_infection_times <- function(deaths, julian_days, infection_to_death_probs) {
  n <- sum(deaths)
  
  # Use the rep function to generate a death date vector of length n
  death_dates <- rep(julian_days, deaths)
  
  # print(head(death_dates, 10))
  
  # n days are sampled representing the number of days from infection to death
  infection_to_death_days <- sample(1:length(infection_to_death_probs), n, replace = TRUE, prob = infection_to_death_probs)
  
  # Calculate the date of infection
  infection_dates <- death_dates - infection_to_death_days
  
  return(infection_dates)
}

# Assume the maximum possible date is 310 days
max_day <- 310

n <- sum(deaths)

# Extend deaths vector
start_day <- min(julian_days)
deaths_extended <- c(
  rep(0, start_day - 1),
  deaths,
  rep(0, max_day - length(deaths) - start_day + 1)
)

# Define the parameters of a lognormal distribution
meanlog <- 3.152
sdlog <- 0.451
max_duration <- 80

# Generates a distribution of duration from infection to death
infection_to_death_probs <- generate_infection_to_death_dist(meanlog, sdlog, max_duration)

# initialize durations
durations <- sample(1:max_duration, size = n, replace = TRUE, prob = infection_to_death_probs)

# Fitness function: Calculates the modified Pearson statistic
calculate_P <- function(actual_deaths, simulated_deaths, max_day) {
  # Extend the deaths vector to length max_day, filling it with zeros
  # # actual_deaths <- c(actual_deaths, rep(0, max_day - length(actual_deaths)))
  P <- sum((actual_deaths - simulated_deaths)^2 / pmax(1, simulated_deaths))
  return(P)
}


# Define optimization function
optimize_t0 <- function(t0, deaths, durations, n.rep, bs = FALSE, max_day, infection_to_death_probs, max_duration) {
  
  
  # Initializes the variable used to store the result
  P_history <- numeric(n.rep)
  inft <- matrix(0, nrow = max_day, ncol = n.rep)
  
  for (rep in 1:n.rep) {
    # 如果 bs = TRUE，在每次迭代中模拟新的死亡数据
    if (bs) {
      deaths_simulated <- rpois(length(deaths_extended), lambda = deaths_extended)
    } else {
      deaths_simulated <- deaths_extended
    }
    
    # Proceed with the rest of the code
    simulated_death_times <- t0 + durations
    simulated_deaths <- tabulate(simulated_death_times, nbins = max_day)
    
    # Calculate the fitness P using simulated deaths
    P <- calculate_P(actual_deaths = deaths_simulated, simulated_deaths = simulated_deaths, max_day)
    
    # Store the P-value and the number of new infections per day
    P_history[rep] <- P
    inft[, rep] <- tabulate(t0, nbins = max_day)
    
    # Randomly shuffle the index order
    indices <- sample(1:n, n)
    durations_index <- 0

    for (i in indices) {
      # Adjust step size selection
      if (rep < n.rep / 2) {
        steps_set <- c(-8, -4, -2, -1, 1, 2, 4, 8)
      } else if (rep < n.rep * 3 / 4) {
        steps_set <- c(-4, -2, -1, 1, 2, 4)
      } else {
        steps_set <- c(-2, -1, 1, 2)
      }
      
      step <- sample(steps_set, size = 1)
      proposed_t0_i <- t0[i] + step
      
      durations_index <- durations_index + 1

      # Calculate the current and proposed date of death
      current_death_time <- t0[i] + durations[durations_index]
      proposed_death_time <- proposed_t0_i + durations[durations_index]
      
      # Make sure the date of death is legal
      if (current_death_time < start_day || current_death_time > max_day) next
      if (proposed_death_time < start_day || proposed_death_time > max_day) next
      
      # Get the actual death toll from simulated deaths
      di_current <- deaths_simulated[current_death_time]
      di_proposed <- deaths_simulated[proposed_death_time]
      
      # Get a simulated death toll
      dsi_current <- simulated_deaths[current_death_time]
      dsi_proposed <- simulated_deaths[proposed_death_time]
      
      # Calculate the adjusted simulated death toll
      dsi_new_current <- dsi_current - 1
      dsi_new_proposed <- dsi_proposed + 1
      
      # Make sure the simulated death toll is not negative
      if (dsi_new_current < 0) next
      
      # Calculate the P-value contribution before and after adjustment
      P_before_current <- (di_current - dsi_current)^2 / pmax(1, dsi_current)
      P_before_proposed <- (di_proposed - dsi_proposed)^2 / pmax(1, dsi_proposed)
      P_before <- P_before_current + P_before_proposed
      
      P_after_current <- (di_current - dsi_new_current)^2 / pmax(1, dsi_new_current)
      P_after_proposed <- (di_proposed - dsi_new_proposed)^2 / pmax(1, dsi_new_proposed)
      P_after <- P_after_current + P_after_proposed
      
      # Calculate ΔP
      delta_P <- P_after - P_before
      
      # Decide whether to accept the offer
      if (delta_P < 0) {
        # Accept proposal to update t0[i]
        t0[i] <- proposed_t0_i
        
        # Updated simulated death toll
        simulated_deaths[current_death_time] <- dsi_new_current
        simulated_deaths[proposed_death_time] <- dsi_new_proposed
        
        # Update the total p-value
        P <- P + delta_P
      }
    }
  }
  
  return(list(P = P_history, inft = inft, t0 = t0))
}


# t: actual date of death vector
# deaths: the number of deaths per day
# n.rep: indicates the total number of iterations
# bs: logical value indicating whether to bootstrapping
# t0: initial infection time vector
# return list(P, inft, t0)
deconv <- function(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL) {
  
  if (is.null(t0)) {
    t0 <- generate_initial_infection_times(deaths, t, infection_to_death_probs)
  }
  
  
  # Call optimization function
  result <- optimize_t0(
    t0 = t0,
    deaths = deaths,
    durations = durations,
    n.rep = n.rep,
    bs = bs,
    max_day = 310,
    infection_to_death_probs = infection_to_death_probs,
    max_duration = max_duration
  )
  
  return(result)
}


results <- deconv(t = julian_days, deaths = deaths, n.rep = 100, bs = FALSE)
t0_converged <- results$t0
# print(results$P)

# Bootstrap with t0 = t0_converged
results <- deconv(t = julian_days, deaths = deaths, n.rep = 100, bs = TRUE, t0 = t0_converged)
# print(results$P)


# Calculate the predicted date of death
Predict_death_data <- results$t0 + durations 
# print(Predict_death_data)

# Estimated daily death toll
Predict_death_number <- tabulate(Predict_death_data, nbins = max_day)


# Calculate the average number of new infections per day
average_infections <- rowMeans(results$inft)

# Calculated confidence interval
ci_upper <- apply(results$inft, 1, quantile, probs = 0.975)
ci_lower <- apply(results$inft, 1, quantile, probs = 0.025)


# Create the X-axis vector
Day <- seq(1, max_day)

# Create a drawing data frame
plot_data <- data.frame(
  Day = Day,
  EstimatedInfections = average_infections,
  CI_upper = ci_upper,
  CI_lower = ci_lower,
  ExtendedDeaths = Predict_death_number,
  ActualDeaths = deaths_extended
)

library(ggplot2)

lockdown_day <- 84

ggplot(plot_data, aes(x = Day)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = EstimatedInfections, color = "Estimated Infections")) +
  geom_line(aes(y = ActualDeaths, color = "Actual Deaths")) +
  geom_line(aes(y = ExtendedDeaths, color = "Predicted Deaths"), linetype = "solid") +
  geom_vline(xintercept = lockdown_day, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("Estimated Infections" = "blue", "Actual Deaths" = "red", "Extended Predicted Deaths" = "black")) +
  labs(title = "Estimated Infection Trajectory vs Actual Deaths",
       subtitle = "Uncertainty and lockdown illustrated, UK lockdown on Day 84",
       x = "Day",
       y = "Number of Cases",
       color = "Legend") +
  annotate("text", x = Inf, y = Inf, label = "Blue area represents confidence interval", hjust = 1.1, vjust = 1.1, color = "blue", size = 2.5) +
  theme_minimal()

# increase the width of the output image
# ggsave("output.png", width = 9.6, height = 6)



# nolint end