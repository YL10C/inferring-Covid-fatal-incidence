## Xiaotian Xing s2661110;

# nolint start

# Set the working directory to the path of the local folder
setwd("D:/OneDrive/文档/Yilin/Edinburgh/Statistical Programming/Assignment2")

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


# Fitness function: Calculates the modified Pearson statistic
calculate_P <- function(actual_deaths, simulated_deaths) {
  # Extend the deaths vector to length max_day, filling it with zeros
  actual_deaths <- c(actual_deaths, rep(0, max_day - length(actual_deaths)))
  P <- sum((actual_deaths - simulated_deaths)^2 / pmax(1, simulated_deaths))
  return(P)
}


# Assume the maximum possible date is 310 days
max_day <- 310


# t: actual date of death vector
# deaths: the number of deaths per day
# n.rep: indicates the total number of iterations
# bs: logical value indicating whether to bootstrapping
# t0: initial infection time vector
# return list(P, inft, t0)
deconv <- function(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL) {
  n <- sum(deaths)

  # Define the parameters of a lognormal distribution
  meanlog <- 3.152
  sdlog <- 0.451
  max_duration <- 80 # maximum number of days in the distribution

  # Generates a distribution of duration from infection to death
  infection_to_death_probs <- generate_infection_to_death_dist(meanlog, sdlog, max_duration)

  if (is.null(t0)) {
    t0 <- generate_initial_infection_times(deaths, julian_days, infection_to_death_probs)
    # t0 has a very small probability of being negative, deal with it later!
  }

  # Initializes a vector of length n.rep
  P_history <- numeric(n.rep)
  # Initialize the matrix of max_day * n.rep
  inft <- matrix(0, nrow = max_day, ncol = n.rep)

  # main loop: n.rep times
  for (rep in 1:n.rep) {
    # Sample a new duration for each individual and calculate the simulated time of death
    # put this in the loop can decrease P value further
    durations <- sample(1:max_duration, size = n, replace = TRUE, prob = infection_to_death_probs)

    simulated_death_times <- t0 + durations

    # Calculate the simulated daily death toll
    simulated_deaths <- tabulate(simulated_death_times, nbins = max_day)

    # Calculate the fitness P
    P <- calculate_P(actual_deaths = deaths, simulated_deaths = simulated_deaths)

    # Store the P-value and the distribution of infected people before optimization
    P_history[rep] <- P
    inft[, rep] <- tabulate(t0, nbins = max_day) # matrix: max_day * rep

    # Randomly shuffle the index order
    indices <- sample(1:n, n)

    for (i in indices) {
      # Choose an adjustment step
      if (rep < n.rep / 2) {
        steps_set <- c(-8, -4, -2, -1, 1, 2, 4, 8)
      } else if (rep < n.rep * 3 / 4) {
        steps_set <- c(-4, -2, -1, 1, 2, 4)
      } else {
        steps_set <- c(-2, -1, 1, 2)
      }

      step <- sample(steps_set, size = 1)
      proposed_t0_i <- t0[i] + step

      # Make sure the proposed time of infection is legal
      if (proposed_t0_i < 1 || proposed_t0_i > max_day) next

      # Calculate the current and proposed date of death
      current_death_time <- t0[i] + durations[i]
      proposed_death_time <- proposed_t0_i + durations[i]

      # Make sure the date of death is legal
      if (current_death_time < 1 || current_death_time > max_day) next
      if (proposed_death_time < 1 || proposed_death_time > max_day) next

      # Extend the deaths vector to length max_day, filling it with zeros
      deaths_extended <- c(deaths, rep(0, max_day - length(deaths)))

      # Get the actual number of deaths on the affected dates
      di_current <- deaths_extended[current_death_time]
      di_proposed <- deaths_extended[proposed_death_time]

      # Get simulated deaths for affected dates (before adjustment)
      dsi_current <- simulated_deaths[current_death_time]
      dsi_proposed <- simulated_deaths[proposed_death_time]

      # Calculate the adjusted simulated death toll
      dsi_new_current <- dsi_current - 1
      dsi_new_proposed <- dsi_proposed + 1

      # Make sure the simulated death toll is not negative
      if (dsi_new_current < 0) next

      # Calculate the P-value contribution before adjustment
      P_before_current <- (di_current - dsi_current)^2 / pmax(1, dsi_current)
      P_before_proposed <- (di_proposed - dsi_proposed)^2 / pmax(1, dsi_proposed)
      P_before <- P_before_current + P_before_proposed

      # Calculate the adjusted P-value contribution
      P_after_current <- (di_current - dsi_new_current)^2 / pmax(1, dsi_new_current)
      P_after_proposed <- (di_proposed - dsi_new_proposed)^2 / pmax(1, dsi_new_proposed)
      P_after <- P_after_current + P_after_proposed

      # Calculate ΔP
      delta_P <- P_after - P_before

      # Decide whether to accept the offer
      if (delta_P < 0) {
        # Accept proposal to update t0[i]
        t0[i] <- proposed_t0_i

        # Update the affected dates in simulated_deaths
        simulated_deaths[current_death_time] <- dsi_new_current
        simulated_deaths[proposed_death_time] <- dsi_new_proposed

        # Update the total p-value
        P <- P + delta_P
      }
    }
  }

  return(list(P = P_history, inft = inft, t0 = t0))
}

ls = deconv(julian_days, deaths)

ls

# meanlog <- 3.152
# sdlog <- 0.451
# max_duration <- 80 # maximum number of days in the distribution

# # Generates a distribution of duration from infection to death
# infection_to_death_probs <- generate_infection_to_death_dist(meanlog, sdlog, max_duration)

# t0 <- generate_initial_infection_times(deaths, julian_days, infection_to_death_probs)

# t0
# nolint end