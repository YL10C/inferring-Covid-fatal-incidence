## Xiaotian Xing s2661110;
##

## 1
data <- read.table("engcov.txt", header = TRUE)

deaths <- data$nhs[1:150]
julian_days <- data$julian[1:150]

generate_infection_to_death_dist <- function(meanlog, sdlog, max_duration) {
  probabilities <- dlnorm(1:max_duration, meanlog, sdlog)
  probabilities / sum(probabilities)
}

meanlog <- 3.152
sdlog <- 0.451
max_duration <- 80
infection_to_death_probs <- generate_infection_to_death_dist(meanlog, sdlog, max_duration)


## 2
generate_initial_infection_times <- function(deaths, julian_days, infection_to_death_probs) {
  n <- sum(deaths)
  death_dates <- rep(julian_days, deaths)
  infection_to_death_days <- sample(1:length(infection_to_death_probs), n, replace = TRUE, prob = infection_to_death_probs)
  infection_dates <- death_dates - infection_to_death_days
  return(infection_dates)
}

initial_infection_times <- generate_initial_infection_times(deaths, julian_days, infection_to_death_probs)

