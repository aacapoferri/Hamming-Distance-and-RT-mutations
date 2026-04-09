#############################################################
# Monte Carlo simulation of mutation accumulation over time #
#                        Version 1.0                        #
#                Author: Adam A. Capoferri, PhD             #
#                Contact: adam.capoferri@nih.gov            #
#############################################################

# Introduction: This performs a Monte Carlo simulation to estimate the accumulation of mutations in sequenced genomes over time. It follows a Poisson distribution to model daily mutations, accounts for assay error, and visualizes the results as mutation counts over different time spans. Simulate mutation accumulation and observe how mutation counts vary over time with different genome sizes. This can be useful for understanding mutation dynamics and sequencing error impacts over time. Acute HIV-1 infection is the pathogen in mind when developing. Detailed in Capoferri and Boltz et al. (manuscript in prep.).

# Informational inputs are as follows:
  # "Sequenced" genome length, or G, allows you to explore multiple lengths to compare estimated lambdas of 500bp vs 3000bp for instance. It is agnostic to what region you are sequencing.
  # The time range (1-30 days), is to capture the acute phase of infection during the early Fiebig stages. Around 30 days, the adaptive immune response starts to kick in, and thus will nullify the model here. This script is designed mostly past the eclipse phase, starting with Fiebig I to Fiebig III and probably IV (11-30 days)
  # Number of genomes as the starting population size
  # The single round of HIV-1 RT mutation rate, epsilon, is a scaling factor (Mansky and Temin)
  # Number of Monte Carlo Iterations to run (currently set to 100,000)
  # Assay error-rate noise of 0.03% based on Palmer et al. Coffin (Journal of Clinical Microbiology, 2005)


# Load libraries
library(ggplot2)
library(viridis)
library(dplyr)

# Function to calculate expected mutations with Monte Carlo simulation
calculate_mutations_monte_carlo <- function(genome_size_bp, num_genomes, mutation_rate_per_bp, days, assay_error_rate) {
  # Calculate total number of bases sequenced
  total_bases_sequenced <- genome_size_bp * num_genomes
  
  # Calculate the expected number of mutations per day
  expected_mutations_per_day <- mutation_rate_per_bp * total_bases_sequenced
  
  # Simulate mutations for each day using Poisson distribution
  total_mutations <- sum(rpois(days, expected_mutations_per_day))
  
  # Apply assay error rate (once added)
  total_mutations <- total_mutations + round(total_bases_sequenced * assay_error_rate)
  
  # Return total mutations and total bases sequenced
  return(list(total_mutations = total_mutations, total_bases_sequenced = total_bases_sequenced))
}

# Function to run Monte Carlo simulation for n iterations
run_monte_carlo_simulation <- function(n_iterations, genome_size_bp, num_genomes, mutation_rate_per_bp, days, assay_error_rate) {
  # Initialize vector to store mutation counts for each iteration
  mutations_results <- numeric(n_iterations)
  
  # Loop through n iterations
  for (i in 1:n_iterations) {
    # Calculate mutations for this iteration using Monte Carlo simulation
    results <- calculate_mutations_monte_carlo(genome_size_bp, num_genomes, mutation_rate_per_bp, days, assay_error_rate)
    
    # Store total mutations
    mutations_results[i] <- results$total_mutations
  }
  
  # Return mutation counts
  return(mutations_results)
}

# Parameters for multiple genome sizes and days
genome_sizes <- c(234)  # Different genome sizes in base pairs
days_list <- c(1:30)         # Different days of infection
num_genomes <- 520        # Number of genomes sequenced
mutation_rate_per_bp <- 3E-05 # RT mutation rate (3 mutations per 100,000 bases per day)
n_iterations <- 100000       # Number of Monte Carlo iterations to run
assay_error_rate <- 0.0003   # Assay error rate of 0.03%

# Initialize data frame to store results
all_results <- data.frame()

# Run the simulation for each combination of genome size and days
for (genome_size in genome_sizes) {
  daily_means <- numeric()  # Store daily mean mutations for each day
  
  for (days in days_list) {
    # Run the Monte Carlo simulation
    mutation_counts <- run_monte_carlo_simulation(n_iterations, genome_size, num_genomes, mutation_rate_per_bp, days, assay_error_rate)
    
    # Calculate mean mutation count for the current day value
    daily_mean_mutations <- mean(mutation_counts)
    daily_means <- c(daily_means, daily_mean_mutations)
  }
  
  # Calculate slope of mean mutations over time using linear regression
  lm_result <- lm(daily_means ~ days_list)
  slope <- round(coef(lm_result)["days_list"], 0)  # Extract and round slope to 1 decimal place
  
  # Append each combination’s summary statistics and slope
  summary_stats <- data.frame(
    Genome_Size = genome_size,
    Days = days_list,
    Mean_Mutations = daily_means,
    Slope = slope  # Include the calculated slope
  )
  
  # Append to main results data frame
  all_results <- rbind(all_results, summary_stats)
}

# Plotting the results without grid lines but with axes lines
mutation_plot <- ggplot(all_results, aes(x = Days, y = Mean_Mutations, color = factor(Genome_Size))) +
  geom_line(size = 1) +
  geom_point(size = 0) +
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30), 
                     labels = c("1", "5", "10", "15", "20", "25", "30")) +
  scale_y_continuous(breaks = seq(0, max(6000, max(all_results$Mean_Mutations)), by = 500)) + # Ensure inclusion of 6000
  labs(title = "Expected Number of Mutations Over Time",
       x = "Time since transmission (days)",
       y = "Expected Number of Mutations",
       color = "Genome Size (bp)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black")) +  # Retain axis lines
  geom_text(data = all_results %>% group_by(Genome_Size) %>% slice(1),
            aes(label = paste(Slope, "mut/day")),
            hjust = -0.1, vjust = 1.5, size = 3.5, color = "black")

# Save plot as PDF, change as needed
ggsave("Expected_Mutations_vs_Time.pdf", mutation_plot, width = 8, height = 6)

# Print the plot
print(mutation_plot)
print(all_results)

# Save the results to a CSV file, change as needed
write.csv(all_results, "Monte_Carlo_Mutations_Results.csv", row.names = FALSE)

### END ###

# Additional Information:
  # Expected Mutations per Day: expected_mutations_per_day=mutation_rate_per_bp×(genome_size_bp×num_genomes).This calculates the average number of mutations expected per day based on the mutation rate per base pair and the total number of bases sequenced across all genomes.
  # Poisson Simulation of Daily Mutations: Each day, the script simulates the actual mutation count using a Poisson distribution with the calculated expected mutations per day as the mean. This simulates the natural variability in mutation counts around the expected average.
  # Assay Error Adjustment: total_mutations=total_mutations+round(total_bases_sequenced×assay_error_rate). After calculating the mutations over time, a fraction is added to account for sequencing errors based on a specified assay error rate.
  # Monte Carlo Simulation: To obtain a stable estimate of mutation accumulation, the script repeats this daily simulation multiple times (n_iterations) for each parameter combination. It then calculates the average mutation count across all iterations for each time point.
  # Slope Calculation (Mutation Accumulation Rate): The script performs a linear regression on the mean mutation counts over days to find the slope. Where slope= delta(mean mutation count)/delta(days) with the slope indicate the rate of mutation accumulation per day.
