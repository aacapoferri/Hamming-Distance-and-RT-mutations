########################################################################
# Lambda of time for different viral genome lengths, epsilon, or time #
#                          Version 1.0                                 #
#                Author: Adam A. Capoferri, PhD                        #
#               Contact: adam.capoferri@nih.gov                        #
########################################################################


##### PART I: lambda(t) where G changes #####
##### PART II: lambda(t) where epsilon changes #####
##### PART III: lambda(t) where time changes (add' info) #####

# Introduction information:
# Reproductive number here is set to 6 from Stafford et al. 2000 based on the timing of HIV acquisition. Note this number may change based on timing.
# The single round of HIV-1 RT mutation rate, epsilon, is a scaling factor (Mansky and Temin)
# "Sequenced" genome length, or G, allows to explore multiple lengths to compare estimated lambdas of 500bp vs 3000bp for instance. It is agnostic to what region you are sequencing
# The time range (1-30) is in "days", is to capture the acute phase of HIV infection during the early Fiebig stages. Around 30 days, the adaptive immune response starts to kick in, and thus will nullify the model here. This script is designed mostly past the eclipse phase, starting with Fiebig I to Fiebig III and probably IV (11-30 days)


##### PART I: lambda(t) where G changes #####

# Load necessary libraries
library(ggplot2)
library(viridis)  # For color palettes

# Define constants
R_0 <- 6  # Reproductive number as per Stafford et al. 2000
epsilon_values <- c(3e-5)  # Multiple scaling factors (Mansky and Temin)

# Define the function for lambda(t)
calculate_lambda <- function(G, t, epsilon) {
  phi <- sqrt(1 + 8 / R_0)
  lambda_t <- epsilon * G * (5/8 * t * (1 + phi) / phi + (1 - phi) / (phi^2)) # (Keele 2008)
  return(lambda_t)
}

# Hardcoded input values for testing
G_values <- c(234, 500, 1000)  # Vector of G values ("sequenced" genome lengths)
t_input <- "1:30"  # Time range input of infection (days). Max of 30 days due to adaptive immune response.

# Helper function to parse range or list input
parse_input <- function(input_str) {
  if (grepl(":", input_str)) {
    parsed_values <- eval(parse(text = input_str))
    return(parsed_values)
  } else {
    parsed_values <- as.numeric(unlist(strsplit(input_str, ",")))
    return(parsed_values)
  }
}

# Parse the time values (t)
t_values <- parse_input(t_input)

# Calculate lambda(t) for each G, each t, and each epsilon, and store in a data frame
lambda_data <- data.frame()

for (epsilon in epsilon_values) {
  for (G in G_values) {
    temp_data <- data.frame(
      Time = t_values,
      Lambda = sapply(t_values, function(t) calculate_lambda(G, t, epsilon)),
      G = as.factor(G),
      Epsilon = as.factor(epsilon)  # Include epsilon as a grouping factor
    )
    lambda_data <- rbind(lambda_data, temp_data)  # Combine data
  }
}

# Print lambda values for debugging
print(lambda_data)

# Create the plot with Viridis "D" colors and solid lines
lambda_plot <- ggplot(lambda_data, aes(x = Time, y = Lambda, color = G)) +  # Color by G
  geom_line(linewidth = 1, linetype = "solid") +  # All lines solid
  scale_y_log10(limits = c(1e-3, 10), breaks = scales::log_breaks(n = 5)) +  # Use log10 ticks on y-axis
  scale_x_continuous(breaks = c(1,5,10,15,20,25,30),
                     labels = c("1","5","10","15","20","25","30")) +  # Set x-axis ticks in increments of 5
  scale_color_viridis(discrete = TRUE, option = "D") +  # Use Viridis D palette
  labs(
    title = "Lambda(t) for Different Genome Sizes",
    x = "Time since transmission (days)",
    y = "Log10 Estimated Lambda",
    color = "Genome size, G (bp)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.line = element_line(color = "black")  # Set x and y axes to black
  )

# Save the plot as PDF
ggsave("lambda_plot_multi_epsilon.pdf", plot = lambda_plot, width = 8, height = 6)

# Save the lambda_data data frame to a CSV file
write.csv(lambda_data, "lambda_data_multi_epsilon.csv", row.names = FALSE)

### END ###




##### PART II: lambda(t) where epsilon changes #####

# Load necessary libraries
library(ggplot2)
library(viridis)  # For color palettes

# Define constants
R_0 <- 6  # Reproductive number as per Stafford et al. 2000
epsilon_values <- c(3.0e-5, 6.29e-6, 4.66e-6, 6.1e-6, 9.49e-6, 5.89e-6, 4.7e-6, 3.74e-6, 4.72e-6, 5.47e-6, 3.52e-6)  # Multiple scaling factors

# Define the function for lambda(t)
calculate_lambda <- function(G, t, epsilon) {
  phi <- sqrt(1 + 8 / R_0)
  lambda_t <- epsilon * G * (5/8 * t * (1 + phi) / phi + (1 - phi) / (phi^2))
  return(lambda_t)
}

# Hardcoded input values for testing
G_values <- c(234)  # Vector of G values ("sequenced" genome lengths)
t_input <- "1:30"  # Time range input of infection (days). Max of 30 days due to adaptive immune response.

# Helper function to parse range or list input
parse_input <- function(input_str) {
  if (grepl(":", input_str)) {
    parsed_values <- eval(parse(text = input_str))
    return(parsed_values)
  } else {
    parsed_values <- as.numeric(unlist(strsplit(input_str, ",")))
    return(parsed_values)
  }
}

# Parse the time values (t)
t_values <- parse_input(t_input)

# Calculate lambda(t) for each G, each t, and each epsilon, and store in a data frame
lambda_data <- data.frame()

for (epsilon in epsilon_values) {
  for (G in G_values) {
    temp_data <- data.frame(
      Time = t_values,
      Lambda = sapply(t_values, function(t) calculate_lambda(G, t, epsilon)),
      G = as.factor(G),
      Epsilon = as.factor(epsilon)  # Include epsilon as a grouping factor
    )
    lambda_data <- rbind(lambda_data, temp_data)  # Combine data
  }
}

# Print lambda values for debugging
print(lambda_data)

# Create the plot, including Epsilon as a grouping and coloring factor
lambda_plot <- ggplot(lambda_data, aes(x = Time, y = Lambda, color = Epsilon)) +  # Color by Epsilon
  geom_line(linewidth = 1, linetype = "solid") +  # Use solid lines for all curves
  scale_y_log10(
    breaks = 10^seq(-4, 1, by = 1),  # Set breaks from 10^-4 to 10^1
    labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use math formatting
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20, 25, 30), 
    labels = as.character(c(1, 5, 10, 15, 20, 25, 30))
  ) +
  scale_color_viridis(discrete = TRUE, option = "D", name = "Epsilon") +  # Use Viridis palette for Epsilon
  labs(
    title = expression(Lambda(t) ~ "for Different Epsilon Values"),  # Enhanced title formatting
    x = "Time since transmission (days)",
    y = expression(log[10] ~ "Estimated Lambda"),  # Y-axis with log base 10 label
    color = "Epsilon (scaling factor)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(hjust = 0.5),  # Center-align the title
    axis.line = element_line(color = "black"),  # Set x and y axes to black
    axis.text.x = element_text(size = 10)  # Adjust x-axis label size
  )

# Save the updated plot as PDF
ggsave("lambda_plot_multi_epsilon_pol.pdf", plot = lambda_plot, width = 8, height = 6)

# Save the lambda_data data frame to a CSV file
write.csv(lambda_data, "lambda_data_multi_epsilon_pol.csv", row.names = FALSE)


##### PART III: lambda(t) where time changes #####

# Calculate the estimated lambda given time AND genome length sequenced. This is a bit more of a "combined" version of the prior two parts where you can account for the genome length and time component in one figure.

#The lambda is then used within a Poisson distribution of: Pr(k>/=0) = (lambda^k*exp(-lambda))/k! This states that given an estimated lambda(G,t) with a Hamming distance value of k, the probability may be calculated. The script is written to calculated these individual estimated lambda (mean) values, then take range of Hamming Distance that are wished to be looked at, and calculate the Probability with respect to the estimated lambda and each k-term.

# Information needed:
  # Reproductive number here is set to 6 from Stafford et al. 2000 based on the timing of infection. Note this number may change, consult relevant articles before embarking.
  # The single round of HIV-1 RT mutation rate, epsilon, is a scaling factor
  # "Sequenced" genome length, or G, allows you to explore multiple lengths to compare estimated lambdas of 500bp vs 3000bp for instance. It is agnostic to what region you are sequencing.
  # The time range (1-30 days), is to capture the acute phase of infection during the early Fiebig stages. Around 30 days, the adaptive immune response starts to kick in, and thus will nullify the model here. This script is designed mostly past the eclipse phase, starting with Fiebig I to Fiebig III and probably IV (11-30 days)
  # How far out with Hamming Distance (0-n) are you interested in calculating the probability.


# Load necessary libraries
library(ggplot2)
library(viridis)

# Define constants
R_0 <- 6  # Reproductive number as per Stafford et al.
epsilon <-  3e-5 # Scaling factor 3e-5 (Mansky and Temin)

# Define the function for lambda(t)
calculate_lambda <- function(G, t) {
  phi <- sqrt(1 + 8 / R_0)
  lambda_t <- epsilon * G * (5/8 * t * (1 + phi) / phi + (1 - phi) / (phi^2)) (Keele 2008)
  return(lambda_t)
}

# Hardcoded input values for testing
G_values <- c(234, 500)  # Vector of G values
t_values <- c(1,10, 12, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 28, 30)  # Specific time points

# Calculate lambda(t) for each G and each t and store in a data frame
lambda_data <- data.frame()
for (G in G_values) {
  temp_data <- data.frame(Time = t_values, Lambda = sapply(t_values, function(t) calculate_lambda(G, t)))
  temp_data$G <- as.factor(G)  # Add a column for G values
  lambda_data <- rbind(lambda_data, temp_data)  # Combine data
}

# Print lambda values for debugging
print(lambda_data)

# Plot lambda(t) over the range of t using log10 transformation on y-axis
lambda_plot <- ggplot(lambda_data, aes(x = Time, y = Lambda, color = as.factor(Time), shape = G)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +  # Increase point size for visibility
  scale_y_log10(limits = c(1e-3, 10), breaks = scales::log_breaks(n = 5)) +  # Use log10 ticks on y-axis
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +  # Set x-axis ticks in increments of 5. If changing the t_input to less than 30 days, ensure to scale this properly.
  scale_color_viridis(discrete = TRUE) +  # Use viridis color palette for time
  labs(title = "",
       x = "λ(t) (days)",
       y = "Log10 Estimated λ",
       color = "Time (days)",  # Legend title for color
       shape = "Genome size, G (bp)") +  # Legend title for shape
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        plot.title = element_text(hjust = 0.5),  # Center the title
        axis.line = element_line(color = "black")) +  # Set x and y axes to black
  guides(shape = guide_legend(override.aes = list(size = 3)))  # Adjust shape size in legend


# Define a function to calculate the PDF of the Hamming distance
calculate_pdf <- function(lambda, k) {
  return((lambda^k * exp(-lambda)) / factorial(k))
}

# Define the Hamming distances (as a vector)
hamming_distances <- 0:10  # Limited to 0 to n; change as necessary.

# Create a data frame to store results for PDF vs Hamming distance
pdf_data <- data.frame()

# Calculate the PDF for each combination of lambda and Hamming distance
for (G in G_values) {
  for (t in t_values) {
    lambda <- calculate_lambda(G, t)
    pdf_values <- sapply(hamming_distances, calculate_pdf, lambda = lambda)  # Calculate PDF for current lambda
    
    # Create a temporary data frame for the current G and t
    temp_data <- data.frame(Hamming_Distance = hamming_distances,
                            PDF = pdf_values,
                            Lambda = as.factor(lambda), 
                            Time = as.factor(t),  # Include time in the data
                            G = as.factor(G))  # Include genome size as a factor
    
    # Combine the temporary data frame with the main data frame
    pdf_data <- rbind(pdf_data, temp_data)
  }
}

# Plot the results with the Viridis "D" color palette
pdf_plot <- ggplot(pdf_data, aes(x = Hamming_Distance, y = PDF, color = Time, shape = G)) +
  geom_line(size = .5) +
  geom_point(size = 2) +  # Add points to represent genome size
  labs(title = "Probability Density Function of Hamming Distance",
       x = "Hamming Distance (k)",
       y = "Probability Density Function P(k ≥ 0)",
       color = "Time (days)",  # Legend title for color
       shape = "Genome size, G (bp)") +  # Legend title for shape
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +  # Set y-axis limit to 0-1
  scale_x_continuous(breaks = 0:10, limits = c(0, 10), expand = c(0, 0)) +  # Show x-axis from 0 to n. Change the limits and breaks based on the Hamming distance you wish to examine.
  scale_color_viridis_d(option = "D") +  # Use the viridis D color palette
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black")) +  # Set axes to black
  guides(shape = guide_legend(override.aes = list(linewidth = 3)))  # Adjust shape size in legend

# Save the PDF plot as a PDF, change file name as needed
ggsave("filename.pdf", plot = pdf_plot, width = 8, height = 6)
write.csv(pdf_data, "pdf_and_hamming_distance test.csv", row.names = FALSE)

### END ###
###########
