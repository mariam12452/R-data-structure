#DNA to protein ttranslation function
translate_dna_to_protein <- function(dna_seq) {
  codon_table <- list(
    "UUU"="F", "UUC"="F",  # Phenylalanine
    "UUA"="L", "UUG"="L",  # Leucine
    "UCU"="S", "UCC"="S", "UCA"="S", "UCG"="S",  # Serine
    "UAU"="Y", "UAC"="Y",  # Tyrosine
    "UAA"="STOP", "UAG"="STOP", "UGA"="STOP",  # Stop codons
    "UGU"="C", "UGC"="C",  # Cysteine
    "UGG"="W",  # Tryptophan
    "CUU"="L", "CUC"="L", "CUA"="L", "CUG"="L",  # Leucine
    "CCU"="P", "CCC"="P", "CCA"="P", "CCG"="P",  # Proline
    "CAU"="H", "CAC"="H",  # Histidine
    "CAA"="Q", "CAG"="Q",  # Glutamine
    "CGU"="R", "CGC"="R", "CGA"="R", "CGG"="R",  # Arginine
    "AUU"="I", "AUC"="I", "AUA"="I",  # Isoleucine
    "AUG"="M",  # Methionine (Start codon)
    "ACU"="T", "ACC"="T", "ACA"="T", "ACG"="T",  # Threonine
    "AAU"="N", "AAC"="N",  # Asparagine
    "AAA"="K", "AAG"="K",  # Lysine
    "AGU"="S", "AGC"="S",  # Serine
    "AGA"="R", "AGG"="R",  # Arginine
    "GUU"="V", "GUC"="V", "GUA"="V", "GUG"="V",  # Valine
    "GCU"="A", "GCC"="A", "GCA"="A", "GCG"="A",  # Alanine
    "GAU"="D", "GAC"="D",  # Aspartic acid
    "GAA"="E", "GAG"="E",  # Glutamic acid
    "GGU"="G", "GGC"="G", "GGA"="G", "GGG"="G"   # Glycine
  )
  # Convert DNA to RNA (replace T with U)
  rna_seq <- gsub("T", "U", dna_seq)
  
  # Split RNA sequence into codons (triplets)
  codons <- substring(rna_seq, seq(1, nchar(rna_seq) - 2, 3), seq(3, nchar(rna_seq), 3))
  
  # Translate codons into amino acids
  protein <- unlist(lapply(codons, function(codon) codon_table[[codon]]))
  
  # Stop translation at the first "STOP" codon
  if ("STOP" %in% protein) {
    stop_index <- which(protein == "STOP")[1]  # Get the first stop codon index
    protein <- protein[1:(stop_index - 1)]  # Keep only amino acids before stop
  }
  
  # Return the final protein sequence as a string
  return(paste(protein, collapse = ""))
}
translate_dna_to_protein("AUGGCUAGGAUCGCCAUGUAA")

#logistic growth simulation function
library(ggplot2)
library(dplyr)

# Define function for generating a logistic growth curve
simulate_logistic_growth <- function(time_steps = 100, P0 = 10, K = 1000, r = 0.2) {
  lag_time <- sample(1:10, 1)   # Random lag phase duration
  exp_time <- sample(10:30, 1)  # Random exponential phase duration
  
  # Initialize population vector
  population <- numeric(time_steps)
  population[1] <- P0  
  
  # Simulate logistic growth over time
  for (t in 2:time_steps) {
    if (t <= lag_time) {
      population[t] <- population[t-1] + r * population[t-1] * 0.1
    } else if (t <= exp_time) {
      population[t] <- population[t-1] + r * population[t-1] * (1 - (population[t-1] / K))
    } else {
      population[t] <- population[t-1] + r * population[t-1] * (1 - (population[t-1] / K)) * 0.5
    }
  }
  
  # Return a dataframe
  return(data.frame(Time = 1:time_steps, Population = population))
}

# Generate a dataframe with 100 different growth curves
generate_multiple_growth_curves <- function(num_simulations = 100, time_steps = 100) {
  all_growth_data <- data.frame()  # Initialize an empty dataframe
  
  for (i in 1:num_simulations) {
    growth_data <- simulate_logistic_growth(time_steps = time_steps, 
                                            P0 = sample(5:15, 1),   # Varying initial population
                                            K = sample(900:1100, 1), # Slightly varying carrying capacity
                                            r = runif(1, 0.1, 0.3))  # Random growth rate
    
    growth_data$Simulation_ID <- i  # Add a column to label each simulation
    all_growth_data <- bind_rows(all_growth_data, growth_data)  # Combine results
  }
  
  return(all_growth_data)
}

# Run function to generate dataset
growth_curves_df <- generate_multiple_growth_curves()

# Select a single growth curve to plot (e.g., Simulation_ID = 1)
selected_curve_id <- 1
selected_curve <- growth_curves_df %>% filter(Simulation_ID == selected_curve_id)

# Plot only one selected growth curve
ggplot(selected_curve, aes(x = Time, y = Population)) +
  geom_line(color = "blue", size = 1) +
  labs(title = paste("Growth Curve for Simulation ID:", selected_curve_id), 
       x = "Time", y = "Population") +
  theme_minimal()

find_threshold_time <- function(growth_data, K = NULL, threshold_percent = 0.8) {
  if (is.null(K)) {
    K <- max(growth_data$Population)  # Assume K is the max population reached
  }
  
  threshold <- threshold_percent * K  # Calculate 80% of K
  
  # Find the first time where the population reaches or exceeds the threshold
  threshold_time <- min(growth_data$Time[growth_data$Population >= threshold], na.rm = TRUE)
  
  return(threshold_time)  # Returns the first time step where the threshold is met
}
# Select a single growth curve (e.g., Simulation_ID = 1)
selected_curve <- growth_curves_df %>% filter(Simulation_ID == 1)

# Find the time when it reaches 80% of carrying capacity
threshold_time <- find_threshold_time(selected_curve)

print(paste("Time to reach 80% of K:", threshold_time))

head(growth_curves_df)

#function for calculating Hamming distance
library(stringr)

calculate_hamming_distance <- function(string1, string2) {
  # Step 1: Make both strings the same length by padding with spaces
  max_length <- max(nchar(string1), nchar(string2))  # Find the longest string length
  string1 <- str_pad(string1, max_length, side = "right", pad = " ")  # Pad shorter string
  string2 <- str_pad(string2, max_length, side = "right", pad = " ")
  
  # Step 2: Convert strings to individual characters
  char1 <- unlist(strsplit(string1, ""))  
  char2 <- unlist(strsplit(string2, ""))  
  
  # Step 3: Compare characters and count differences
  hamming_distance <- sum(char1 != char2)  
  
  return(hamming_distance)  # Return the computed Hamming Distance
}

calculate_hamming_distance("mariamfolake","maraimfolakemi")



