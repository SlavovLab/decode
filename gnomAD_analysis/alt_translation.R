library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(tidyverse)
library(magrittr)
library(viridis)

#############
# Decide here if you want to plot 
#
# All mutations:
# exact_mutation = T
# other_mutation = T
#
# Exact mutations:
# exact_mutation = T
# other_mutation = F
#
# All but exact mutations:
# exact_mutation = F
# other_mutation = T
#############

exact_mutation = T # missense mutation that results in observed AAS
other_mutation = T # other missense mutations in codon


# Which plot to run
plot_sfs = T # figures 5f and extended data 10e
plot_constraint = T # figures 5e and extended data 10d

# Create a named vector where each amino acid corresponds to a list of codons
codon_map <- list(
  A = c("GCT", "GCC", "GCA", "GCG"),
  R = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
  N = c("AAT", "AAC"),
  D = c("GAT", "GAC"),
  C = c("TGT", "TGC"),
  Q = c("CAA", "CAG"),
  E = c("GAA", "GAG"),
  G = c("GGT", "GGC", "GGA", "GGG"),
  H = c("CAT", "CAC"),
  I = c("ATT", "ATC", "ATA"),
  L = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
  K = c("AAA", "AAG"),
  M = c("ATG"),
  F = c("TTT", "TTC"),
  P = c("CCT", "CCC", "CCA", "CCG"),
  S = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
  T = c("ACT", "ACC", "ACA", "ACG"),
  W = c("TGG"),
  Y = c("TAT", "TAC"),
  V = c("GTT", "GTC", "GTA", "GTG")
)


compute_ratio_and_ci <- function(filtered_table) {
  obs <- sum(filtered_table$polymorphic)
  exp <- sum(filtered_table$exp)
  ratio <- obs / exp
  ci <- generate_ci(obs = obs, exp = exp, alpha = 0.025)
  
  return(list(ratio = ratio, ci = ci))
}

generate_ci = function(obs, exp, alpha=0.025) {
  distr = data.frame(x = seq(0, 2, by=0.001)) %>%
    mutate(y = dpois(obs, x*exp),
           y_cumulative = cumsum(y)) %>%
    mutate(y_cumulative_norm = y_cumulative/max(y_cumulative))
  
  low =  distr %>% dplyr::filter(y_cumulative_norm < alpha) %$% max(x)
  high = distr %>% dplyr::filter(y_cumulative_norm > 1 - alpha) %$% min(x)
  return(data.frame(lower=low, upper=high))
}

# Function to convert string to vector
convert_to_vector <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  x <- strsplit(x, ",")[[1]]
  as.numeric(x)
}

str_replace_at <- function(string, position, replacement) {
  substr(string, position, position) <- replacement
  return(string)
}




raas = read.delim("raas_table.tsv", header = TRUE, sep = "\t")
case_missense <- read.delim("missense_table.tsv.bgz", header = TRUE, sep = "\t")

raas <- raas %>%
  separate(ID, into = c("Chromosome", "Position", "Strand"), sep = ":", remove = FALSE)

# Convert Chromosome and Position columns to numeric
raas <- raas %>%
  mutate(Chromosome = as.numeric(Chromosome),
         Position = as.numeric(Position))

# Reorder the dataframe by Chromosome and then by Position
raas <- raas %>%
  arrange(Chromosome, Position)

# If you want to drop the separate columns and keep only the original ID column
raas <- raas %>%
  dplyr::select(ID, everything())


case_missense$MMSeq2 = as.numeric(case_missense$MMSeq2)


case_missense$exp <- 1-as.numeric(sub("^.*0\\.0:([0-9\\.]+),.*$", "\\1", case_missense$mu_annotation))
case_missense$RAAS_quartile <- cut(case_missense$RAAS.median, breaks = quantile(case_missense$RAAS.median, probs = seq(0, 1, 0.25)), include.lowest = TRUE)

df <- case_missense %>%
  separate_rows(fromto, sep = ";") %>%
  mutate(second_aa = sub(".*:(.*)", "\\1", fromto)) %>%
  mutate(codons = lapply(second_aa, function(aa) codon_map[[aa]])) 

df2 <- df %>%
  mutate(alt2 = case_when(
    strand == "+" ~ alt,
    strand == "-" ~ chartr("ACGT", "TGCA", alt)
  ))%>%
  mutate(pos_codon2 = case_when(
    strand == "+" ~ pos_codon,
    strand == "-" & pos_codon == 1 ~ 3,
    strand == "-" & pos_codon == 3 ~ 1,
    strand == "-" & pos_codon == 2 ~ 2
  ))

df3 <- df2 %>%
  mutate(codon_mut = str_replace_at(codon, pos_codon2, alt2))



raas_not_in_case_missense <- raas %>%
  anti_join(case_missense, by = "ID")



df4 <- df3 %>%
  distinct(locus, ref, alt, fromto, .keep_all = TRUE)

raas_in_case_missense2 <- df3 %>%
  distinct(ID, .keep_all = TRUE)


raas_in_case_missense_non_poly <- df3 %>%
  filter(polymorphic == 0) %>%
  distinct(locus,ref, alt, fromto, .keep_all = TRUE)

raas_in_case_missense_poly <- df3 %>%
  filter(polymorphic > 0) %>%
  distinct(locus,ref, alt, fromto, .keep_all = TRUE)


df_test <- df4 %>%
  rowwise() %>%
  mutate(exact_mutation = codon_mut %in% codons) %>%
  ungroup() %>%
  distinct() %>%
  mutate(mutated_pos_codon = pos_codon2, codons_matching_alt_translation = codons,
         mutated_base = alt2, mutated_codon = codon_mut, position_in_codon_to_mutate = pos_codon2) %>%
  dplyr::select(ID, locus, transcript, RAAS.median, codon, fromto, 
         MMSeq2, ref, alt, mutated_pos_codon, AC_gnomADv4, AF_gnomADv4,polymorphic, exp,
         codons_matching_alt_translation, mutated_base, position_in_codon_to_mutate, 
         mutated_codon, exact_mutation, mu_annotation)

df_test <- df_test %>%
  separate(ID, into = c("Chromosome", "Position", "Strand"), sep = ":", remove = FALSE)

df_test <- df_test %>%
  mutate(Chromosome = as.numeric(Chromosome),
         Position = as.numeric(Position))

df_test <- df_test %>%
  arrange(Chromosome, Position) %>%
  distinct(locus,ref, alt, fromto, .keep_all = TRUE)

df_test2 = df_test[df_test$exact_mutation,]

df_test3 = df_test2 %>%
  filter(polymorphic == 1)

df_test4 = df_test[df_test$exact_mutation == F,]



if(exact_mutation & other_mutation)
{
  print("Using all data - exact mutation and others")
  data_plot = df_test
}else if(exact_mutation & !other_mutation)
{
  print("Using only exact mutations")
  data_plot = df_test2
}else if(other_mutation & !exact_mutation)
{
  print("Using only mutations other than exact mutations")
  data_plot = df_test4
}else
{
  stop("No mutation included")
  
}






mcol <- viridis::inferno(5)
vcol <- viridis::viridis(5)
mcol[5] <- vcol[5]
arno <- colorRampPalette(mcol)

quartile_colors <- rev(arno(4))

assign_bin_af <- function(value) {
  if (value == 0) {
    return("0")
  } else if (value > 0 & value <= 1e-6) {
    return("0 to 10⁻⁶")
  } else if (value > 1e-6 & value <= 1e-5) {
    return("10⁻⁶ to 10⁻⁵")
  } else if (value > 1e-5 & value <= 1e-4) {
    return("10⁻⁵ to 10⁻⁴")
  } else if (value > 1e-4 & value <= 1e-3) {
    return("10⁻⁴ to 10⁻³")
  } else if (value > 1e-3 & value <= 1e-2) {
    return("10⁻³ to 10⁻²")
  } else if (value > 1e-2 & value <= 1e-1) {
    return("10⁻² to 10⁻¹")
  } else {
    return("10⁻¹ to 0.5")
  }
}


# Calculate the quartiles for the RAAS.median column
quartiles <- quantile(data_plot$RAAS.median, probs = c(0.25, 0.5, 0.75))

# Assign quartiles to RAAS.median
assign_quartile <- function(value) {
  if (value <= quartiles[1]) {
    return("Lowest RAAS Quartile")
  } else if (value <= quartiles[2]) {
    return("Second RAAS Quartile")
  } else if (value <= quartiles[3]) {
    return("Third RAAS Quartile")
  } else {
    return("Highest RAAS Quartile")
  }
}

# Add the RAAS quartiles and bins to the dataset
data_plot2 <- data_plot %>%
  mutate(quartile = sapply(RAAS.median, assign_quartile),
         bin = sapply(AF_gnomADv4, assign_bin_af))

# Set the correct order of quartile levels
data_plot2$quartile <- factor(data_plot2$quartile, 
                             levels = c("Highest RAAS Quartile", 
                                        "Third RAAS Quartile", 
                                        "Second RAAS Quartile", 
                                        "Lowest RAAS Quartile"))

# Summarize the data
df_combined_summary <- data_plot2 %>%
  group_by(quartile, bin) %>%
  summarise(count = n()) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  complete(bin, quartile, fill = list(count = 0, percentage = 0)) %>%
  arrange(desc(quartile))

# Set the order of bins so that "0" is at the top, and bins increase toward the bottom
df_combined_summary$bin <- factor(df_combined_summary$bin, 
                                  levels = c("10⁻¹ to 0.5", "10⁻² to 10⁻¹", 
                                             "10⁻³ to 10⁻²", "10⁻⁴ to 10⁻³", 
                                             "10⁻⁵ to 10⁻⁴", "10⁻⁶ to 10⁻⁵", 
                                             "0 to 10⁻⁶", "0"))

if(plot_sfs)
{
  # Plot the data
  ggplot(df_combined_summary, aes(x = percentage, y = bin, fill = quartile)) +
    geom_bar(stat = "identity", color = "black", size = 0.1, alpha = 0.7, 
             position = position_dodge(width = 0.65), width = 0.65) +  # Adjusted dodge and bar width
    labs(title = "Distribution of Allele Frequencies Across RAAS Quartiles",
         x = "Percentage (%)",
         y = "Allele Frequency") +
    scale_fill_manual(values = quartile_colors) +
    theme_classic() +  # Use a classic theme for a clean look
    theme(panel.grid = element_blank(),  # Remove grid lines
          axis.text.x = element_text(size = 14, face = "bold"),  # Adjusted x-axis text size and boldness
          axis.text.y = element_text(size = 14, face = "bold"),  # Adjusted y-axis text size and boldness
          axis.title.x = element_text(size = 16, face = "bold"),  # Bold axis titles
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.line = element_line(size = 0.3, color = "black"),  # Thicker axis lines
          legend.text = element_text(size = 16),  # Consistent size for legend text
          legend.title = element_blank(),  # Remove legend title
          plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),  # Bold and centered title
          legend.position = c(0.7, 0.3))  # Legend position inside the plot
  
  
}



################

# Read data and process initial columns
table_agg <- data_plot
table_agg$exp <- 1 - as.numeric(sub("^.*0\\.0:([0-9\\.]+),.*$", "\\1", table_agg$mu_annotation))

table_agg = table_agg %>% mutate(enst = transcript)

# Read additional tables and rename columns
fha_agg = read_table("fha_agg.tsv.bgz")
names(fha_agg)[1] = "enst"
fha_high_agg = read_table("fha_high_agg.tsv.bgz")
names(fha_high_agg) = c("enst", "obs_high", "exp_high")
fha_low_agg = read_table("fha_low_agg.tsv.bgz")
names(fha_low_agg) = c("enst", "obs_low", "exp_low")

# Merge the tables
case_missense = merge(table_agg, fha_agg, by="enst")
case_missense2 = merge(case_missense, fha_high_agg, by="enst")
case_missense3 = merge(case_missense2, fha_low_agg, by="enst")

# Extract RAAS values and compute quartiles
RAAS = unlist(case_missense3$RAAS)
case_missense3$RAAS_quartile <- cut(RAAS, 
                                    breaks = quantile(RAAS, probs = seq(0, 1, 0.25)), 
                                    include.lowest = TRUE)

# Convert RAAS_quartile to a factor and reverse the levels
case_missense3$RAAS_quartile <- factor(case_missense3$RAAS_quartile, 
                                       levels = rev(levels(case_missense3$RAAS_quartile)))

# Function to compute ratio and confidence interval
compute_ratio_and_ci <- function(filtered_table) {
  obs <- sum(filtered_table$AC_gnomADv4 >= 1)
  exp <- sum(as.numeric(filtered_table$exp))
  ratio <- obs / exp
  ci <- generate_ci(obs = obs, exp = exp, alpha = 0.025)
  
  return(list(ratio = ratio, ci = ci))
}

# Calculate ratios and CI for each quartile
quartile_results <- lapply(levels(case_missense3$RAAS_quartile), function(quartile) {
  quartile_data <- case_missense3 %>% filter(RAAS_quartile == quartile)
  compute_ratio_and_ci(quartile_data)
})

# Step 1: Calculate the average ratio and CI for each quartile for both 'low AM' and 'high AM'
comparison_results <- lapply(levels(case_missense3$RAAS_quartile), function(quartile) {
  quartile_data <- case_missense3 %>% filter(RAAS_quartile == quartile)
  
  avg_ratio_low <- sum(quartile_data$obs_low, na.rm = TRUE) / sum(quartile_data$exp_low, na.rm = TRUE)
  ci_low <- generate_ci(sum(quartile_data$obs_low), sum(quartile_data$exp_low))
  
  avg_ratio_high <- sum(quartile_data$obs_high, na.rm = TRUE) / sum(quartile_data$exp_high, na.rm = TRUE)
  ci_high <- generate_ci(sum(quartile_data$obs_high), sum(quartile_data$exp_high))
  
  result <- data.frame(
    RAAS_quartile = rep(quartile, 2),
    Ratio = c(avg_ratio_low, avg_ratio_high),
    LowerCI = c(ci_low$lower, ci_high$lower),
    UpperCI = c(ci_low$upper, ci_high$upper),
    DataType = factor(c("Low AM missense in corresp. genes", "High AM missense in corresp. genes"), 
                      levels = c("Missense", "Low AM missense in corresp. genes", "High AM missense in corresp. genes"))
  )
  
  return(result)
})

# Combine the results into a single data frame
comparison_plot_data <- do.call(rbind, comparison_results)

# Prepare the Missense Data (assuming 'quartile_results' is already calculated similarly to 'comparison_results')
quartile_plot_data <- data.frame(
  RAAS_quartile = factor(levels(case_missense3$RAAS_quartile), 
                         levels = rev(levels(case_missense3$RAAS_quartile))),
  Ratio = unlist(unname(lapply(quartile_results, function(res) res$ratio))),
  LowerCI = unlist(unname(lapply(quartile_results, function(res) res$ci[1]))),
  UpperCI = unlist(unname(lapply(quartile_results, function(res) res$ci[2]))),
  DataType = factor("Missense", levels = c("Missense", "Low AM missense in corresp. genes", "High AM missense in corresp. genes"))
)

# Combine all the data into one data frame for plotting
full_plot_data <- rbind(quartile_plot_data, comparison_plot_data)





# Generate the color palette using arno(4)
mcol <- viridis::inferno(5)
vcol <- viridis::viridis(5)
mcol[5] <- vcol[5]
arno <- colorRampPalette(mcol)

# Generate the colors and reverse them to match the order of the quartiles
quartile_colors <- rev(arno(4))

# Define the custom color mapping for AlphaMissense categories
alpha_missense_colors <- c("Low AM missense in corresp. genes" = "blue", 
                           "High AM missense in corresp. genes" = "red")


if(plot_constraint)
{
  ggplot(full_plot_data, aes(x = RAAS_quartile, y = Ratio)) +
    # RAAS quartile points with matching error bars
    geom_point(data = full_plot_data %>% filter(DataType == "Missense"), 
               aes(fill = RAAS_quartile), size = 5, shape = 21, color = "black") +  # RAAS_quartile points with black outline
    
    # Add matching error bars using the RAAS_quartile color
    geom_errorbar(data = full_plot_data %>% filter(DataType == "Missense"), 
                  aes(ymin = LowerCI, ymax = UpperCI, color = RAAS_quartile), width = 0.2) +  # Error bars matching colors
    
    # Add points for Low and High AM missense with custom colors
    geom_point(data = full_plot_data %>% filter(DataType != "Missense"), 
               aes(color = DataType), size = 4) +
    
    # Map the RAAS_quartile to fill colors but reverse the labels
    scale_fill_manual(values = quartile_colors, 
                      labels = c(c("Highest RAAS Quartile",
                                   "Third RAAS Quartile", 
                                   "Second RAAS Quartile", 
                                   "Lowest RAAS Quartile")),
                      limits = quartile_plot_data$RAAS_quartile) +  # Corresponds to the lowest value range
    
    # Separate color scale for AlphaMissense categories, apply custom labels with title "Controls"
    scale_color_manual(values = alpha_missense_colors,
                       labels = c("Low AM missense in corresp. genes", 
                                  "High AM missense in corresp. genes"),
                       name = "Controls") +  # Add "Controls" title for AM missense
    
    theme_minimal() +
    
    # Remove the title of the legend for RAAS Quartiles
    labs(title = "Missense AAS vs low and high AM missense, broken by RAAS quartile",
         x = "RAAS Quartile", y = "Ratio (Obs/Exp)", fill = NULL) +
    
    # Adjust text size for all elements
    theme(legend.position = "right",
          legend.text = element_text(size = 14),  # Increase the legend text size
          legend.title = element_text(size = 16), # Increase legend title size
          text = element_text(size = 16),         # Increase overall text size (axis labels, title)
          axis.text = element_text(size = 14),    # Increase axis tick text size
          axis.title = element_text(size = 16),   # Increase axis title size
          plot.title = element_text(size = 18))   # Increase the plot title size
}

