##### READ ME #####

# This code, together with the count data for calcareous nannofossils from Buttenheim
# is part of the master's thesis "Response of calcareous nannofossils to environmental
# changes of the Pliensbachian-Toarcian boundary (Early Jurassic) in the Buttenheim 
# clay pit (Bavaria, Germany)". 

##### Loading of packages #####

library(vegan)
library(ggplot2)
library(patchwork)
library(dplyr)

# Set seed for reproducibility
set.seed(2024)

##### Loading of data #####

# Load quantitative data of calcareous nannofossils from Buttenheim 
count_dat <- read.csv("data.csv", sep=";")
# Set the row names to the sample designations
rownames(count_dat) <- count_dat$sample

# Check for missing values in count_dat
if (any(is.na(count_dat))) {
  cat("There are missing values in count_dat!\nIf a genus is not identified in a sample, it should be recorded as 0!")
}

# Load geochemical data, kindly provided  A. Merkel (unpublished)
geochem_dat <- read.csv("geochemical.csv", sep=";")

# Check for missing values in geochem_dat
if (any(is.na(geochem_dat))) {
  cat("There are missing values in geochem_dat!")
}

##### Calculate absolute and relative abundance of calcareous nannofossils #####

# Standardize the absolute abundance to 100 FOV
count_dat$abundance <- (count_dat$total_nanno/count_dat$total_FOV)*100

# Create list of genera names
genera <- c("Schizo", "Crepido", "Bisc", "Sim", "Solla", "Stau", "Strad", "Loth")

# Create list of genera percentage 
genera_per <- c("Schizo_per", "Crepido_per", "Bisc_per", "Sim_per", 
                "Solla_per", "Stau_per", "Strad_per", "Loth_per")

# Calculate the relative abundance of genera and save it to the count data frame
for (i in seq_along(genera)) {
  col_name <- paste0(genera[i], "_per")
  count_dat[[col_name]] <- (count_dat[[genera[i]]] / count_dat$total_nanno) * 100
}

##### Estimating error using standard error of the proportion ######

# Subset to data frame, that includes the abundance counts of genera and total nannofossils
genera_percent <- count_dat[, c("total_nanno", genera_per)]
 
# Create a empty container for standard errors
se_list <- list()

# Function to calculate standard errors of proportions
calculate_se <- function(percentages, total_n) {
  proportions <- percentages / 100  # Convert percentages to proportions
  se <- sqrt((proportions * (1 - proportions)) / total_n)
  # Convert the standard error from proportions to percentage
  se <- se * 100
  return(se)
}

# Calculate standard error of proportions for each genus and save it into cound_dat
for (i in seq_along(genera_per)) {
  genus <- genera_per[i]
  se_col_name <- paste0(genera[i], "_se_prop")
  count_dat[[se_col_name]] <- calculate_se(genera_percent[[genus]], genera_percent$total_nanno)
}

##### Error estimation using bootstrapping #####

# Define a matrix to save errors for each genus
boot_se_list <- list()

# Define a function to repeat category entries a given number of times
repeach <- function(category, counts) {
  if (length(category) != length(counts)) stop("The length of counts and category don't match!")
  final <- unlist(mapply(rep, category, counts))
  return(final)
}

# Define the number of bootstrap trials
trials <- 10000

# Bootstrapping for each sample, that saves the results into count_dat
for (k in 1:nrow(count_dat)) {
  counts <- count_dat[k, genera]
  
  # Create a long format of data
  long <- repeach(category = genera, counts = as.numeric(counts))
  
  # Initialize results matrix for this sample
  results <- matrix(0, nrow = trials, ncol = length(genera))
  colnames(results) <- genera
  
  for (i in 1:trials) {
    # Create replicate sample
    replicate <- sample(long, length(long), replace = TRUE)
    
    # Tabulate the bootstrapped counts
    bootTable <- table(replicate)
    
    # Store the results of the bootstrapped counts
    results[i, names(bootTable)] <- bootTable
  }
  
  # Calculate standard error of abundances for this sample
  se <- apply(results, 2, sd)
  
  # Store the standard errors
  for (i in seq_along(genera)) {
    se_col_name <- paste0(genera[i], "_boot_se")
    count_dat[k, se_col_name] <- se[i]
  }
}

##### Compare results from Proportions and Bootstrapping #####

# Create a summary data frame to compare Proportions and Bootstrapping
comparison_df <- data.frame(
  Genera = rep(genera, each = nrow(count_dat)),
  SE_Prop = unlist(lapply(genera, function(g) count_dat[[paste0(g, "_se_prop")]])),
  SE_Boot = unlist(lapply(genera, function(g) count_dat[[paste0(g, "_boot_se")]]))
)

# Calculate differences
comparison_df$Difference <- comparison_df$SE_Boot - comparison_df$SE_Prop

# Perform t-Test
t_test_result <- t.test(comparison_df$SE_Boot, comparison_df$SE_Prop, paired = TRUE)

# Bootstrapping captures more variability in the SE, so it is used to calculate 
# the lower and upper confidence intervals 

# Set the Z-score for a 95% confidence interval
Z <- 1.96

# Loop over each genus to calculate lower and upper confidence intervals
for (i in seq_along(genera_per)) {
  genus <- genera_per[i]
  se_col_name <- paste0(genera[i], "_boot_se")
  
  # Calculate lower and upper confidence intervals
  lower_ci_col_name <- paste0(genera[i], "_lower_boot")
  upper_ci_col_name <- paste0(genera[i], "_upper_boot")
  
  count_dat[[lower_ci_col_name]] <- count_dat[[genus]] - (Z * count_dat[[se_col_name]])
  count_dat[[upper_ci_col_name]] <- count_dat[[genus]] + (Z * count_dat[[se_col_name]])
}

##### Calculating Shannon Index and richness #####

# Subset to data frame, that only includes the abundance counts of genera
genera_count <- count_dat[ ,genera]

# Function for calculating the Shannon Index from the vegan package
shannon_index <- vegan::diversity(genera_count, index="shannon")
# Add the Shannon Index as a new column to the data frame
count_dat$Shannon <- shannon_index

# Calculate richness as the number of genera with counts > 0 in each sample (row)
richness <- rowSums(genera_count > 0)
# Add the richness as a new column to the data frame
count_dat$Richness <- richness

##### Subsetting of data for plotting and clustering #####

# Subset to exclude duplicate entries from the plotting
unique_dat <- count_dat[count_dat$duplicate != "yes", ]

# Subset to all soft lithologies
dat_soft <- unique_dat[unique_dat$lithology == 0,]
#Subset to all hard lithologies
dat_hard <- unique_dat[unique_dat$lithology == 1,]

# Subset to Toarcian data with more than 300 total counts
Toarc <- unique_dat[unique_dat$total_nanno > 300,]
# Subset to soft lithlogies in the Toarcian
Toarc_soft <- Toarc[Toarc$lithology == 0,]
# Subset to hard lithologies in the Toarcian
Toarc_hard <- Toarc[Toarc$lithology == 1,]

##### Creating plots for geochemical data, absolute and relative abundance #####

# creating a theme for all the plots of the genera
theme_plot <- ggplot2::theme(
  axis.text.y = element_blank(), 
  axis.ticks.y = element_blank(),
  axis.text.x = element_text(size = 10),
  axis.title.x = element_text(size = 12),
)

# Geochemical data: Carbon 13 isotops
plot_C <- ggplot2::ggplot(geochem_dat, aes(x= X13C, y=height)) +
  geom_point(color="black", size=2) +
  geom_line(color="black", orientation="y") +
  ylim(0, max(Toarc$height)) +
  labs(x="δ13C isotopes", y = "") +
  theme_plot

# Geochemical data: Oxygen 18 isotops
plot_O <- ggplot2::ggplot(geochem_dat, aes(x= X18O, y=height)) +
  geom_point(color="black", size=2) +
  geom_line(color="black", orientation="y") +
  ylim(0, max(Toarc$height)) +
  labs(x="δ18O isotopes", y = "") + 
  theme_plot

# Abundance per 100 FOVs for complete outcrop
plot_abundance <- ggplot(dat_soft, aes(x= abundance, y=height)) +
  geom_point(data=dat_hard, aes(x=abundance, y=height), color="red", size=2) +
  geom_line(data=dat_hard, aes(x=abundance, y=height), color="red", orientation="y") +
  geom_point(color="black", size=2) +
  geom_line(color="black", orientation="y") +
  ylim(0, max(Toarc$height)) +
  scale_x_continuous(limits = c(0, 210)) +
  labs(x="Absolute abundance", y = "") + 
  theme_plot

# Shannon Index
plot_Shannon <- ggplot(Toarc_soft, aes(x= Shannon, y=height)) +
  geom_point(data=Toarc_hard, aes(x=Shannon, y=height), color="red", size=2) +
  geom_line(data=Toarc_hard, aes(x=Shannon, y=height), color="red", orientation="y") +
  geom_point(color="black", size=2) +
  geom_line(color="black", orientation="y") +
  scale_x_continuous(limits = c(1, 1.9)) +
  ylim(0, max(Toarc$height)) +
  labs(x="Shannon Index", y = "") + 
  theme_plot

# Richness
plot_Richness <- ggplot(dat_soft, aes(x= Richness, y=height)) +
  geom_point(data=dat_hard, aes(x=Richness, y=height), color="red", size=2) +
  geom_line(data=dat_hard, aes(x=Richness, y=height), color="red", orientation="y") +
  geom_point(color="black", size=2) +
  geom_line(color="black", orientation="y") +
  scale_x_continuous(limits = c(0, 8)) +
  ylim(0, max(Toarc$height)) +
  labs(x="Richness", y = "") + 
  theme_plot

# Combined plot for the geochemical data, the absolute abundance, Shannon Index and Richness
combined_geo_abund_plot <- plot_C + plot_O + plot_abundance + plot_Richness + plot_Shannon +
  plot_layout(ncol = 5, widths = rep(1, 5), guides = "collect") & 
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

# Print combined plot
windows()
print(combined_geo_abund_plot)


# Create plot for absolute and relative abundance in the Toarcian

# Abundance per 100 FOVs for Toarcian
plot_abundance_Toarc <- ggplot(Toarc_soft, aes(x= abundance, y=height)) +
  geom_point(data=Toarc_hard, aes(x=abundance, y=height), color="red") +
  geom_line(data=Toarc_hard, aes(x=abundance, y=height), color="red", orientation="y") +
  geom_point(color="black") +
  geom_line(color="black", orientation="y") +
  scale_x_continuous(limits = c(0, 210)) +
  ylim(2, max(Toarc$height)) +
  labs(x="Absolute abundance", y = "")


# Create function for plotting the relative abundance of genera
plot_genera <- function(data_soft, data_hard, x_var, x_label, low_var, up_var, x_lim) {
  ggplot(data_soft, aes_string(x = x_var, y = "height")) +
    geom_ribbon(aes_string(xmin = low_var, xmax = up_var), alpha = 0.2) + 
    geom_point(data = data_hard, aes_string(x = x_var, y = "height"), color = "red") +
    geom_line(data = data_hard, aes_string(x = x_var, y = "height"), color = "red", orientation = "y") +
    geom_ribbon(data = data_hard, aes_string(xmin = low_var, xmax = up_var), alpha = 0.1) + 
    geom_point(color = "black") +
    geom_line(color = "black", orientation = "y") +
    scale_x_continuous(limits = c(-6, x_lim)) +
    ylim(2, max(Toarc$height)) +
    labs(x = x_label, y = "") + 
    theme_plot
}

# Create plots for the different genera
plot_Schizo <- plot_genera(Toarc_soft, Toarc_hard, "Schizo_per", "Schizosphaerella sp.", "Schizo_lower_boot", "Schizo_upper_boot", 90)
plot_Crepido <- plot_genera(Toarc_soft, Toarc_hard, "Crepido_per", "Crepidolithus sp.", "Crepido_lower_boot", "Crepido_upper_boot", 35)
plot_Stau <- plot_genera(Toarc_soft, Toarc_hard, "Stau_per", "Staurolithes sp.", "Stau_lower_boot", "Stau_upper_boot", 35)
plot_Bisc <- plot_genera(Toarc_soft, Toarc_hard, "Bisc_per", "Biscutum sp.", "Bisc_lower_boot", "Bisc_upper_boot", 35)
plot_Sim <- plot_genera(Toarc_soft, Toarc_hard, "Sim_per", "Similiscutum sp.", "Sim_lower_boot", "Sim_upper_boot", 35)
plot_Solla <- plot_genera(Toarc_soft, Toarc_hard, "Solla_per", "Sollasites sp.", "Solla_lower_boot", "Solla_upper_boot", 35)
plot_Strad <- plot_genera(Toarc_soft, Toarc_hard, "Strad_per", "Stradnerlithus sp.", "Strad_lower_boot", "Strad_upper_boot", 35)
plot_Loth <- plot_genera(Toarc_soft, Toarc_hard, "Loth_per", "Lotharingus sp.", "Loth_lower_boot", "Loth_upper_boot", 35)

# Combined plot for the absolute and all the relative abundances
combined_abundance_plot <- plot_abundance_Toarc + plot_Schizo + plot_Crepido + plot_Stau + plot_Bisc + plot_Sim + plot_Strad + plot_Loth + plot_Solla +
  plot_layout(nrow = 1, widths = rep(1, 9), guides = "collect") & theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

# Print  plots in separate window
windows()
print(combined_abundance_plot)


##### Hierarchical clustering using the percentage abundance #####

# Isolate the relevant columns
clust_data <- Toarc[, c("Schizo_per", "Crepido_per", "Stau_per", "Bisc_per", 
                        "Sim_per", "Solla_per", "Strad_per", "Loth_per")]

# Extract height for labeling later
height_labels <- Toarc$height

# Calculate the distance matrix using Euclidean distance
eu_dis <- dist(clust_data)
# Euclidean distance with different linkage methods
eu_clus_sin <- hclust(eu_dis, method = "single")
eu_clus_com <- hclust(eu_dis, method = "complete")
eu_clus_war <- hclust(eu_dis, method = "ward.D")
eu_clus_ave <- hclust(eu_dis, method = "average")

# Calculate the distance matrix using Bray-Curtis distance
br_dis <- vegan::vegdist(clust_data, method = "bray")
# Bray-Curtis distance with different linkage methods
br_clus_sin <- hclust(br_dis, method = "single")
br_clus_com <- hclust(br_dis, method = "complete")
br_clus_war <- hclust(br_dis, method = "ward.D")
br_clus_ave <- hclust(br_dis, method = "average")

# Function for plotting the dendrograms 
fun_dend <- function(clus_method, title){
  par(mar = c(2, 6, 2, 1) + 0.1)
  plot(clus_method, main = "", sub = "", xlab = "", ylab = "", lwd = 2, cex = 1.5, cex.axis = 1.5, labels = round(height_labels[clus_method$order], 2))
  mtext(title, side = 2, line = 4.5, cex = 0.8)
}

windows()
# Set up layout for dendrogram plots
par(mfrow = c(4, 2), mar = c(4, 6, 2, 1))  # 4 rows and 2 columns

# Create individual dendrogram plots
fun_dend(eu_clus_sin, "Euclidean - Single Linkage")
fun_dend(eu_clus_com, "Euclidean - Complete Linkage")
fun_dend(eu_clus_war, "Euclidean - Ward's Method")
fun_dend(eu_clus_ave, "Euclidean - Average Linkage")
fun_dend(br_clus_sin, "Bray-Curtis - Single Linkage")
fun_dend(br_clus_com, "Bray-Curtis - Complete Linkage")
fun_dend(br_clus_war, "Bray-Curtis - Ward's Method")
fun_dend(br_clus_ave, "Bray-Curtis - Average Linkage")

##### NMDS by percentage of nannofossils #####

#create vector with full genus names
genus_fullnames <- c(
  Schizo_per = "Schizosphaerella", 
  Crepido_per = "Crepidolithus", 
  Sim_per = "Similiscutum", 
  Bisc_per = "Biscutum", 
  Stau_per = "Staurolithes", 
  Solla_per = "Sollasites", 
  Loth_per = "Lotharingius", 
  Strad_per = "Stradnerlithus"
)

# Set the row names
rownames(Toarc) <- Toarc$height

# Subset to the counts of Toarcian samples with more than 300 specimen
Toarc_per <- Toarc[ ,genera_per]

#Perform NMDS analysis 
NMDS_per=metaMDS(Toarc_per, # data
                 k=2, # number of dimensions
                 distance = "bray") # distance matrix

# Extract site and species scores
nmds_sites <- as.data.frame(scores(NMDS_per, display = "sites"))
nmds_species <- as.data.frame(scores(NMDS_per, display = "species"))

# Rename the rownames of NMDS genera
rownames(nmds_species) <- genus_fullnames[rownames(nmds_species)]

# Add height and lithology to the NMDS site scores
nmds_sites$height <- as.numeric(rownames(Toarc))
nmds_sites$lithology <- Toarc$lithology

# Define colors based on lithology
lithology_colors <- ifelse(nmds_sites$lithology == 0, "green", "orange")

# Create NMDS plot using ggplot
plot_nmds <- ggplot() +
  # Plot sites as points colored by lithology
  geom_point(data = nmds_sites, aes(x = NMDS1, y = NMDS2, color = factor(lithology)), size = 5) +
  # Add arrows between successive samples (indicating progression by height)
  geom_segment(data = nmds_sites %>% arrange(height),
               aes(x = lag(NMDS1), y = lag(NMDS2), xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 1.5) +
  # Add height labels
  geom_text(data = nmds_sites, aes(x = NMDS1, y = NMDS2, label = height), 
            vjust = -1, size = 6, color = "black") +  
  # Customize plot (before adding species text to ensure it's in the foreground)
  scale_color_manual(values = c("green", "orange"), name = "Lithology") +
  labs(x = "NMDS1", y = "NMDS2", title = "NMDS Plot", size = 3) +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 20), 
        axis.text = element_text(size = 14)) +
  # Plot species as text in the foreground (dark red labels)
  geom_text(data = nmds_species, aes(x = NMDS1, y = NMDS2, label = rownames(nmds_species)),
            color = "darkred", size = 5)


### Isolate NMDS1
nmds1 <- scores(NMDS_per, display = "sites")[, 1]
nmds_plot_data <- data.frame(
  height = as.numeric(rownames(Toarc)),
  NMDS1 = nmds1
)

# Determine the common y-axis limits based on both NMDS1 and NMDS2
y_limits <- range(c(nmds1, scores(NMDS_per, display = "sites")[, 2]))

plot_nmds1 <- ggplot(nmds_plot_data, aes(x = height, y = NMDS1)) +
  geom_point(color = "black") +
  geom_line(color = "black") +
  labs(x = "Height (Sequence of Sample Heights)", y = "NMDS1") +
  theme(aspect.ratio = 1) +  # Fix the aspect ratio to 1
  ylim(y_limits)  # Set common y-axis limits

### Isolate NMDS2
nmds2 <- scores(NMDS_per, display = "sites")[, 2]
nmds_plot_data <- data.frame(
  height = as.numeric(rownames(Toarc)),
  NMDS2 = nmds2
)

plot_nmds2 <- ggplot(nmds_plot_data, aes(x = height, y = NMDS2)) +
  geom_point(color = "black") +
  geom_line(color = "black") +
  labs(x = "Height (Sequence of Sample Heights)", y = "NMDS2") +
  theme(aspect.ratio = 1) +  # Fix the aspect ratio to 1
  ylim(y_limits)  # Set common y-axis limits


# Combine plots side by side using 
combined_nmds <- (plot_nmds / (plot_nmds1 + plot_nmds2)) +
  plot_layout(nrow = 2)

# Display the combined plot
windows()
print(combined_nmds)

##### Compare the absolute and relative abundance across both preparation methods #####

# Define colors for each genera
genera_colors <- c(
  "Schizo" = "red",
  "Crepido" = "blue",
  "Bisc" = "green",
  "Sim" = "purple",
  "Solla" = "orange",
  "Stau" = "cyan",
  "Strad" = "magenta",
  "Loth" = "yellow"
)

# Define colors for each genera
genera_colors_per <- c(
  "Schizo_per" = "red",
  "Crepido_per" = "blue",
  "Bisc_per" = "green",
  "Sim_per" = "purple",
  "Solla_per" = "orange",
  "Stau_per" = "cyan",
  "Strad_per" = "magenta",
  "Loth_per" = "yellow"
)

### Absolute abundance for limestone

# Filter the data to include only rows where the sample is "13A" or "13C"
prep_13 <- count_dat[count_dat$sample %in% c("13A", "13C"), ]

# Subset to include the absolute abundance of genera
prep_13_gen <- prep_13[, c("Schizo", "Crepido", "Bisc", "Sim",
                           "Solla", "Stau", "Strad", "Loth")]

# Subset to count data of the samples
prep_A13_gen <- as.numeric(prep_13_gen[1, ])
prep_C13_gen <- as.numeric(prep_13_gen[2, ])

# Perform Spearman correlation
spearman_test_13_gen <- cor.test(prep_A13_gen, prep_C13_gen, method = "spearman")

# perfrom non-parametric wilcox test
wilcox_prep_13_gen <- wilcox.test(prep_A13_gen, prep_C13_gen)

# warp dataframe for plotting
prep_13_gen <- as.data.frame(t(prep_13_gen))

# Add a genera column
prep_13_gen$Genera <- rownames(prep_13_gen)
# Convert Genera to a factor
prep_13_gen$Genera <- factor(prep_13_gen$Genera, levels = names(genera_colors))

# Create the plot
plot_prep_13_gen <- ggplot(prep_13_gen, aes(x = `13C`, y = `13A`, label = Genera, color = Genera)) + 
  geom_abline(intercept = 0, slope =1, col="white", lwd=1) +  
  geom_point(size = 3, show.legend = FALSE) +
  labs(x = "Total counts in fractured sample", y = "Total counts in ground and etched sample", title = "Preparation methods of limestone") + 
  xlim(c(0, 250)) +
  ylim(c(0, 250)) +
  theme_gray() +
  theme(aspect.ratio = 1)

### Relative abundance for limestone

# Subset to include the relative abundance of genera
prep_13_per <- prep_13[, c("Schizo_per", "Crepido_per", "Bisc_per", "Sim_per",
                           "Solla_per", "Stau_per", "Strad_per", "Loth_per")]

# Subset to count data of the samples
prep_A13_per <- as.numeric(prep_13_per[1, ])
prep_C13_per <- as.numeric(prep_13_per[2, ])
 
# Perform Spearman correlation
spearman_test_13_per <- cor.test(prep_A13_per, prep_C13_per, method = "spearman")

# perfrom non-parametric wilcox test
wilcox_prep_13_per <- wilcox.test(prep_A13_per, prep_C13_per)

# warp dataframe for plotting
prep_13_per <- as.data.frame(t(prep_13_per))

# Add a genera column
prep_13_per$Genera <- rownames(prep_13_per)
# Convert Genera to a factor
prep_13_per$Genera <- factor(prep_13_per$Genera, levels = names(genera_colors_per))

# Create the plot
plot_prep_13_per <- ggplot(prep_13_per, aes(x = `13C`, y = `13A`, label = Genera, color = Genera))+ 
  geom_abline(intercept = 0, slope =1, col="white", lwd=1) +  
  geom_point(size = 3, show.legend = FALSE) +
  labs(x = "Relative abundance in fractured sample", y = "Relative abundance in ground and etched sample", title = "Preparation methods of limestone") + 
  xlim(c(0, 100)) +
  ylim(c(0, 100)) +
  theme_gray() +
  theme(aspect.ratio = 1)

### Absolute abundance for marl

# Filter the data to include only rows where the sample is "14A" or "14C"
prep_14 <- count_dat[count_dat$sample %in% c("14A", "14C"), ]

# Subset to include the absolute abundance of genera
prep_14_gen <- prep_14[, c("Schizo", "Crepido", "Bisc", "Sim",
                           "Solla", "Stau", "Strad", "Loth")]

# Subset to count data of the samples
prep_A14_gen <- as.numeric(prep_14_gen[1, ])
prep_C14_gen <- as.numeric(prep_14_gen[2, ])

# Perform Spearman correlation
spearman_test_14_gen <- cor.test(prep_A14_gen, prep_C14_gen, method = "spearman")

# perfrom non-parametric wilcox test
wilcox_prep_14_gen <- wilcox.test(prep_A14_gen, prep_C14_gen)

# warp dataframe for plotting
prep_14_gen <- as.data.frame(t(prep_14_gen))

# Add a genera column
prep_14_gen$Genera <- rownames(prep_14_gen)
# Convert Genera to a factor
prep_14_gen$Genera <- factor(prep_14_gen$Genera, levels = names(genera_colors))

# Create the plot
plot_prep_14_gen <- ggplot(prep_14_gen, aes(x = `14C`, y = `14A`, label = Genera, color = Genera)) + 
  geom_abline(intercept = 0, slope =1, col="white", lwd=1) +  
  geom_point(size = 3, show.legend = TRUE) +
  labs(x = "Total counts in fractured sample", y = "Total counts in ground and etched sample", title = "Preparation methods of marl") + 
  xlim(c(0, 350)) +
  ylim(c(0, 350)) +
  theme_gray() +
  theme(aspect.ratio = 1, legend.position = "bottom")
  
### Relative Abundance for marl

# Subset to include the relative abundance of genera
prep_14_per <- prep_14[, c("Schizo_per", "Crepido_per", "Bisc_per", "Sim_per",
                           "Solla_per", "Stau_per", "Strad_per", "Loth_per")]

# Subset to count data of the samples
prep_A14_per <- as.numeric(prep_14_per[1, ])
prep_C14_per <- as.numeric(prep_14_per[2, ])

# Perform Spearman correlation
spearman_test_14_per <- cor.test(prep_A14_per, prep_C14_per, method = "spearman")

# perfrom non-parametric wilcox test
wilcox_prep_14_per <- wilcox.test(prep_A14_per, prep_C14_per)

# warp dataframe for plotting
prep_14_per <- as.data.frame(t(prep_14_per))

# Add a genera column
prep_14_per$Genera <- rownames(prep_14_per)
# Convert Genera to a factor
prep_14_per$Genera <- factor(prep_14_per$Genera, levels = names(genera_colors_per))

# Create the plot
plot_prep_14_per <- ggplot(prep_14_per, aes(x = `14C`, y = `14A`, label = Genera, color = Genera)) + 
  geom_abline(intercept = 0, slope =1, col="white", lwd=1) +  
  geom_point(size = 3, show.legend = FALSE) +
  labs(x = "Relative abundance in fractured sample", y = "Relative abundance in ground and etched sample", title = "Preparation methods of marl") + 
  xlim(c(0, 100)) +
  ylim(c(0, 100)) +
  theme_gray() +
  theme(aspect.ratio = 1)


# Combine the plots
compare_prep_plot <- plot_prep_13_gen +  plot_prep_13_per + plot_prep_14_gen + plot_prep_14_per
  plot_layout(nrow = 1, widths = rep(2, 2), guides = "collect") & theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

windows()
compare_prep_plot





##### Assess the Reproducibility of the counting method #####


### Absolute abundance of limestone

# Filter the data to include only rows where the sample is "13C" or "13D"
dup_13 <- count_dat[count_dat$sample %in% c("13C", "13D"), ]

# Subset to include the total abundance and count data of the genera 
dup_13_gen <- dup_13[, c("Schizo", "Crepido", "Bisc", "Sim",
                         "Solla", "Stau", "Strad", "Loth")]

# perfrom non-parametric wilcox test
dup_C13 <- as.numeric(dup_13_gen[1, ])
dup_D13 <- as.numeric(dup_13_gen[2, ])
wilcox_13 <- wilcox.test(dup_C13, dup_D13)

# warp dataframe for plotting
dup_13_gen <- as.data.frame(t(dup_13_gen))

# Add a genera column
dup_13_gen$Genera <- rownames(dup_13_gen)
# Convert Genera to a factor
dup_13_gen$Genera <- factor(dup_13_gen$Genera, levels = names(genera_colors))

# perfom correlation test
cor.test(dup_13_gen[, 1], dup_13_gen[, 2], method = "spearman")

# perform linear regression
dup_13_reg <-lm(formula = `13D` ~ `13C`, data=dup_13_gen)
# Save coefficients of linear regression
dup_13_coeff <- coefficients(dup_13_reg)
dup_13_intercept<-dup_13_coeff[1] 
dup_13_slope<- dup_13_coeff[2]

# Create the plot
plot_dup_13 <- ggplot(dup_13_gen, aes(x = `13C`, y = `13D`, label = Genera, color = Genera)) + 
  geom_abline(intercept = 0, slope =1, col="white", lwd=1) +
  geom_point(size = 3, show.legend = TRUE) +
  geom_abline(intercept = dup_13_intercept, slope = dup_13_slope, col = "black", linetype = "dashed") +
  labs(x = "Total counts in sample 13C", y = "Total counts in sample 13D", title = "Duplicate counts of limestone") + 
  xlim(c(-20, 250)) +
  ylim(c(-20, 250)) +
  theme_gray() +
  theme(aspect.ratio = 1, legend.position = "bottom") +  # Fix the aspect ratio to 1
  annotate("text", x = 150, y = 35, label = paste("Intercept:", format(dup_13_intercept, digits = 2, nsmall = 1)), hjust = 0, color = "black") +
  annotate("text", x = 150, y = 15, label = paste("Slope:", format(dup_13_slope, digits = 2, nsmall = 2)), hjust = 0, color = "black")

### Absolute abundance of marl 

# Filter the data to include only rows where the sample is "14C" or "14D"
dup_14 <- count_dat[count_dat$sample %in% c("14C", "14D"), ]

# Subset to genus_per
# Subset columns directly by their names
dup_14_gen <- dup_14[, c("Schizo", "Crepido", "Bisc", "Sim",
                         "Solla", "Stau", "Strad", "Loth")]

# Perfrom non-parametric wilcox.test
dup_C14 <- as.numeric(dup_14_gen[1, ])
dup_D14 <- as.numeric(dup_14_gen[2, ])
wilcox_14 <- wilcox.test(dup_C14, dup_D14)

# warp dataframe for plotting
dup_14_gen <- as.data.frame(t(dup_14_gen))

# Add a genera column
dup_14_gen$Genera <- rownames(dup_14_gen)
# Convert Genera to a factor
dup_14_gen$Genera <- factor(dup_14_gen$Genera, levels = names(genera_colors))

# perfom correlation test
cor.test(dup_14_gen[, 1], dup_14_gen[, 2], method = "spearman")

# perform linear regression
dup_14_reg <-lm(formula = `14D` ~ `14C`, data=dup_14_gen)
# Save coefficients of linear regression
dup_14_coeff <- coefficients(dup_14_reg)
dup_14_intercept<-dup_14_coeff[1] 
dup_14_slope<- dup_14_coeff[2]

# Create the plot
plot_dup_14 <- ggplot(dup_14_gen, aes(x = `14C`, y = `14D`, label = Genera, color = Genera)) + 
  geom_abline(intercept = 0, slope =1, col="white", lwd=1) +  
  geom_point(size = 3, show.legend = FALSE) +
  geom_abline(intercept = dup_14_intercept, slope = dup_14_slope, col= "black", linetype="dashed") +
  labs(x = "Total counts in sample 14C", y = "Total counts in sample 14D", title = "Duplicate counts of marl") + 
  xlim(c(0, 400)) +
  ylim(c(0, 400)) +
  theme_gray() +
  theme(aspect.ratio = 1) + # Fix the aspect ratio to 1
  annotate("text", x = 250, y = 75, label = paste("Intercept:", format(dup_14_intercept, digits = 3, nsmall = 1)), hjust = 0, color = "black") +
  annotate("text", x = 250, y = 45, label = paste("Slope:", format(dup_14_slope, digits = 2, nsmall = 2)), hjust = 0, color = "black")

# Combined plot
compare_plot <- plot_dup_13 +  plot_dup_14 
  plot_layout(nrow = 1, widths = rep(1, 2), guides = "collect") & theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

windows()
compare_plot
