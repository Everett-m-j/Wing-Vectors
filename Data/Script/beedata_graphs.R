#beedata graphing
#Everett Mathiason
#08/04/2025

rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggcorrplot)
library(ggplot2)


#Get folder where data is
set_wd <- ("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data")
read.csv("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/occurrences.csv")
bee_occurences <- read.csv("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/occurrences.csv")
colnames(bee_occurences)
bee_occurences <- bee_occurences %>%
  filter(!is.na(percent_asymmetry))
bee_occurences <- bee_occurences %>%
  select(-Notes., -Total.Complete.Agr, -Total.Complete.Dev, -Total.Complete.For, -ImageName)
#bee_occurences$WingSide <- ifelse(bee_occurences$WingSide == "L", 1,
                                 #ifelse(bee_occurences$WingSide == "R", 2, NA))
#bee_occurences$Land_cover <- ifelse(bee_occurences$Land_cover == "Agricultural", 1,
                                  #ifelse(bee_occurences$Land_cover == "Developed", 2,
                                         #ifelse(bee_occurences$Land_cover == "Forest", 3, NA)))
bee_occurences_clean <- bee_occurences[-c(8, 7, 51, 52, 93, 94, 91, 92, 47, 48), ]

pca_data <- bee_occurences_clean %>%
  select(area, circularity, length, width, percent_asymmetry, percent_overlap) %>%
  na.omit()
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

summary(pca_result)         # To see variance explained
pca_result$rotation         # Loadings (how variables contribute)
pca_result$x                # The transformed PCA scores

# Create a data frame with PCA scores and Land Cover
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Land_cover <- bee_occurences_clean$Land_cover[complete.cases(pca_data)]

# Plot PC1 vs PC2 by land cover
ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(Land_cover))) +
  geom_point(size = 3) +
  labs(title = "PCA of Wing Metrics by Land Cover",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Land_Cover") +
  theme_minimal()
# Create convex hulls for each group
hulls <- pca_scores %>%
  group_by(Land_cover) %>%
  slice(chull(PC1, PC2))

ggplot(pca_scores, aes(x = PC1, y = PC2, color = Land_cover)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_polygon(data = hulls, aes(fill = Land_cover), alpha = 0.2, color = NA) +
  labs(title = "PCA with Convex Hulls by Land Cover",
       x = "PC1", y = "PC2") +
  theme_minimal()
# Extract and scale loadings
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$varname <- rownames(loadings)
loadings <- loadings %>%
  mutate(PC1 = PC1 * 3,  # scale arrows
         PC2 = PC2 * 3)

# Add to the plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Land_cover)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "t", level = 0.95, geom = "polygon", alpha = 0.2, aes(fill = Land_cover), color = NA) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = loadings, aes(x = PC1, y = PC2, label = varname), color = "black", vjust = 1.2) +
  labs(title = "PCA with Ellipses and Variable Loadings",
       x = "PC1", y = "PC2") +
  theme_minimal()

#linear graph overlap by land cover
# Calculate means
overlap_summary <- bee_occurences_clean %>%
  group_by(Land_cover) %>%
  summarise(mean_overlap = mean(percent_overlap, na.rm = TRUE))


# plot for anova
ggplot(bee_occurences_clean, aes(x = Land_cover, y = percent_asymmetry, color = Land_cover)) +
  geom_point() +
  labs(title = "Percent Asymmetry by Land Cover",
       x = "Land Cover",
       y = "Percent Asymmetry") +
  theme_minimal() +
  theme(legend.position = "none")

anova <- aov(percent_asymmetry ~ Land_cover, data = Data1)
summary(anova)
TukeyHSD(anova)

ggplot(Data1, aes(x = Land_cover, y = percent_asymmetry, color = Land_cover)) +
  geom_boxplot() +
  geom_jitter() +
  labs(title = "Percent Asymmetry by Land Cover",
       x = "Land Cover",
       y = "Percent Asymmetry") +
  theme_minimal() +
  theme(legend.position = "none")

Data1 <- subset(bee_occurences_clean, WingSide == "L")


# Measure percent overlap
data_high_overlap <- Data1 %>% filter(percent_overlap > 80)

ggplot(data_high_overlap, aes(x = Land_cover, y = percent_overlap, fill = Land_cover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(title = "High Overlap (>80%) by Land Cover",
       x = "Land Cover", y = "Percent Overlap") +
  theme_minimal()

#measure percent overlap by cell 1 or 4
ggplot(data_high_overlap, aes(x = factor(Cell), y = percent_overlap, fill = Land_cover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.4) +
  labs(
    title = "Percent Overlap by Wing Cell and Land Cover",
    x = "Wing Cell",
    y = "Percent Overlap"
  ) +
  theme_minimal()


#scatter plot length vs width
ggplot(bee_occurences_clean, aes(x = width, y = length, color = Land_cover)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(title = "Wing Cell Length vs. Wing Cell Width",
       x = "Cell Width", y = "Cell Length") +
  theme_minimal()

# box plot by area
ggplot(Data1, aes(x = Land_cover, y = area, color = Land_cover)) +
  geom_boxplot() +
  geom_jitter() +
  labs(title = "Percent Asymmetry by Land Cover",
       x = "Land Cover",
       y = "Percent Asymmetry") +
  theme_minimal() +
  theme(legend.position = "none")

#summarizing data
  data_summarized <- bee_occurences_clean %>%
    group_by(catalogNumber, WingSide) %>%
    summarise(across(c(area, circularity, length, width, percent_overlap), mean, na.rm = TRUE), .groups = "drop")
  
  # fluctuaing asymmetry
  data_wide <- bee_occurences_clean %>%
    pivot_wider(
      names_from = WingSide,
      values_from = c(area, circularity, length, width),
      names_glue = "{.value}_{WingSide}"
    )
  
  data_wide_1 <- data_wide %>%
    filter(Cell == "1")
  
  data_wide_2 <- data_wide %>% 
    filter(Cell == "4")
####################################################
  data_wide_FA <- data_wide %>% mutate(
    FA_area = abs(area_L - area_R),
    FA_length = abs(length_L - length_R),
    FA_width = abs(width_L - width_R),
    FA_circ = abs(circularity_L - circularity_R),
    FA_mean = rowMeans(select(., starts_with("FA_")), na.rm = TRUE)  # Optional: overall FA score
  )
  combined_data <- data_wide_FA %>%
    left_join(regression_results_all, 
              by = c("Land_cover", "Cell"))
  
  pca_data_FA <- data_wide_FA %>%
    select(FA_area, FA_circ, FA_length, FA_width) %>%
    na.omit()
  pca_result <- prcomp(pca_data_FA, center = TRUE, scale. = TRUE)
  
  summary(pca_result)         # To see variance explained
  pca_result$rotation         # Loadings (how variables contribute)
  pca_result$x                # The transformed PCA scores
  
  # Create a data frame with PCA scores and Land Cover
  pca_scores <- as.data.frame(pca_result$x)
  pca_scores$Land_cover <- data_wide$Land_cover[complete.cases(pca_data_FA)]
  
  # Plot PC1 vs PC2 by land cover
  ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(Land_cover))) +
    geom_point(size = 3) +
    labs(title = "PCA of Wing Metrics by Land Cover",
         x = "Principal Component 1",
         y = "Principal Component 2",
         color = "Land_Cover") +
    theme_minimal()
  
  # Plot PC1 vs nrmse results by land cover
  ggplot(pca_scores, aes(x = PC1, y = nrmse, color = as.factor(Land_cover))) +
    geom_point(size = 3) +
    labs(title = "PCA of Wing Metrics by Nrmse",
         x = "Principal Component 1",
         y = "NRMSE",
         color = "Land_Cover") +
    theme_minimal()
  
  # Extract and scale loadings
  loadings <- as.data.frame(pca_result$rotation[, 1:2])
  loadings$varname <- rownames(loadings)
  loadings <- loadings %>%
    mutate(PC1 = PC1 * 3,  # scale arrows
           PC2 = PC2 * 3)
  
  # Add to the plot
  ggplot(pca_scores, aes(x = PC1, y = PC2, color = Land_cover)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(type = "t", level = 0.95, geom = "polygon", alpha = 0.2, aes(fill = Land_cover), color = NA) +
    geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                 arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_text(data = loadings, aes(x = PC1, y = PC2, label = varname), color = "black", vjust = 1.2) +
    labs(title = "PCA with Ellipses and Variable Loadings",
         x = "PC1", y = "PC2") +
    theme_minimal()
  
  
  ggplot(data_wide, aes(x = area_L, y = area_R)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(title = "Left vs Right Area", x = "Left Area", y = "Right Area") +
    theme_minimal()
  
  
  # Bubble Plot
  # Select only the FA variables and Land_cover
  data_long_FA <- data_wide %>%
    select(Land_cover, 
           FA_area, FA_length, FA_width, FA_circ, FA_overlap) %>%
    pivot_longer(
      cols = starts_with("FA_"),
      names_to = "FA_variable",
      values_to = "FA_value"
    )
  
  # Summarize by mean FA per land cover × variable
  summary_data <- data_long_FA %>%
    group_by(Land_cover, FA_variable) %>%
    summarise(mean_FA = mean(FA_value, na.rm = TRUE),
              .groups = "drop")
  
 
   ggplot(summary_data, aes(x = Land_cover, y = FA_variable)) +
    geom_point(aes(size = mean_FA, fill = mean_FA), shape = 21, color = "black") +
    scale_size_continuous(range = c(3, 15)) +
    scale_fill_viridis_c(option = "plasma") +
    theme_minimal() +
    labs(
      x = "Land Cover",
      y = "Fluctuating Asymmetry Variable",
      size = "Mean FA",
      fill = "Mean FA",
      title = "Mean Fluctuating Asymmetry by Land Cover and Shape Variable"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  
   
   #plot by land cover
   ggplot(data_wide, aes(x = Land_cover, y = FA_circ, color = Land_cover)) +
     geom_point() +
     labs(title = "FA by Land Cover",
          x = "Land Cover",
          y = "Asymmetry") +
     theme_minimal() +
     theme(legend.position = "none")
   
   # Using NRMSE and Regression to find Asymmetry
   #NMSRE
   # Function to calculate NRMSE
   nrmse <- function(L, R) {
     model <- lm(R ~ L)
     rmse <- sqrt(mean(residuals(model)^2, na.rm = TRUE))
     mean_val <- mean(c(L, R), na.rm = TRUE)
     rmse / mean_val
   }
   
   
   # Calculate NRMSEs
   nrmse_area <- nrmse(data_wide$area_L, data_wide$area_R)
   nrmse_length <- nrmse(data_wide$length_L, data_wide$length_R)
   nrmse_width <- nrmse(data_wide$width_L, data_wide$width_R)
   nrmse_circularity <- nrmse(data_wide$circularity_L, data_wide$circularity_R)
   
   # regression
   model_area <- lm(area_R ~ area_L, data = data_wide)
   model_length <- lm(length_R ~ length_L, data = data_wide)
   model_width <- lm(width_R ~ width_L, data = data_wide)
   model_circularity <- lm(circularity_R ~ circularity_L, data = data_wide)
   
   #forloop for regression and nrmse, along with plots:
   # Define variables to analyze
   variables <- c("area", "circularity", "length", "width")

   # Initialize empty data frames
   regression_results <- data.frame()
   nrmse_results <- data.frame()
   
   # Loop through land covers
   for (lc in unique(data_wide$Land_cover)) {
     data_lc <- filter(data_wide, Land_cover == lc)
     
     for (var in variables) {
       L <- data_lc[[paste0(var, "_L")]]
       R <- data_lc[[paste0(var, "_R")]]
       
       # Remove NA and ensure numeric
       valid <- complete.cases(L, R)
       L <- as.numeric(L[valid])
       R <- as.numeric(R[valid])
       
       if (length(L) > 1) {
         # Linear regression
         model <- lm(R ~ L)
         r_squared <- summary(model)$r.squared
         
         # NRMSE (based on regression residuals)
         rmse <- sqrt(mean(residuals(model)^2))
         mean_val <- mean(c(L, R), na.rm = TRUE)
         nrmse <- rmse / mean_val
         
         # Save results
         regression_results <- rbind(regression_results, data.frame(
           Land_cover = lc,
           variable = var,
           r_squared = r_squared
         ))
         
         nrmse_results <- rbind(nrmse_results, data.frame(
           Land_cover = lc,
           variable = var,
           nrmse = nrmse
         ))
       }
     }
   }
   # -----------------------------
   # Plot: Regression R² values
   # -----------------------------
   ggplot(regression_results, aes(x = Land_cover, y = r_squared, fill = variable)) +
     geom_bar(stat = "identity", position = "dodge") +
     labs(
       title = "Regression R² by Variable and Land Cover",
       x = "Land Cover",
       y = "R² (Left vs Right Regression)",
       fill = "Variable"
     ) +
     theme_minimal()
   
   
   # -----------------------------
   # Plot: NRMSE values
   # -----------------------------
   ggplot(nrmse_results, aes(x = Land_cover, y = nrmse, fill = variable)) +
     geom_bar(stat = "identity", position = "dodge") +
     labs(
       title = "Normalized RMSE by Variable and Land Cover",
       x = "Land Cover",
       y = "NRMSE (Left vs Right)",
       fill = "Variable"
     ) +
     theme_minimal()


   # -----------------------------   
   ####By Cell 1 or 4
   # -----------------------------
   
   # Initialize empty data frames
   regression_results_1 <- data.frame()
   nrmse_results_1 <- data.frame()
   
   for (lc in unique(data_wide_1$Land_cover)) {
     data_lc <- filter(data_wide_1, Land_cover == lc)
     
     for (var in variables) {
       L <- data_lc[[paste0(var, "_L")]]
       R <- data_lc[[paste0(var, "_R")]]
       
       # Remove NA and ensure numeric
       valid <- complete.cases(L, R)
       L <- as.numeric(L[valid])
       R <- as.numeric(R[valid])
       
       if (length(L) > 1) {
         # Linear regression
         model <- lm(R ~ L)
         r_squared <- summary(model)$r.squared
         
         # NRMSE (based on regression residuals)
         rmse <- sqrt(mean(residuals(model)^2))
         mean_val <- mean(c(L, R), na.rm = TRUE)
         nrmse <- rmse / mean_val
         
         # Save results
         regression_results_1 <- rbind(regression_results_1, data.frame(
           Land_cover = lc,
           variable = var,
           r_squared = r_squared
         ))
         
         nrmse_results_1 <- rbind(nrmse_results_1, data.frame(
           Land_cover = lc,
           variable = var,
           nrmse = nrmse
         ))
       }
     }
   }
  
   
   # Initialize empty data frames
   regression_results_4 <- data.frame()
   nrmse_results_4 <- data.frame()
   
   for (lc in unique(data_wide_2$Land_cover)) {
     data_lc <- filter(data_wide_2, Land_cover == lc)
     
     for (var in variables) {
       L <- data_lc[[paste0(var, "_L")]]
       R <- data_lc[[paste0(var, "_R")]]
       
       # Remove NA and ensure numeric
       valid <- complete.cases(L, R)
       L <- as.numeric(L[valid])
       R <- as.numeric(R[valid])
       
       if (length(L) > 1) {
         # Linear regression
         model <- lm(R ~ L)
         r_squared <- summary(model)$r.squared
         
         # NRMSE (based on regression residuals)
         rmse <- sqrt(mean(residuals(model)^2))
         mean_val <- mean(c(L, R), na.rm = TRUE)
         nrmse <- rmse / mean_val
         
         # Save results
         regression_results_4 <- rbind(regression_results_4, data.frame(
           Land_cover = lc,
           variable = var,
           r_squared = r_squared
         ))
         
         nrmse_results_4 <- rbind(nrmse_results_4, data.frame(
           Land_cover = lc,
           variable = var,
           nrmse = nrmse
         ))
       }
     }
   }
   
   # Add Cell identifier
   regression_results_1 <- regression_results_1 %>% mutate(Cell = "Cell 1")
   regression_results_4 <- regression_results_4 %>% mutate(Cell = "Cell 4")
   
   # Combine
   regression_results_all <- bind_rows(regression_results_1, regression_results_4)
   regression_results_all <- regression_results_all %>%
     mutate(Cell = (gsub("cell ", "", Cell)))
   
   # Plot side-by-side bars grouped by Cell
   ggplot(regression_results_all, aes(x = Land_cover, y = r_squared, fill = variable)) +
     geom_bar(stat = "identity", position = "dodge") +
     facet_wrap(~Cell) +  # separate panels per variable if you have multiple
     labs(
       title = "Regression R² by Land Cover and Wing Cell",
       x = "Land Cover",
       y = "R² (Left vs Right Regression)",
       fill = "Wing Cell"
     ) +
     theme_minimal()

   # Add Cell identifier
   nrmse_results_1 <- nrmse_results_1 %>% mutate(Cell = "Cell 1")
   nrmse_results_4 <- nrmse_results_4 %>% mutate(Cell = "Cell 4")
   
   # Combine
   nrmse_results_all <- bind_rows(nrmse_results_1, nrmse_results_4)
   
   # Plot
   ggplot(nrmse_results_all, aes(x = Land_cover, y = nrmse, fill = variable)) +
     geom_bar(stat = "identity", position = "dodge") +
     facet_wrap(~Cell) +
     labs(
       title = "Normalized RMSE by Land Cover and Wing Cell",
       x = "Land Cover",
       y = "NRMSE (Left vs Right)",
       fill = "Wing Cell"
     ) +
     theme_minimal()   
   
   
   #PCA by NRMSE value
   
   pca_data <- nrmse_results_all %>%
     select(nrmse, variable, Land_cover) %>%
     pivot_wider(names_from = variable, values_from = nrmse) %>%
     na.omit()
   pca_result <- prcomp(pca_data %>% select(-Land_cover), center = TRUE, scale. = TRUE)
   
   summary(pca_result)         # To see variance explained
   pca_result$rotation         # Loadings (how variables contribute)
   pca_result$x                # The transformed PCA scores
   
   # Create a data frame with PCA scores and Land Cover
   pca_scores <- as.data.frame(pca_result$x)
   pca_scores$Land_cover <- nrmse_results_all$Land_cover[complete.cases(pca_data)]
   
   # Plot PC1 vs PC2 by land cover
   ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(Land_cover))) +
     geom_point(size = 3) +
     labs(title = "PCA of Wing Metrics by Land Cover",
          x = "Principal Component 1",
          y = "Principal Component 2",
          color = "Land_Cover") +
     theme_minimal()
   
   
   