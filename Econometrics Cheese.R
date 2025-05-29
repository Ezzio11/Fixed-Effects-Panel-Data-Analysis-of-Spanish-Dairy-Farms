# 1. INSTALL AND LOAD REQUIRED PACKAGES --------------------------------------

# Install packages if not already installed
install.packages(c("plm", "ggplot2", "dplyr", "GGally", "lmtest", "pbgtest", "patchwork"))

# Load required libraries
library(plm)       # Panel data econometrics
library(ggplot2)   # Data visualization
library(dplyr)     # Data manipulation
library(GGally)    # Extended ggplot functionality
library(lmtest)    # Hypothesis testing for linear models
library(car)       # Companion to Applied Regression
library(sandwich)  # Robust covariance matrix estimators
library(patchwork) # Combine multiple ggplot plots
library(magrittr)
library(tidyr)

# Function to detect outliers using IQR
detect_outliers_iqr <- function(data, var) {
  Q1 <- quantile(data[[var]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[var]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower <- Q1 - 1.5 * IQR
  upper <- Q3 + 1.5 * IQR
  data$outlier <- ifelse(data[[var]] < lower | data[[var]] > upper, TRUE, FALSE)
  data
}

# 2. LOAD AND INSPECT DATA ---------------------------------------------------

# Load dataset from CSV file
data <- read.csv("C:/Users/MSI/Downloads/dairy.csv")

# Quick check of data structure and summary statistics
head(data)     # View first few rows
summary(data)  # Summary statistics for all variables

# Check for missing data in all columns
data %>% 
  summarise(across(everything(), ~sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "Variable", values_to = "NA_Count")


# 3. DATA PREPARATION --------------------------------------------------------

# Convert YEAR to factor for proper categorical treatment
data$YEAR <- as.factor(data$YEAR)

# Define key variables for analysis
key_vars <- c("MILK", "FEED", "COWS", "LAND", "LABOR")
log_key_vars <- c("YIT","X1", "X2", "X3", "X4")


# 4. EXPLORATORY DATA ANALYSIS -----------------------------------------------

## 4.1 Pairwaise Correlation Matrix
ggpairs(data[log_key_vars])

## 4.1 Histograms of Key Variables ----
# Create histograms for all key variables
plot_list <- lapply(key_vars, function(var) {
  ggplot(data, aes_string(x = var)) +
    geom_histogram(bins = 30, fill = "grey50", color = "white") +
    labs(title = var, x = var, y = "Count") +
    theme_minimal()
})

plot_list_log <- lapply(log_key_vars, function(var) {
  ggplot(data, aes_string(x = var)) +
    geom_histogram(bins = 30, fill = "grey50", color = "white") +
    labs(title = var, x = var, y = "Count") +
    theme_minimal()
})

# Combine histograms using patchwork (3 columns)
wrap_plots(plotlist = plot_list, ncol = 3)
wrap_plots(plotlist = plot_list_log, ncol = 3) # Log variables

## 4.2 Boxplots of Key Variables ----
# Boxplots of original variables
boxplot_list <- lapply(key_vars, function(var) {
  ggplot(data, aes_string(y = var)) +
    geom_boxplot(fill = "gray70", color = "grey30", outlier.colour = "darkred") +
    labs(title = paste("Boxplot of", var), y = var) +
    theme_minimal()
})

# Boxplots of log-transformed variables
boxplot_list_log <- lapply(key_vars, function(var) {
  ggplot(data, aes_string(y = paste0("log(", var, ")"))) +
    geom_boxplot(fill = "gray70", color = "grey30", outlier.colour = "darkred") +
    labs(title = paste("Boxplot of log(", var, ")", sep = ""), y = paste("log(", var, ")", sep = "")) +
    theme_minimal()
})

# Combine and show boxplots in 3 columns
wrap_plots(plotlist = boxplot_list, ncol = 3)
wrap_plots(plotlist = boxplot_list_log, ncol = 3)

## 4.3 Yearly Boxplots ----
# Create boxplots showing distribution by year
boxplot_list <- lapply(key_vars, function(var) {
  ggplot(data, aes(x = factor(YEAR), y = var)) +
    geom_boxplot(fill = "gray70", color = "grey30", outlier.colour = "darkred") +
    labs(title = paste("Boxplot of", var), y = var) +
    theme_minimal()
})

# Combine yearly boxplots (3 columns)
wrap_plots(plotlist = boxplot_list, ncol = 3)

## 4.4 Time Series Plots ----
# Create line plots showing average values over time
lineplot_list <- lapply(key_vars, function(var) {
  data %>%
    group_by(YEAR) %>%
    summarise(avg_value = mean(.data[[var]], na.rm = TRUE)) %>%
    ggplot(aes(x = YEAR, y = avg_value, group = 1)) +
    geom_line(color = "grey70", size = 1) +
    geom_point(color = "steelblue") +
    labs(title = paste("Average", var, "Over Time"),
         x = "Year", y = paste("Average", var)) +
    theme_minimal()
})

# Combine time series plots (3 columns)
wrap_plots(plotlist = lineplot_list, ncol = 3)


# 5. PANEL DATA ANALYSIS -----------------------------------------------------

## 5.1 Prepare Panel Data Structure ----
# Convert to pdata.frame with FARM and YEAR as panel indices
pdata <- pdata.frame(data, index = c("FARM", "YEAR"))

## 5.2 Time Trend Visualization ----
# Plot average milk production over time
data %>%
  group_by(YEAR) %>%
  summarise(mean_YIT = mean(YIT, na.rm=TRUE)) %>%
  ggplot(aes(x=YEAR, y=mean_YIT)) + 
  geom_line(group=1) + geom_point() +
  labs(title="Average Log Milk Production Over Years")

## 5.3 Correlation Analysis ----
# Create correlation plot of input variables
ggpairs(data[,c("X1","X2","X3","X4")])

for (var in log_key_vars) {
  outliers <- detect_outliers_iqr(data, var)
  n_outliers <- sum(outliers$outlier)
  total <- nrow(data)
  pct_outliers <- round(100 * n_outliers / total, 2)
  
  cat(sprintf("Outliers in %s: %d data points (%.2f%% of total)\n", var, n_outliers, pct_outliers))
}

# 6. MODEL ESTIMATION and Dignosis --------------------------------------------------------

# Define translog production function formula

run_panel_models <- function(formula, pdata) {
  #--------------------------------------------------------------------------
  # Function: run_panel_models
  # Purpose:  Run and compare panel data models (Pooled OLS, FE, RE) with diagnostics
  # Inputs:   formula - model formula
  #           pdata   - panel data frame
  # Output:   Returns fixed effects model object
  #--------------------------------------------------------------------------
  
  # ====================== MODEL ESTIMATION =================================
  
  # 1. POOLED OLS MODEL
  cat("\n--- Running Pooled OLS ---\n")
  model_pooled <- plm(formula, data = pdata, model = "pooling")
  print(summary(model_pooled))
  
  # 2. RANDOM EFFECTS TEST (Breusch-Pagan LM test)
  cat("\n--- Lagrange Multiplier Test for Random Effects ---\n")
  print(plmtest(model_pooled, type = "bp"))
  
  # 3. FIXED EFFECTS MODEL
  cat("\n--- Running Fixed Effects Model ---\n")
  model_fe <- plm(formula, data = pdata, model = "within")
  print(summary(model_fe))
  
  # 4. RANDOM EFFECTS MODEL
  cat("\n--- Running Random Effects Model ---\n")
  model_re <- plm(formula, data = pdata, model = "random")
  print(summary(model_re))
  
  # 5. HAUSMAN TEST (FE vs RE comparison)
  cat("\n--- Hausman Test (FE vs RE) ---\n")
  print(phtest(model_fe, model_re))
  
  
  # ===================== MODEL DIAGNOSTICS ================================
  cat("\n--- Model Diagnostics (Fixed Effects Model) ---\n")
  
  # 1. RESIDUAL ANALYSIS
  fitted_values <- fitted(model_fe)
  residuals <- resid(model_fe)
  
  # 1.1 Normality check (QQ plot)
  qqnorm(residuals)
  qqline(residuals, col = "red", lwd = 2)
  
  # 1.2 Residuals vs fitted plot
  print(
    ggplot(data.frame(fitted = fitted_values, residuals = residuals), 
         aes(x = fitted, y = residuals)) +
      geom_point(alpha = 0.4) +
      geom_hline(yintercept = 0, color = "red") +
      labs(title = "Residuals vs Fitted (FE Model)", 
           x = "Fitted Values", 
           y = "Residuals") +
      theme_minimal()
  )
  
  # 2. HETEROSKEDASTICITY TEST
  cat("\nBreusch-Pagan Test for heteroskedasticity:\n")
  print(bptest(model_fe))
  
  # 3. SERIAL CORRELATION TEST
  cat("\nBreusch-Godfrey Test for serial correlation:\n")
  print(pbgtest(model_fe))
  
  # 4. MULTICOLLINEARITY CHECK
  cat("\nVariance Inflation Factors (VIF):\n")
  print(vif(lm(formula, data = pdata)))
  
  # 5. ROBUST STANDARD ERRORS
  cat("\nCoefficients with Cluster-Robust Standard Errors:\n")
  print(coeftest(model_fe, vcov = vcovHC(model_fe, type = "HC3", cluster = c("group", "time"))))
  
  
  # ===================== RETURN RESULTS ===================================
  
  # Return fixed effects model for further analysis
  return(model_fe)
}
# RQ1. What is the impact of core inputs (COWS, LAND, LABOR, FEED) on milk production over time?
formula_translog <- YIT ~ X1 + X2 + X3 + X4
run_panel_models(formula_translog, pdata)

# RQ2. Are there diminishing or increasing returns to scale in milk production?
formula_translog2 <- YIT ~ X1 + X2 + X3 + X4 + X11 + X22 + X33 + X44
run_panel_models(formula_translog2, pdata)

# RQ3. How do interactions between inputs affect milk production efficiency?
formula_translog3 <- YIT ~ X1 + X2 + X3 + X4 + X12 + X13 + X14 + X23 + X24 + X34
run_panel_models(formula_translog3, pdata)

# RQ4. Are there time-fixed effects that impact milk production?
formula_translog4 <- YIT ~ X1 + X2 + X3 + X4 + factor(YEAR)
run_panel_models(formula_translog4, pdata)