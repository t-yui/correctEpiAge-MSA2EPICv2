library(tidyverse)
library(ggplot2)
library(knitr)
library(gridExtra)
library(cowplot)

# settings

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

clocks <- c("DNAmHorvathClock", "DNAmHannumClock", "DNAmPhenoAge", 
            "DNAmGrimAge", "DNAmGrimAgev2", "DNAmDunedinPACE")

clock_labels <- c("Horvath", "Hannum", "PhenoAge", 
                  "GrimAge", "GrimAge v2", "DunedinPACE")

# functions for metrics

rmse <- function(truth, estimate) sqrt(mean((estimate - truth)^2, na.rm = TRUE))
mae <- function(truth, estimate) mean(abs(estimate - truth), na.rm = TRUE)
bias <- function(truth, estimate) mean(estimate - truth, na.rm = TRUE)

# functions for plot

bland_altman_plot <- function(truth, estimate, title = "", xlabel = NULL, ylabel = NULL) {
  df <- data.frame(
    mean_val = (truth + estimate) / 2,
    diff = estimate - truth
  )
  
  m <- mean(df$diff, na.rm = TRUE)
  sd <- sd(df$diff, na.rm = TRUE)
  ul <- m + 1.96 * sd
  ll <- m - 1.96 * sd
  
  if(is.null(xlabel)) xlabel <- "Mean of EPICv2 and MSA"
  if(is.null(ylabel)) ylabel <- "Difference (MSA - EPICv2)"
  
  p <- ggplot(df, aes(x = mean_val, y = diff)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray50", size = 0.5) +
    geom_hline(yintercept = m, linetype = 1, color = "blue", size = 0.5) +
    geom_hline(yintercept = ul, linetype = 3, color = "red", size = 0.5) +
    geom_hline(yintercept = ll, linetype = 3, color = "red", size = 0.5) +
    labs(title = title, x = xlabel, y = ylabel) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

scatter_plot <- function(truth, estimate, title = "", xlabel = NULL, ylabel = NULL) {
  r_val <- cor(truth, estimate, use = "complete.obs")
  
  if(is.null(xlabel)) xlabel <- "EPICv2"
  if(is.null(ylabel)) ylabel <- "MSA"
  
  p <- ggplot(data.frame(truth = truth, estimate = estimate), 
              aes(x = truth, y = estimate)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50", size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", size = 0.5) +
    labs(
      title = title,
      subtitle = sprintf("r = %.3f", r_val),
      x = xlabel,
      y = ylabel
    ) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}


# perform analysis

all_train_data <- list()
all_test_data <- list()
all_models <- list()
all_plots <- list()

for(i in 1:length(clocks)) {
  clock <- clocks[i]
  clock_label <- clock_labels[i]
  
  train_file <- sprintf("inputs/251104_%s_MSA_vs_EPICv2_training.csv", clock)
  train_data <- read.csv(train_file, stringsAsFactors = FALSE)
  
  test_file <- sprintf("inputs/251104_%s_MSA_vs_EPICv2_testing.csv", clock)
  test_data <- read.csv(test_file, stringsAsFactors = FALSE)
  
  train_data <- train_data %>%
    transmute(
      Name = as.character(Name),
      MSA = as.numeric(MSA),
      EPICv2 = as.numeric(EPICv2),
      Age = as.numeric(Chronological_Age),
      Sex = ifelse(Gender == "Male", 1, 0)
    ) %>%
    filter(!is.na(MSA) & !is.na(EPICv2))
  
  test_data <- test_data %>%
    transmute(
      Name = as.character(Name),
      MSA = as.numeric(MSA),
      EPICv2 = as.numeric(EPICv2),
      Age = as.numeric(Chronological_Age),
      Sex = ifelse(Gender == "Male", 1, 0)
    ) %>%
    filter(!is.na(MSA) & !is.na(EPICv2))
  
  all_train_data[[clock]] <- train_data
  all_test_data[[clock]] <- test_data
  
  ## fit models to taining data
  model1 <- lm(EPICv2 ~ 1 + offset(MSA), data = train_data)
  model2 <- lm(EPICv2 ~ MSA, data = train_data)
  model3 <- lm(EPICv2 ~ MSA + Age + Sex, data = train_data)
  
  all_models[[clock]] <- list(model1 = model1, model2 = model2, model3 = model3)
  
  ## apply corrections to test data
  test_data$corrected_model1 <- predict(model1, newdata = test_data)
  test_data$corrected_model2 <- predict(model2, newdata = test_data)
  test_data$corrected_model3 <- predict(model3, newdata = test_data)
  
  all_test_data[[clock]] <- test_data
}

# figure 1: EPICv2 vs MSA scatter and Bland-Altman plots

plot_list_fig1 <- list()
for(i in 1:length(clocks)) {
  clock_label <- clock_labels[i]
  train_data <- all_train_data[[clocks[i]]]
  
  # Scatter plot
  p_scatter <- scatter_plot(train_data$EPICv2, train_data$MSA, 
                            clock_label,
                            xlabel = "EPICv2",
                            ylabel = "MSA")
  
  # Bland-Altman plot
  p_ba <- bland_altman_plot(train_data$EPICv2, train_data$MSA,
                            "",
                            xlabel = "Mean of EPICv2 and MSA",
                            ylabel = "Difference (MSA - EPICv2)")
  
  plot_list_fig1[[2*(i-1)+1]] <- p_scatter
  plot_list_fig1[[2*(i-1)+2]] <- p_ba
}

pdf("results/figures/Figure1_EPICv2_vs_MSA_comparison.pdf", width = 8, height = 12)
grid.arrange(grobs = plot_list_fig1, ncol = 2, nrow = 6,
             top = "EPICv2 vs MSA: Scatter plots and Bland-Altman plots")
dev.off()

# figure 2: EPICv2 vs Corrected MSA scatter plots

plot_list_fig2 <- list()
for(i in 1:length(clocks)) {
  clock_label <- clock_labels[i]
  test_data <- all_test_data[[clocks[i]]]
  p1 <- scatter_plot(test_data$EPICv2, test_data$corrected_model1,
                     paste(clock_label, "- Model 1"),
                     xlabel = "EPICv2",
                     ylabel = "Corrected MSA (Model 1)")
  p2 <- scatter_plot(test_data$EPICv2, test_data$corrected_model2,
                     paste(clock_label, "- Model 2"),
                     xlabel = "EPICv2",
                     ylabel = "Corrected MSA (Model 2)")
  p3 <- scatter_plot(test_data$EPICv2, test_data$corrected_model3,
                     paste(clock_label, "- Model 3"),
                     xlabel = "EPICv2",
                     ylabel = "Corrected MSA (Model 3)")
  
  plot_list_fig2[[3*(i-1)+1]] <- p1
  plot_list_fig2[[3*(i-1)+2]] <- p2
  plot_list_fig2[[3*(i-1)+3]] <- p3
}

pdf("results/figures/Figure2_EPICv2_vs_CorrectedMSA_scatter.pdf", width = 9, height = 12)
grid.arrange(grobs = plot_list_fig2, ncol = 3, nrow = 6,
             top = "EPICv2 vs Corrected MSA: Scatter plots")
dev.off()

# figure 3: EPICv2 vs Corrected MSA Bland-Altman plots

plot_list_fig3 <- list()
for(i in 1:length(clocks)) {
  clock_label <- clock_labels[i]
  test_data <- all_test_data[[clocks[i]]]
  p1 <- bland_altman_plot(test_data$EPICv2, test_data$corrected_model1,
                          paste(clock_label, "- Model 1"),
                          xlabel = "Mean of EPICv2 and Corrected MSA",
                          ylabel = "Difference (Corrected - EPICv2)")
  p2 <- bland_altman_plot(test_data$EPICv2, test_data$corrected_model2,
                          paste(clock_label, "- Model 2"),
                          xlabel = "Mean of EPICv2 and Corrected MSA",
                          ylabel = "Difference (Corrected - EPICv2)")
  p3 <- bland_altman_plot(test_data$EPICv2, test_data$corrected_model3,
                          paste(clock_label, "- Model 3"),
                          xlabel = "Mean of EPICv2 and Corrected MSA",
                          ylabel = "Difference (Corrected - EPICv2)")
  
  plot_list_fig3[[3*(i-1)+1]] <- p1
  plot_list_fig3[[3*(i-1)+2]] <- p2
  plot_list_fig3[[3*(i-1)+3]] <- p3
}

pdf("results/figures/Figure3_EPICv2_vs_CorrectedMSA_BlandAltman.pdf", width = 9, height = 12)
grid.arrange(grobs = plot_list_fig3, ncol = 3, nrow = 6,
             top = "EPICv2 vs Corrected MSA: Bland-Altman plots")
dev.off()


# table 1: training data characteristics

first_train <- all_train_data[[1]]
n_train <- nrow(first_train)
n_male <- sum(first_train$Sex == 1)
n_female <- sum(first_train$Sex == 0)
tbl1 <- data.frame(
  Characteristic = c(
    "Sample size, n",
    "Male, n (%)",
    "Female, n (%)",
    "Chronological age, mean (SD)",
    "Chronological age, range"
  ),
  Value = c(
    sprintf("%d", n_train),
    sprintf("%d (%.1f%%)", n_male, 100 * n_male / n_train),
    sprintf("%d (%.1f%%)", n_female, 100 * n_female / n_train),
    sprintf("%.1f (%.1f)", mean(first_train$Age), sd(first_train$Age)),
    sprintf("%.1f - %.1f", min(first_train$Age), max(first_train$Age))
  ),
  stringsAsFactors = FALSE
)
for (i in 1:length(clocks)) {
  clock_label <- clock_labels[i]
  train_data <- all_train_data[[clocks[i]]]
  
  tbl1 <- rbind(
    tbl1,
    data.frame(
      Characteristic = sprintf("%s (EPICv2), mean (SD)", clock_label),
      Value = sprintf("%.2f (%.2f)", mean(train_data$EPICv2), sd(train_data$EPICv2)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Characteristic = sprintf("%s (MSA), mean (SD)", clock_label),
      Value = sprintf("%.2f (%.2f)", mean(train_data$MSA), sd(train_data$MSA)),
      stringsAsFactors = FALSE
    )
  )
}

print(tbl1, row.names = FALSE)


# table 2: Test data characteristics

first_test <- all_test_data[[1]]
n_test <- nrow(first_test)
n_male_test <- sum(first_test$Sex == 1)
n_female_test <- sum(first_test$Sex == 0)
tbl2 <- data.frame(
  Characteristic = c(
    "Sample size, n",
    "Male, n (%)",
    "Female, n (%)",
    "Chronological age, mean (SD)",
    "Chronological age, range"
  ),
  Value = c(
    sprintf("%d", n_test),
    sprintf("%d (%.1f%%)", n_male_test, 100 * n_male_test / n_test),
    sprintf("%d (%.1f%%)", n_female_test, 100 * n_female_test / n_test),
    sprintf("%.1f (%.1f)", mean(first_test$Age), sd(first_test$Age)),
    sprintf("%.1f - %.1f", min(first_test$Age), max(first_test$Age))
  ),
  stringsAsFactors = FALSE
)
for (i in 1:length(clocks)) {
  clock_label <- clock_labels[i]
  test_data <- all_test_data[[clocks[i]]]
  
  tbl2 <- rbind(
    tbl2,
    data.frame(
      Characteristic = sprintf("%s (EPICv2), mean (SD)", clock_label),
      Value = sprintf("%.2f (%.2f)", mean(test_data$EPICv2), sd(test_data$EPICv2)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Characteristic = sprintf("%s (MSA), mean (SD)", clock_label),
      Value = sprintf("%.2f (%.2f)", mean(test_data$MSA), sd(test_data$MSA)),
      stringsAsFactors = FALSE
    )
  )
}

print(tbl2, row.names = FALSE)


# table 3: MSA - EPICv2 systematic differences

mean_diff3 <- numeric(length(clocks))
sd_diff3   <- numeric(length(clocks))

for (i in 1:length(clocks)) {
  train_data <- all_train_data[[clocks[i]]]
  diff <- train_data$MSA - train_data$EPICv2
  mean_diff3[i] <- mean(diff)
  sd_diff3[i]   <- sd(diff)
}

tbl3 <- data.frame(
  `Epigenetic Clock` = clock_labels,
  `Mean Difference`  = sprintf("%.3f", mean_diff3),
  SD                 = sprintf("%.3f", sd_diff3),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

print(tbl3, row.names = FALSE)


# table 4: regression coefficients for all models

rows4 <- list()

for (i in 1:length(clocks)) {
  clock_label <- clock_labels[i]
  models <- all_models[[clocks[i]]]
  
  ## Model 1
  s1 <- summary(models$model1)
  coef1 <- s1$coefficients
  rows4[[length(rows4) + 1]] <- data.frame(
    Clock    = clock_label,
    Model    = "Model 1",
    Term     = "Intercept",
    Estimate = sprintf("%.4f", coef1[1, 1]),
    Std.Error= sprintf("%.4f", coef1[1, 2]),
    t.value  = sprintf("%.3f", coef1[1, 3]),
    p.value  = sprintf("%.3e", coef1[1, 4]),
    stringsAsFactors = FALSE
  )
  
  ## Model 2
  s2 <- summary(models$model2)
  coef2 <- s2$coefficients
  rows4[[length(rows4) + 1]] <- data.frame(
    Clock    = clock_label,
    Model    = "Model 2",
    Term     = "Intercept",
    Estimate = sprintf("%.4f", coef2[1, 1]),
    Std.Error= sprintf("%.4f", coef2[1, 2]),
    t.value  = sprintf("%.3f", coef2[1, 3]),
    p.value  = sprintf("%.3e", coef2[1, 4]),
    stringsAsFactors = FALSE
  )
  rows4[[length(rows4) + 1]] <- data.frame(
    Clock    = clock_label,
    Model    = "Model 2",
    Term     = "MSA",
    Estimate = sprintf("%.4f", coef2[2, 1]),
    Std.Error= sprintf("%.4f", coef2[2, 2]),
    t.value  = sprintf("%.3f", coef2[2, 3]),
    p.value  = sprintf("%.3e", coef2[2, 4]),
    stringsAsFactors = FALSE
  )
  
  ## Model 3
  s3 <- summary(models$model3)
  coef3 <- s3$coefficients
  rows4[[length(rows4) + 1]] <- data.frame(
    Clock    = clock_label,
    Model    = "Model 3",
    Term     = "Intercept",
    Estimate = sprintf("%.4f", coef3[1, 1]),
    Std.Error= sprintf("%.4f", coef3[1, 2]),
    t.value  = sprintf("%.3f", coef3[1, 3]),
    p.value  = sprintf("%.3e", coef3[1, 4]),
    stringsAsFactors = FALSE
  )
  rows4[[length(rows4) + 1]] <- data.frame(
    Clock    = clock_label,
    Model    = "Model 3",
    Term     = "MSA",
    Estimate = sprintf("%.4f", coef3[2, 1]),
    Std.Error= sprintf("%.4f", coef3[2, 2]),
    t.value  = sprintf("%.3f", coef3[2, 3]),
    p.value  = sprintf("%.3e", coef3[2, 4]),
    stringsAsFactors = FALSE
  )
  rows4[[length(rows4) + 1]] <- data.frame(
    Clock    = clock_label,
    Model    = "Model 3",
    Term     = "Age",
    Estimate = sprintf("%.4f", coef3[3, 1]),
    Std.Error= sprintf("%.4f", coef3[3, 2]),
    t.value  = sprintf("%.3f", coef3[3, 3]),
    p.value  = sprintf("%.3e", coef3[3, 4]),
    stringsAsFactors = FALSE
  )
  rows4[[length(rows4) + 1]] <- data.frame(
    Clock    = clock_label,
    Model    = "Model 3",
    Term     = "Sex",
    Estimate = sprintf("%.4f", coef3[4, 1]),
    Std.Error= sprintf("%.4f", coef3[4, 2]),
    t.value  = sprintf("%.3f", coef3[4, 3]),
    p.value  = sprintf("%.3e", coef3[4, 4]),
    stringsAsFactors = FALSE
  )
}

tbl4 <- do.call(rbind, rows4)
print(tbl4, row.names = FALSE)


# table 5: comprehensive model performance evaluation

rows5 <- list()

for (i in 1:length(clocks)) {
  clock_label <- clock_labels[i]
  test_data <- all_test_data[[clocks[i]]]
  
  ## Uncorrected
  mae_uncorr <- mae(test_data$EPICv2, test_data$MSA)
  rmse_uncorr <- rmse(test_data$EPICv2, test_data$MSA)
  bias_uncorr <- bias(test_data$EPICv2, test_data$MSA)
  cor_uncorr <- cor(test_data$EPICv2, test_data$MSA, use = "complete.obs")
  
  rows5[[length(rows5) + 1]] <- data.frame(
    Clock      = clock_label,
    Model      = "Uncorrected",
    MAE        = sprintf("%.3f", mae_uncorr),
    RMSE       = sprintf("%.3f", rmse_uncorr),
    Bias       = sprintf("%.3f", bias_uncorr),
    Correlation= sprintf("%.3f", cor_uncorr),
    stringsAsFactors = FALSE
  )
  
  ## Model 1
  mae1 <- mae(test_data$EPICv2, test_data$corrected_model1)
  rmse1 <- rmse(test_data$EPICv2, test_data$corrected_model1)
  bias1 <- bias(test_data$EPICv2, test_data$corrected_model1)
  cor1 <- cor(test_data$EPICv2, test_data$corrected_model1, use = "complete.obs")
  
  rows5[[length(rows5) + 1]] <- data.frame(
    Clock      = clock_label,
    Model      = "Model 1",
    MAE        = sprintf("%.3f", mae1),
    RMSE       = sprintf("%.3f", rmse1),
    Bias       = sprintf("%.3f", bias1),
    Correlation= sprintf("%.3f", cor1),
    stringsAsFactors = FALSE
  )
  
  ## Model 2
  mae2 <- mae(test_data$EPICv2, test_data$corrected_model2)
  rmse2 <- rmse(test_data$EPICv2, test_data$corrected_model2)
  bias2 <- bias(test_data$EPICv2, test_data$corrected_model2)
  cor2 <- cor(test_data$EPICv2, test_data$corrected_model2, use = "complete.obs")
  
  rows5[[length(rows5) + 1]] <- data.frame(
    Clock      = clock_label,
    Model      = "Model 2",
    MAE        = sprintf("%.3f", mae2),
    RMSE       = sprintf("%.3f", rmse2),
    Bias       = sprintf("%.3f", bias2),
    Correlation= sprintf("%.3f", cor2),
    stringsAsFactors = FALSE
  )
  
  ## Model 3
  mae3 <- mae(test_data$EPICv2, test_data$corrected_model3)
  rmse3 <- rmse(test_data$EPICv2, test_data$corrected_model3)
  bias3 <- bias(test_data$EPICv2, test_data$corrected_model3)
  cor3 <- cor(test_data$EPICv2, test_data$corrected_model3, use = "complete.obs")
  
  rows5[[length(rows5) + 1]] <- data.frame(
    Clock      = clock_label,
    Model      = "Model 3",
    MAE        = sprintf("%.3f", mae3),
    RMSE       = sprintf("%.3f", rmse3),
    Bias       = sprintf("%.3f", bias3),
    Correlation= sprintf("%.3f", cor3),
    stringsAsFactors = FALSE
  )
}

tbl5 <- do.call(rbind, rows5)
print(tbl5, row.names = FALSE)

