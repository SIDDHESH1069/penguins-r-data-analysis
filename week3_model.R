# ==========================================================================
# Week 3 - Statistical Analysis and Predictive Modeling using R
# Dataset: Palmer Archipelago (Antarctica) Penguins (CC0)
# Goal   : (A) formal hypothesis testing, (B) a cross-validated classifier
#          that predicts a penguin's SEX from its body measurements,
#          (C) full diagnostics and optimisation.
# Author : Siddhesh Mhatre
# ==========================================================================
library(ggplot2)
library(MASS)

set.seed(42)
OUT <- "/tmp/w3/figs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sv <- function(p, file, w = 9, h = 5.6, dpi = 200)
  ggsave(file.path(OUT, file), p, width = w, height = h, dpi = dpi, bg = "white")

pal  <- c(Adelie = "#E1701A", Chinstrap = "#7B3FA0", Gentoo = "#1B7F79")
psex <- c(female = "#C2185B", male = "#1565C0")
theme_pg <- theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(colour = "grey35", size = 11),
        plot.caption = element_text(colour = "grey45", size = 9),
        panel.grid.minor = element_blank(), legend.position = "top")
theme_set(theme_pg)

# ---- 1. Load and prepare (Week 1 cleaning pipeline) ----------------------
penguins <- read.csv("/tmp/w3/penguins.csv", stringsAsFactors = FALSE)
penguins[penguins == ""] <- NA
df <- penguins[!is.na(penguins$bill_length_mm), ]
df$species <- factor(df$species, levels = c("Adelie", "Chinstrap", "Gentoo"))
df$island  <- factor(df$island)
df$year    <- factor(df$year)
# Modelling frame: sex is the response, so the 9 sex-unknown rows are held out
model_df <- df[!is.na(df$sex), ]
model_df$sex <- factor(model_df$sex, levels = c("female", "male"))

cat("== Dimensions: cleaned / modelling ==\n"); print(dim(df)); print(dim(model_df))
cat("\n== str(model_df) ==\n"); str(model_df)
cat("\n== summary of predictors ==\n")
print(summary(model_df[, c("bill_length_mm", "bill_depth_mm",
                           "flipper_length_mm", "body_mass_g")]))
cat("\n== Class balance (sex) ==\n"); print(table(model_df$sex))
cat("\n== sex x species ==\n"); print(table(model_df$sex, model_df$species))

num_vars <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")

# ==========================================================================
# 2. EXPLORATORY STATISTICAL ANALYSIS / HYPOTHESIS TESTING
# ==========================================================================
cat("\n\n########## 2. HYPOTHESIS TESTING ##########\n")

cat("\n== H0-1 Shapiro-Wilk normality tests (whole sample) ==\n")
sw <- t(sapply(num_vars, function(v) {
  s <- shapiro.test(model_df[[v]]); c(W = unname(s$statistic), p = s$p.value)
}))
print(round(sw, 4))

cat("\n== Shapiro-Wilk within species x sex strata (body mass) ==\n")
strata <- interaction(model_df$species, model_df$sex, drop = TRUE)
print(round(sapply(split(model_df$body_mass_g, strata),
                   function(x) shapiro.test(x)$p.value), 4))

cat("\n== H0-2 Equality of variances (Bartlett, mass ~ sex) ==\n")
print(bartlett.test(body_mass_g ~ sex, data = model_df))

cat("\n== H0-3 Welch two-sample t-tests: predictor ~ sex ==\n")
tt <- t(sapply(num_vars, function(v) {
  f  <- as.formula(paste(v, "~ sex"))
  tr <- t.test(f, data = model_df)
  wr <- wilcox.test(f, data = model_df)
  d  <- diff(rev(tr$estimate))
  sp <- sd(model_df[[v]])
  c(mean_female = unname(tr$estimate[1]), mean_male = unname(tr$estimate[2]),
    diff = unname(-d), t = unname(tr$statistic), df = unname(tr$parameter),
    p_t = tr$p.value, p_wilcox = wr$p.value, cohens_d = unname(-d) / sp)
}))
print(round(tt, 4))

cat("\n== H0-4 Two-way ANOVA: body mass ~ species * sex ==\n")
print(summary(aov(body_mass_g ~ species * sex, data = model_df)))

cat("\n== H0-5 Chi-square: is sex independent of species? ==\n")
print(chisq.test(table(model_df$species, model_df$sex)))

cat("\n== H0-6 Pearson correlations among predictors ==\n")
cm <- cor(model_df[, num_vars]); print(round(cm, 3))
cat("\nPairwise correlation tests (r, p):\n")
cmb <- combn(num_vars, 2)
print(round(t(apply(cmb, 2, function(p) {
  ct <- cor.test(model_df[[p[1]]], model_df[[p[2]]])
  c(r = unname(ct$estimate), p = ct$p.value)
})), 4))
cat("\nVariance inflation proxy - correlation of predictors within species:\n")
print(round(sapply(split(model_df, model_df$species),
                   function(d) cor(d$flipper_length_mm, d$body_mass_g)), 3))

# ---- Figure 1: density of each predictor by sex --------------------------
long <- do.call(rbind, lapply(num_vars, function(v)
  data.frame(variable = v, value = model_df[[v]], sex = model_df$sex)))
long$variable <- factor(long$variable, levels = num_vars,
                        labels = c("Bill length (mm)", "Bill depth (mm)",
                                   "Flipper length (mm)", "Body mass (g)"))
p1 <- ggplot(long, aes(value, fill = sex, colour = sex)) +
  geom_density(alpha = .35, linewidth = .8) +
  facet_wrap(~ variable, scales = "free") +
  scale_fill_manual(values = psex) + scale_colour_manual(values = psex) +
  labs(title = "Every body measurement separates the sexes - but none perfectly",
       subtitle = "Kernel densities by sex; all four Welch t-tests are significant at p < 0.001",
       x = NULL, y = "Density", fill = "Sex", colour = "Sex",
       caption = "Data: palmerpenguins (CC0), n = 333")
sv(p1, "fig01_predictor_densities.png", h = 6.2)

# ---- Figure 2: QQ plots (normality assumption) ---------------------------
qq <- do.call(rbind, lapply(num_vars, function(v) {
  x <- model_df[[v]]; q <- qqnorm(x, plot.it = FALSE)
  data.frame(variable = v, theoretical = q$x, sample = scale(q$y)[, 1])
}))
qq$variable <- factor(qq$variable, levels = num_vars,
                      labels = c("Bill length", "Bill depth",
                                 "Flipper length", "Body mass"))
p2 <- ggplot(qq, aes(theoretical, sample)) +
  geom_abline(slope = 1, intercept = 0, colour = "#B2182B", linewidth = .8) +
  geom_point(size = 1.3, alpha = .6, colour = "grey25") +
  facet_wrap(~ variable) +
  labs(title = "Normal Q-Q plots of the four candidate predictors",
       subtitle = "Points track the reference line closely; mild right skew in body mass",
       x = "Theoretical quantiles", y = "Standardised sample quantiles")
sv(p2, "fig02_qq_normality.png", h = 6.2)

# ---- Figure 3: boxplots by species and sex -------------------------------
p3 <- ggplot(model_df, aes(species, body_mass_g, fill = sex)) +
  geom_boxplot(width = .65, alpha = .75, outlier.alpha = .4) +
  scale_fill_manual(values = psex) +
  labs(title = "The sex gap in body mass holds inside every species",
       subtitle = "Two-way ANOVA: species F = 758.4, sex F = 387.5, interaction p = 0.0002",
       x = NULL, y = "Body mass (g)", fill = "Sex")
sv(p3, "fig03_mass_species_sex_box.png")

# ---- Figure 4: correlation heat map --------------------------------------
cml <- data.frame(as.table(round(cm, 2))); names(cml) <- c("v1", "v2", "r")
p4 <- ggplot(cml, aes(v1, v2, fill = r)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  labs(title = "Predictor correlation structure",
       subtitle = "Flipper length and body mass are strongly collinear (r = 0.87)",
       x = NULL, y = NULL, fill = "r") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
sv(p4, "fig04_corr_heatmap.png", w = 8.6, h = 6.4)

# ==========================================================================
# 3. MODEL BUILDING - logistic regression, 10-fold stratified CV
# ==========================================================================
cat("\n\n########## 3. MODEL BUILDING ##########\n")

forms <- list(
  M1_mass_only  = sex ~ body_mass_g,
  M2_all_num    = sex ~ bill_length_mm + bill_depth_mm + flipper_length_mm + body_mass_g,
  M3_plus_spp   = sex ~ bill_length_mm + bill_depth_mm + flipper_length_mm +
                        body_mass_g + species,
  M4_spp_inter  = sex ~ species * (bill_depth_mm + body_mass_g) +
                        bill_length_mm + flipper_length_mm
)

# stratified fold assignment
make_folds <- function(y, k = 10) {
  fold <- integer(length(y))
  for (lv in levels(y)) {
    idx <- which(y == lv)
    fold[sample(idx)] <- rep_len(1:k, length(idx))
  }
  fold
}
folds <- make_folds(model_df$sex, 10)
cat("\n== Fold sizes ==\n"); print(table(folds, model_df$sex))

metrics <- function(truth, prob, thr = 0.5) {
  pred <- factor(ifelse(prob >= thr, "male", "female"), levels = levels(truth))
  cmx  <- table(Predicted = pred, Actual = truth)
  TP <- cmx["male", "male"]; TN <- cmx["female", "female"]
  FP <- cmx["male", "female"]; FN <- cmx["female", "male"]
  prec <- TP / (TP + FP); rec <- TP / (TP + FN)
  o <- order(prob); r <- rank(prob, ties.method = "average")
  n1 <- sum(truth == "male"); n0 <- sum(truth == "female")
  auc <- (sum(r[truth == "male"]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  c(accuracy = (TP + TN) / length(truth), sensitivity = rec,
    specificity = TN / (TN + FP), precision = prec,
    f1 = 2 * prec * rec / (prec + rec), auc = auc,
    brier = mean((prob - (truth == "male")) ^ 2))
}

cv_glm <- function(form, data, folds) {
  prob <- numeric(nrow(data))
  for (k in sort(unique(folds))) {
    fit <- glm(form, data = data[folds != k, ], family = binomial())
    prob[folds == k] <- predict(fit, newdata = data[folds == k, ], type = "response")
  }
  prob
}

cv_prob <- lapply(forms, cv_glm, data = model_df, folds = folds)
cv_tab  <- t(sapply(cv_prob, metrics, truth = model_df$sex))
cat("\n== 10-fold cross-validated performance (threshold 0.5) ==\n")
print(round(cv_tab, 4))

# LDA benchmark on the same folds
cv_lda <- numeric(nrow(model_df))
for (k in sort(unique(folds))) {
  fit <- lda(forms$M3_plus_spp, data = model_df[folds != k, ])
  cv_lda[folds == k] <- predict(fit, newdata = model_df[folds == k, ])$posterior[, "male"]
}
cat("\n== LDA benchmark (same folds, M3 predictors) ==\n")
print(round(metrics(model_df$sex, cv_lda), 4))

# per-fold accuracy of the champion model
best <- "M3_plus_spp"
fold_acc <- sapply(sort(unique(folds)), function(k)
  mean((cv_prob[[best]][folds == k] >= .5) == (model_df$sex[folds == k] == "male")))
cat("\n== Per-fold accuracy of", best, "==\n"); print(round(fold_acc, 4))
cat("mean =", round(mean(fold_acc), 4), " sd =", round(sd(fold_acc), 4), "\n")

# final model fitted on all data
final <- glm(forms[[best]], data = model_df, family = binomial())
cat("\n== summary(final model) ==\n"); print(summary(final))
cat("\n== Odds ratios with 95% profile CI ==\n")
or <- exp(cbind(OR = coef(final), suppressMessages(confint(final))))
print(round(or, 4))
cat("\n== Likelihood-ratio test vs null, and vs M2 (no species) ==\n")
m2 <- glm(forms$M2_all_num, data = model_df, family = binomial())
print(anova(m2, final, test = "Chisq"))
cat("\nMcFadden pseudo R2 =",
    round(1 - final$deviance / final$null.deviance, 4),
    " AIC:", round(AIC(m2), 1), "->", round(AIC(final), 1), "\n")

# ---- Figure 5: model comparison bars -------------------------------------
mc <- data.frame(model = rep(rownames(cv_tab), 2),
                 metric = rep(c("Accuracy", "AUC"), each = nrow(cv_tab)),
                 value = c(cv_tab[, "accuracy"], cv_tab[, "auc"]))
p5 <- ggplot(mc, aes(model, value, fill = metric)) +
  geom_col(position = position_dodge(.75), width = .7) +
  geom_text(aes(label = sprintf("%.3f", value)),
            position = position_dodge(.75), vjust = -.4, size = 3.4) +
  scale_fill_manual(values = c(Accuracy = "#1B7F79", AUC = "#E1701A")) +
  coord_cartesian(ylim = c(0.5, 1.03)) +
  labs(title = "Cross-validated accuracy and AUC for four candidate specifications",
       subtitle = "Adding species to the four measurements gives the largest single gain",
       x = NULL, y = "10-fold CV score", fill = NULL)
sv(p5, "fig05_model_comparison.png")

# ---- Figure 6: ROC curves ------------------------------------------------
roc_pts <- function(truth, prob, label) {
  thr <- sort(unique(c(0, prob, 1)), decreasing = TRUE)
  do.call(rbind, lapply(thr, function(t) {
    pr <- prob >= t
    data.frame(model = label,
               fpr = sum(pr & truth == "female") / sum(truth == "female"),
               tpr = sum(pr & truth == "male") / sum(truth == "male"))
  }))
}
roc_df <- rbind(roc_pts(model_df$sex, cv_prob$M1_mass_only,
                        sprintf("M1 mass only (AUC %.3f)", cv_tab["M1_mass_only", "auc"])),
                roc_pts(model_df$sex, cv_prob$M2_all_num,
                        sprintf("M2 measurements (AUC %.3f)", cv_tab["M2_all_num", "auc"])),
                roc_pts(model_df$sex, cv_prob$M3_plus_spp,
                        sprintf("M3 + species (AUC %.3f)", cv_tab["M3_plus_spp", "auc"])),
                roc_pts(model_df$sex, cv_lda,
                        sprintf("LDA benchmark (AUC %.3f)", metrics(model_df$sex, cv_lda)["auc"])))
p6 <- ggplot(roc_df, aes(fpr, tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 1) + coord_equal() +
  scale_colour_manual(values = c("#9E9E9E", "#E1701A", "#1B7F79", "#7B3FA0")) +
  labs(title = "ROC curves from out-of-fold predictions",
       subtitle = "Every curve is built only from predictions the model never trained on",
       x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)",
       colour = NULL) +
  theme(legend.position = "right")
sv(p6, "fig06_roc_curves.png", w = 9, h = 6)

# ---- Figure 7: confusion matrix of the champion --------------------------
pred_best <- factor(ifelse(cv_prob[[best]] >= .5, "male", "female"),
                    levels = levels(model_df$sex))
cmx <- table(Predicted = pred_best, Actual = model_df$sex)
cat("\n== Confusion matrix (out-of-fold,", best, ") ==\n"); print(cmx)
cmd <- data.frame(cmx); cmd$rate <- cmd$Freq / rep(colSums(cmx), each = 2)
p7 <- ggplot(cmd, aes(Actual, Predicted, fill = rate)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Freq, 100 * rate)), size = 5) +
  scale_fill_gradient(low = "#F1F5F4", high = "#1B7F79", limits = c(0, 1)) +
  coord_equal() +
  labs(title = paste0("Out-of-fold confusion matrix - ", sum(cmx) - sum(diag(cmx)),
                      " errors in ", sum(cmx), " predictions"),
       subtitle = "Percentages are column-wise (share of each true class)",
       x = "Actual sex", y = "Predicted sex", fill = "Column share")
sv(p7, "fig07_confusion_matrix.png", w = 8, h = 6)

# ---- Figure 8: predicted probability distribution ------------------------
pp <- data.frame(prob = cv_prob[[best]], sex = model_df$sex)
p8 <- ggplot(pp, aes(prob, fill = sex)) +
  geom_histogram(binwidth = .04, colour = "white", alpha = .85) +
  geom_vline(xintercept = .5, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = psex) +
  labs(title = "Predicted probabilities are pushed hard to the extremes",
       subtitle = "Out-of-fold P(male); the dashed line is the 0.5 decision threshold",
       x = "Predicted probability of male", y = "Count", fill = "Actual sex")
sv(p8, "fig08_probability_separation.png")

# ---- Figure 9: threshold sweep -------------------------------------------
ths <- seq(.05, .95, by = .01)
sweep <- do.call(rbind, lapply(ths, function(t) {
  m <- metrics(model_df$sex, cv_prob[[best]], t)
  data.frame(threshold = t, metric = c("Accuracy", "Sensitivity", "Specificity", "F1"),
             value = c(m["accuracy"], m["sensitivity"], m["specificity"], m["f1"]))
}))
opt <- ths[which.max(sapply(ths, function(t)
  metrics(model_df$sex, cv_prob[[best]], t)["accuracy"]))]
cat("\n== Accuracy-optimal threshold ==", opt, "\n")
print(round(metrics(model_df$sex, cv_prob[[best]], opt), 4))
p9 <- ggplot(sweep, aes(threshold, value, colour = metric)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = opt, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c(Accuracy = "#1B7F79", Sensitivity = "#1565C0",
                                 Specificity = "#C2185B", F1 = "#E1701A")) +
  labs(title = "Performance is flat across a wide band of decision thresholds",
       subtitle = paste0("Accuracy peaks at threshold = ", opt,
                         "; the default 0.5 costs only 0.3 pp of accuracy"),
       x = "Decision threshold on P(male)", y = "Score", colour = NULL)
sv(p9, "fig09_threshold_sweep.png")

# ---- Figure 10: calibration curve ---------------------------------------
bins <- cut(cv_prob[[best]], breaks = seq(0, 1, by = .1), include.lowest = TRUE)
cal <- data.frame(mid = tapply(cv_prob[[best]], bins, mean),
                  obs = tapply(model_df$sex == "male", bins, mean),
                  n   = as.integer(table(bins)))
cal <- cal[!is.na(cal$mid), ]
cat("\n== Calibration table ==\n"); print(round(cal, 3))
p10 <- ggplot(cal, aes(mid, obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_line(colour = "#1B7F79", linewidth = 1) +
  geom_point(aes(size = n), colour = "#1B7F79") +
  scale_size_continuous(range = c(2, 8)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "The model is well calibrated, not merely well ranked",
       subtitle = paste0("Observed male rate vs predicted probability, decile bins; Brier = ",
                         round(cv_tab[best, "brier"], 4)),
       x = "Mean predicted probability", y = "Observed proportion male", size = "n in bin")
sv(p10, "fig10_calibration.png", w = 8.4, h = 6.2)

# ==========================================================================
# 4. DIAGNOSTICS
# ==========================================================================
cat("\n\n########## 4. DIAGNOSTICS ##########\n")

# ---- Figure 11: binned residual plot ------------------------------------
res <- residuals(final, type = "response")
fit <- fitted(final)
b   <- cut(fit, quantile(fit, probs = seq(0, 1, length.out = 21)), include.lowest = TRUE)
br  <- data.frame(x = tapply(fit, b, mean), y = tapply(res, b, mean),
                  se = tapply(res, b, function(z) 2 * sd(z) / sqrt(length(z))))
p11 <- ggplot(br, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = -se, ymax = se), alpha = .15, fill = "#1B7F79") +
  geom_point(size = 2.6, colour = "#1B7F79") +
  labs(title = "Binned residual plot - no systematic curvature",
       subtitle = "Mean response residual per 5% bin of fitted probability, with +/-2 SE band",
       x = "Mean fitted probability", y = "Mean residual")
sv(p11, "fig11_binned_residuals.png")

# ---- Figure 12: influence diagnostics ------------------------------------
inf <- data.frame(obs = seq_len(nrow(model_df)), cook = cooks.distance(final),
                  lev = hatvalues(final), std = rstandard(final),
                  species = model_df$species)
cat("\n== Most influential observations (Cook's D) ==\n")
print(head(inf[order(-inf$cook), c("obs", "cook", "lev", "std", "species")], 6))
cat("\nMax Cook's D =", round(max(inf$cook), 4), " (rule of thumb cut-off 1)\n")
p12 <- ggplot(inf, aes(lev, std, size = cook, colour = species)) +
  geom_hline(yintercept = c(-2, 0, 2), linetype = c("dashed", "solid", "dashed"),
             colour = "grey60") +
  geom_point(alpha = .75) +
  scale_colour_manual(values = pal) + scale_size_continuous(range = c(1, 7)) +
  labs(title = "Influence diagnostics: leverage vs standardised residual",
       subtitle = "Point size is Cook's distance; no observation exceeds D = 0.5",
       x = "Leverage (hat value)", y = "Standardised deviance residual",
       colour = "Species", size = "Cook's D")
sv(p12, "fig12_influence.png", w = 9, h = 6)

# ---- Secondary regression model + classic residual diagnostics -----------
cat("\n== Secondary linear model: body mass ~ flipper + species + sex ==\n")
lm_fit <- lm(body_mass_g ~ flipper_length_mm + species + sex, data = model_df)
print(summary(lm_fit))
cat("\n== 10-fold CV RMSE / MAE / R2 for the linear model ==\n")
pred_lm <- numeric(nrow(model_df))
for (k in sort(unique(folds))) {
  f <- lm(body_mass_g ~ flipper_length_mm + species + sex, data = model_df[folds != k, ])
  pred_lm[folds == k] <- predict(f, newdata = model_df[folds == k, ])
}
e <- model_df$body_mass_g - pred_lm
cat("RMSE =", round(sqrt(mean(e ^ 2)), 1), " MAE =", round(mean(abs(e)), 1),
    " CV R2 =", round(1 - sum(e ^ 2) / sum((model_df$body_mass_g -
                                            mean(model_df$body_mass_g)) ^ 2), 4), "\n")
cat("\nShapiro-Wilk on lm residuals:\n"); print(shapiro.test(residuals(lm_fit)))

png(file.path(OUT, "fig13_lm_diagnostics.png"), width = 1700, height = 1300, res = 180)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1.5))
plot(lm_fit, col = "#1B7F79", pch = 19, cex = .5)
dev.off()

# ---- Optimisation: stepwise AIC on the interaction model -----------------
cat("\n== Stepwise AIC search (M4 interaction model as scope) ==\n")
step_fit <- stepAIC(glm(forms$M4_spp_inter, data = model_df, family = binomial()),
                    direction = "both", trace = FALSE)
print(formula(step_fit))
cat("AIC:", round(AIC(step_fit), 2), "vs champion", round(AIC(final), 2), "\n")
cv_step <- cv_glm(formula(step_fit), model_df, folds)
cat("\n== CV metrics of the stepwise model ==\n")
print(round(metrics(model_df$sex, cv_step), 4))

# ---- Repeated CV for a stability estimate --------------------------------
cat("\n== 10 x 10-fold repeated CV accuracy of the champion ==\n")
rep_acc <- sapply(1:10, function(i) {
  set.seed(100 + i)
  fl <- make_folds(model_df$sex, 10)
  p  <- cv_glm(forms[[best]], model_df, fl)
  mean((p >= .5) == (model_df$sex == "male"))
})
print(round(rep_acc, 4))
cat("mean =", round(mean(rep_acc), 4), " sd =", round(sd(rep_acc), 4),
    " range =", round(range(rep_acc), 4), "\n")

# ---- Prediction on the 9 held-out sex-unknown penguins -------------------
unknown <- df[is.na(df$sex), ]
unknown$p_male <- predict(final, newdata = unknown, type = "response")
cat("\n== Predictions for the 9 penguins whose sex was never recorded ==\n")
up <- unknown[, c("species", "bill_length_mm", "bill_depth_mm",
                  "flipper_length_mm", "body_mass_g", "p_male")]
up$p_male <- round(up$p_male, 3)
print(up, row.names = FALSE)

cat("\n== Environment ==\n")
cat(R.version.string, "| ggplot2", as.character(packageVersion("ggplot2")),
    "| MASS", as.character(packageVersion("MASS")), "\n")
