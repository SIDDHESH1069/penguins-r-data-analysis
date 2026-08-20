# ==========================================================================
# Week 2 - Data Visualization and Insight Communication using R
# Dataset: Palmer Archipelago (Antarctica) Penguins (CC0)
# Author : Siddhesh Mhatre
# ==========================================================================
library(ggplot2)
library(lattice)

set.seed(42)
OUT <- "/tmp/w2/figs"
sv <- function(p, file, w = 9, h = 5.6, dpi = 200) {
  ggsave(file.path(OUT, file), p, width = w, height = h, dpi = dpi, bg = "white")
}

# ---- 1. Load the cleaned data (Week 1 pipeline, re-applied) ---------------
penguins <- read.csv("/tmp/w2/penguins.csv", stringsAsFactors = FALSE)
penguins[penguins == ""] <- NA
df <- penguins[!is.na(penguins$bill_length_mm), ]          # drop 2 empty records
df$sex[is.na(df$sex)] <- "unrecorded"                      # keep 9 measured rows
df$species <- factor(df$species, levels = c("Adelie", "Chinstrap", "Gentoo"))
df$island  <- factor(df$island)
df$sex     <- factor(df$sex, levels = c("female", "male", "unrecorded"))
df$year    <- factor(df$year)

cat("== Dimensions ==\n"); print(dim(df))
cat("\n== str() ==\n");    str(df)
cat("\n== summary() ==\n"); print(summary(df))

cat("\n== Mean body mass by species and sex (g) ==\n")
print(round(tapply(df$body_mass_g, list(df$species, df$sex), mean, na.rm = TRUE), 1))

cat("\n== Counts: species x island ==\n")
print(table(df$species, df$island))

# ---- Shared theme --------------------------------------------------------
pal <- c(Adelie = "#E1701A", Chinstrap = "#7B3FA0", Gentoo = "#1B7F79")
theme_pg <- theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(colour = "grey35", size = 11),
        plot.caption  = element_text(colour = "grey45", size = 9),
        panel.grid.minor = element_blank(),
        legend.position = "top")
theme_set(theme_pg)

# ---- Figure 1: Bar chart - counts by species and island ------------------
p1 <- ggplot(df, aes(island, fill = species)) +
  geom_bar(width = .7, colour = "white") +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_stack(vjust = .5), colour = "white", size = 4) +
  scale_fill_manual(values = pal) +
  labs(title = "Where each penguin species was observed",
       subtitle = "Stacked counts of 342 observations across three islands",
       x = "Island", y = "Number of penguins", fill = "Species",
       caption = "Data: palmerpenguins (CC0)")
sv(p1, "fig01_species_island_bar.png")

# ---- Figure 2: Histogram + density of flipper length --------------------
p2 <- ggplot(df, aes(flipper_length_mm, fill = species)) +
  geom_histogram(binwidth = 3, colour = "white", alpha = .85, position = "identity") +
  geom_vline(xintercept = mean(df$flipper_length_mm), linetype = "dashed", colour = "grey30") +
  annotate("text", x = mean(df$flipper_length_mm) + 1, y = 26,
           label = paste0("overall mean = ", round(mean(df$flipper_length_mm), 1), " mm"),
           hjust = 0, size = 3.4, colour = "grey30") +
  scale_fill_manual(values = pal) +
  labs(title = "Flipper length is bimodal - and species explains it",
       subtitle = "Histogram (3 mm bins) of flipper length, filled by species",
       x = "Flipper length (mm)", y = "Count", fill = "Species")
sv(p2, "fig02_flipper_histogram.png")

# ---- Figure 3: Density by species ---------------------------------------
p3 <- ggplot(df, aes(body_mass_g, fill = species, colour = species)) +
  geom_density(alpha = .35, linewidth = .8) +
  scale_fill_manual(values = pal) + scale_colour_manual(values = pal) +
  labs(title = "Body-mass distributions barely overlap for Gentoo",
       subtitle = "Kernel density estimates of body mass by species",
       x = "Body mass (g)", y = "Density", fill = "Species", colour = "Species")
sv(p3, "fig03_mass_density.png")

# ---- Figure 4: Boxplot + jitter, mass by species and sex ----------------
p4 <- ggplot(subset(df, sex != "unrecorded"),
             aes(species, body_mass_g, fill = species)) +
  geom_boxplot(outlier.shape = NA, alpha = .55, width = .6) +
  geom_jitter(width = .16, alpha = .5, size = 1.4, aes(colour = species)) +
  facet_wrap(~ sex) +
  scale_fill_manual(values = pal) + scale_colour_manual(values = pal) +
  labs(title = "Males are consistently heavier within every species",
       subtitle = "Boxplots with raw observations, faceted by sex",
       x = NULL, y = "Body mass (g)") +
  theme(legend.position = "none")
sv(p4, "fig04_mass_box_sex.png")

# ---- Figure 5: Scatter with per-species regression lines ----------------
p5 <- ggplot(df, aes(flipper_length_mm, body_mass_g, colour = species)) +
  geom_point(size = 2, alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = .9) +
  scale_colour_manual(values = pal) +
  labs(title = "Flipper length predicts body mass in every species",
       subtitle = paste0("Linear fits per species; pooled Pearson r = ",
                         round(cor(df$flipper_length_mm, df$body_mass_g), 3)),
       x = "Flipper length (mm)", y = "Body mass (g)", colour = "Species")
sv(p5, "fig05_flipper_mass_scatter.png")

# ---- Figure 6: Simpson's paradox - bill depth vs length ----------------
r_pool <- cor(df$bill_length_mm, df$bill_depth_mm)
p6 <- ggplot(df, aes(bill_length_mm, bill_depth_mm, colour = species)) +
  geom_point(size = 2, alpha = .8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .9) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey25",
              linetype = "dashed", aes(group = 1)) +
  scale_colour_manual(values = pal) +
  labs(title = "Simpson's paradox: the pooled trend reverses within species",
       subtitle = paste0("Dashed grey = pooled fit (r = ", round(r_pool, 2),
                         "); coloured = within-species fits (all positive)"),
       x = "Bill length (mm)", y = "Bill depth (mm)", colour = "Species")
sv(p6, "fig06_simpsons_paradox.png")

# ---- Figure 7: Line chart - mean measurements over study years ---------
agg <- aggregate(cbind(body_mass_g, flipper_length_mm, bill_length_mm) ~ species + year,
                 data = df, FUN = mean)
cat("\n== Mean measurements by species and year ==\n")
print(transform(agg, body_mass_g = round(body_mass_g, 0),
                flipper_length_mm = round(flipper_length_mm, 1),
                bill_length_mm = round(bill_length_mm, 1)))
p7 <- ggplot(agg, aes(as.integer(as.character(year)), body_mass_g,
                      colour = species, group = species)) +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(breaks = 2007:2009) +
  labs(title = "Mean body mass is stable across the three breeding seasons",
       subtitle = "Annual species means, 2007-2009",
       x = "Study year", y = "Mean body mass (g)", colour = "Species")
sv(p7, "fig07_year_line.png")

# ---- Figure 8: Heat map of correlation matrix --------------------------
num <- df[, c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")]
cm  <- round(cor(num), 2)
cat("\n== Correlation matrix ==\n"); print(cm)
cml <- data.frame(as.table(cm)); names(cml) <- c("v1", "v2", "r")
p8 <- ggplot(cml, aes(v1, v2, fill = r)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  labs(title = "Correlation structure of the four body measurements",
       subtitle = "Pearson coefficients; bill depth is negatively correlated with the rest",
       x = NULL, y = NULL, fill = "r") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
sv(p8, "fig08_corr_heatmap.png", w = 8.6, h = 6.4)

# ---- Figure 9: lattice xyplot matrix ----------------------------------
png(file.path(OUT, "fig09_lattice_splom.png"), width = 1700, height = 1450, res = 190)
print(splom(~ num, groups = df$species, pch = 19, cex = .45,
            col = pal[levels(df$species)],
            key = list(space = "top", columns = 3,
                       points = list(pch = 19, col = pal[levels(df$species)]),
                       text = list(levels(df$species))),
            main = "Pairwise scatterplot matrix (lattice::splom)"))
dev.off()

# ---- Figure 10: base R mosaic of species x island ---------------------
png(file.path(OUT, "fig10_base_mosaic.png"), width = 1700, height = 1050, res = 190)
par(mar = c(4, 4, 4, 2))
mosaicplot(table(df$island, df$species), color = pal, las = 1,
           main = "Species composition of each island (base R mosaicplot)",
           xlab = "Island", ylab = "Species")
dev.off()

# ---- Statistical confirmation ----------------------------------------
cat("\n== ANOVA: body mass ~ species ==\n")
print(summary(aov(body_mass_g ~ species, data = df)))
cat("\n== Linear model: body mass ~ flipper length ==\n")
print(summary(lm(body_mass_g ~ flipper_length_mm, data = df))$coefficients)
cat("\nR-squared:", round(summary(lm(body_mass_g ~ flipper_length_mm, data = df))$r.squared, 3), "\n")
cat("\n== Within-species bill correlations ==\n")
print(round(sapply(split(df, df$species), function(d) cor(d$bill_length_mm, d$bill_depth_mm)), 3))
cat("\nPooled bill correlation:", round(r_pool, 3), "\n")
cat("\n== sessionInfo (abridged) ==\n")
cat(R.version.string, "| ggplot2", as.character(packageVersion("ggplot2")),
    "| lattice", as.character(packageVersion("lattice")), "\n")
