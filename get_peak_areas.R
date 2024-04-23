library(dplyr)
library(baseline)
setwd("~/Documents/Research - Others/Shu MYB99")
# metadata <- read.csv("MYB99_metadata_FINAL.csv")
# metadata <- metadata %>% mutate(sample = sample2)
# metadata <- metadata %>% select(id, sample)
# metadata <- metadata %>% mutate(expression = case_when(sample == "2099" ~ 5,
#                                                        sample == "4899" ~ 4,
#                                                        sample == "991W" ~ -1,
#                                                        sample == "995W" ~ -3,
#                                                        sample == "9999SR" ~ -4,
#                                                        sample == "4M" ~ 1,
#                                                        sample == "COL" ~ 3,
#                                                        .default = NA))
# metadata$expression_positive <- ifelse(metadata$expression > 0, "Positive", "Negative")
# metadata$expression_positive[is.na(metadata$expression_positive)] <- "Unknown"
# metadata$has1605 <- "No_1605"
# metadata$has1605[metadata$sample %in% c("2099", "4899", "COL", "991W", "995W")] <- "Yes_1605"
# 
# 
# spectra <- read.csv("MYB99_spectra_filtered_FINAL.csv")
# 
# isolation_forest_dodgy_obs <- c(55, 63, 122, 153, 161, 138, 139, 157)
# spectra <- spectra[, -1]
# wavenumbers <- strip_non_numbers(names(spectra))
# 
# exclude_obs <- c(1, 7, 9, 12, 53, 55, 63, 65, 68, 97, 100, 112, 117, 119, 122)
# exclude_obs <- c(exclude_obs, isolation_forest_dodgy_obs) %>% unique %>% as.numeric %>% sort
# metadata <- metadata[-exclude_obs, ]
# spectra <- spectra[-exclude_obs, ]
# 
# ############################################################
# ############################################################
# prepro_method <- "MPF"
# ############################################################
# ############################################################
# 
# if(prepro_method == "EMSC"){
#   spectra <- EMSC::EMSC(spectra, degree = 2)$corrected %>% data.frame()
# } else if(prepro_method == "MPF"){
#   spectra <- spectra %>%
#     as.matrix %>%
#     baseline(method = "modpolyfit", degree = 3) %>%
#     getCorrected %>%
#     t %>%
#     scale %>%
#     t %>%
#     data.frame %>%
#     setNames(wavenumbers)
# 
# }
# 
# all_data <- data.frame(metadata, spectra)
# all_data <- all_data %>% filter(!sample %in% c("991", "995", "MLP8"))
# all_data$sample <- all_data$sample %>% str_replace_all("W", "")
# all_data <- all_data %>% arrange(expression)

#saveRDS(all_data, "Shi MYB99 Data.RDS")

all_data <- readRDS("Shi MYB99 Data.RDS")
spectra <- all_data[-c(1:7)]
metadata <- all_data[1:7]
# write.csv(spectra, "Shi MYB99 spectra.csv", row.names = F)
# write.csv(metadata, "Shi MYB99 metadata.csv", row.names = F)

plotting_data <- unique(metadata[, c("sample", "expression")])
plotting_data$plot_pch <- 1:nrow(plotting_data)
plotting_data$plot_col <- paste0(c(rev(scico::scico(7, palette = "batlow")), rep("#000000", 4)), "CC")
plotting_data$plot_col_exp <- paste0(c(rev(scico::scico(11, palette = "vik"))[c(2, 3, 5, 7, 9, 10, 11)], rep("#000000", 4)), "CC")
plotting_data <- plotting_data %>% left_join(all_data %>% count(sample), by = "sample")
#plotting_data$plot_pch_alt <- c(0:2, 15:18, 3, 4, 8, 11)
plotting_data$plot_pch_alt <- c(18, 17, 16, 16, 17, 18, 15, 3, 4, 8, 11)

subset_myb <- function(samples){
  all_data_subset <<- all_data %>% filter(sample %in% samples)
  metadata_subset <<- all_data_subset[1:5]
  spectra_subset <<- all_data_subset[-c(1:5)]
  plotting_data_subset <<- plotting_data %>% filter(sample %in% samples)
}

# all_data %>% select(-id, -expression) %>% group_by(sample) %>% summarise_all("mean")

############################ Principal component analysis



plot_pca_scores <- function(apca, plot_pch = "plot_pch", plot_col = "plot_col", plot_cex, ...){
  pc1 <- 1
  pc2 <- 2
  par(mfrow = c(1, 1), mar = c(5.6, 4.6, 1.1, 1.1), pty = "s")
  plot(apca$x[, pc1], apca$x[, pc2], type = "n", axes = F, xlab = "", ylab = "")
  abline(h = 0, lty = 1, col = "darkgrey")
  abline(v = 0, lty = 1, col = "darkgrey")
  #centroids <- apca$x[, 1:2] %>% bind_cols(sample = metadata_subset$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
  
  for(i in 1:nrow(plotting_data_subset)){
    points(apca$x[metadata_subset$sample == plotting_data_subset$sample[i], pc1], apca$x[metadata_subset$sample == plotting_data_subset$sample[i], pc2], col = plotting_data_subset[[plot_col]][i], pch = plotting_data_subset[[plot_pch]][i], cex = plot_cex, lwd = 3)
  }
  
  legend(legend = plotting_data_subset$sample, pch = plotting_data_subset[[plot_pch]], col = plotting_data_subset[[plot_col]], cex = 1.5, ...)
  axis(1, cex.axis = 1.5, tck = 0.01)
  axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
  mtext(side = 1, line = 4, at = inRange(apca$x[, pc1], 0.5), cex = 2, text = paste0("PC", pc1, " scores (", pca_importance(apca, pc1), "%)"))
  mtext(side = 2, line = 3, at = inRange(apca$x[, pc2], 0.5), cex = 2, text = paste0("PC", pc2, " scores (", pca_importance(apca, pc2), "%)"))
  box()
}

plot_pca_loadings <- function(apca, pc1 = 1, pc2 = 2){
  apca_wns <- apca$rotation %>% row.names %>% strip_non_numbers
  mfrow(2, 1)
  par(pty = "m")
  plot(apca_wns, apca$rotation[, pc1], type = "n", xlim = rev(range(apca_wns)), axes = F,
       ylab = c(paste0("PC", pc1, " (", pca_importance(apca, pc1), "%)")),
       xlab = expression(paste("Wavenumber cm"^"-1"*"")),
       cex.lab = 1.5)
  abline(h = 0, col = "darkgrey")
  abline(v = c(1514, 1605, 1169), col = "red")
  lines(apca_wns, apca$rotation[, pc1], lwd = 2)
  internal_axis_ticks(1:2, 1:2, cex.axis = 1.2)
  plot(apca_wns, apca$rotation[, pc2], type = "n", xlim = rev(range(apca_wns)), axes = F,
       ylab = c(paste0("PC", pc2, " (", pca_importance(apca, pc2), "%)")),
       xlab = expression(paste("Wavenumber cm"^"-1"*"")),
       cex.lab = 1.5)
  abline(h = 0, col = "darkgrey")
  abline(v = c(1514, 1605, 1169), col = "red")
  lines(apca_wns, apca$rotation[, pc2], lwd = 2)
  internal_axis_ticks(1:2, 1:2, cex.axis = 1.2)
}

plot_wn_box <- function(tgt_wn, denom_wn = NULL, notch = T){
  wavenumbers <- strip_non_numbers(names(spectra_subset))
  abs_wn <- apply(spectra_subset[(nearest(wavenumbers, tgt_wn)$index - 1):(nearest(wavenumbers, tgt_wn)$index + 1)], 1, max)
  if(!is.null(denom_wn)){
    abs_wn2 <- apply(spectra_subset[(nearest(wavenumbers, denom_wn)$index - 1):(nearest(wavenumbers, denom_wn)$index + 1)], 1, max)
    abs_wn <- abs_wn / abs_wn2
    my_ylab <- paste0(tgt_wn, " cm-\u00b9 / ", denom_wn, " cm-\u00b9")
  } else {
    my_ylab <- paste0("Relative absorbance @ ", tgt_wn, " cm-\u00b9")
  }
  bp <- data.frame(sample = factor(metadata_subset$sample, levels = plotting_data_subset$sample), 
                   absorbance = abs_wn)
  plot(absorbance ~ sample, bp, col = plotting_data_subset$plot_col_exp, outline = F, pars = list(whisklty = 1, staplelwd = 0), las = 1, axes = F,
       xlab = "",
       ylab = my_ylab,
       ylim = range(absorbance),
       cex.lab = 1.5,
       notch = notch)
  stripchart(absorbance ~ sample, bp, vertical = T, method = "jitter", add = T, pch = 1, cex = 2, col = "#00000055", lwd = 2)
  internal_axis_ticks(2, 2, cex.axis = 1.2)
  mtext(1, 2, at = 1:nrow(plotting_data_subset), text = plotting_data_subset$sample, cex = 1.2)
}

plot_wn_against <- function(wn_a, wn_b = NULL, plot_pch = "plot_pch", plot_col = "plot_col", plot_cex, ...){
  wavenumbers <- strip_non_numbers(names(spectra_subset))
  x <- apply(spectra_subset[(nearest(wavenumbers, wn_a)$index - 1):(nearest(wavenumbers, wn_a)$index + 1)], 1, max)
  y <- apply(spectra_subset[(nearest(wavenumbers, wn_b)$index - 1):(nearest(wavenumbers, wn_b)$index + 1)], 1, max)
  
  ord <- data.frame(x = x, y = y, sample = factor(metadata_subset$sample, levels = plotting_data_subset$sample))
  
  plot(y ~ x, ord, col = plotting_data_subset$plot_col_exp, las = 1, axes = F,
       xlab = paste0("Relative absorbance @ ", wn_a, " cm-\u00b9"),
       ylab = paste0("Relative absorbance @ ", wn_b, " cm-\u00b9"),
       cex.lab = 1.5,
       type = "n")
  for(i in 1:nrow(plotting_data_subset)){
    points(ord[metadata_subset$sample == plotting_data_subset$sample[i], 1], 
           ord[metadata_subset$sample == plotting_data_subset$sample[i], 2], 
           col = plotting_data_subset[[plot_col]][i], 
           pch = plotting_data_subset[[plot_pch]][i], 
           cex = plot_cex, 
           lwd = 3)
  }
  legend(legend = plotting_data_subset$sample, pch = plotting_data_subset[[plot_pch]], col = plotting_data_subset[[plot_col]], cex = 1.5, ...)
  internal_axis_ticks(1:2, 1:2, cex.axis = 1.2)
}

plot_wn_against_alt <- function(wn_a, wn_b = NULL, plot_col = "plot_col", plot_pch = rep(16, 6), plot_cex, ...){
  wavenumbers <- strip_non_numbers(names(spectra_subset))
  x <- apply(spectra_subset[(nearest(wavenumbers, wn_a)$index - 1):(nearest(wavenumbers, wn_a)$index + 1)], 1, max)
  y <- apply(spectra_subset[(nearest(wavenumbers, wn_b)$index - 1):(nearest(wavenumbers, wn_b)$index + 1)], 1, max)
  
  ord <- data.frame(x = x, y = y, sample = factor(metadata_subset$sample, levels = plotting_data_subset$sample))
  
  plot(y ~ x, ord, col = plotting_data_subset$plot_col_exp, las = 1, axes = F,
       xlab = paste0("Relative absorbance @ ", wn_a, " cm-\u00b9"),
       ylab = paste0("Relative absorbance @ ", wn_b, " cm-\u00b9"),
       cex.lab = 1.5,
       type = "n")
  for(i in 1:nrow(plotting_data_subset)){
    points(ord[metadata_subset$sample == plotting_data_subset$sample[i], 1], 
           ord[metadata_subset$sample == plotting_data_subset$sample[i], 2], 
           col = plotting_data_subset[[plot_col]][i], 
           pch = plot_pch[i], 
           cex = plot_cex[i], 
           lwd = 3)
  }
  express <- plotting_data_subset$expression
  legend_symbols <- ifelse(express < 0, "-", "+")
  legend_text <- character(length(express))
  for(i in 1:length(express)){
    legend_text[i] <- paste0(rep(legend_symbols[i], abs(express[i])), collapse = "")
  }
  
  legend(legend = legend_text, pch = 16, col = plotting_data_subset[[plot_col]], cex = 2, ...)
  internal_axis_ticks(1:2, 1:2, cex.axis = 1.2)
}

plot_wn_against_expression <- function(wn_a, wn_b = NULL, plot_col = "plot_col", plot_cex, ...){
  wavenumbers <- strip_non_numbers(names(spectra_subset))
  x <- apply(spectra_subset[(nearest(wavenumbers, wn_a)$index - 1):(nearest(wavenumbers, wn_a)$index + 1)], 1, max)
  y <- apply(spectra_subset[(nearest(wavenumbers, wn_b)$index - 1):(nearest(wavenumbers, wn_b)$index + 1)], 1, max)
  
  ord <- data.frame(x = x, y = y, sample = factor(metadata_subset$sample, levels = plotting_data_subset$sample), expression = metadata_subset$expression) %>% mutate(ratio = x / y)
  
  plot(ratio ~ expression, ord, col = plotting_data_subset$plot_col_exp, las = 1, axes = F,
       xlab = paste0("Expression level"),
       ylab = paste0(wn_a, " cm-\u00b9 / ", wn_b, " cm-\u00b9"),
       cex.lab = 1.5,
       type = "n")
  for(i in 1:nrow(plotting_data_subset)){
    points(ord[metadata_subset$sample == plotting_data_subset$sample[i], ]$expression, 
           ord[metadata_subset$sample == plotting_data_subset$sample[i], ]$ratio, 
           col = plotting_data_subset[[plot_col]][i], 
           pch = 16, 
           cex = plot_cex, 
           lwd = 3)
  }
  express <- plotting_data_subset$expression
  legend_symbols <- ifelse(express < 0, "-", "+")
  legend_text <- character(length(express))
  for(i in 1:length(express)){
    legend_text[i] <- paste0(rep(legend_symbols[i], abs(express[i])), collapse = "")
  }
  
  legend(legend = legend_text, pch = 16, col = plotting_data_subset[[plot_col]], cex = 2, ...)
  internal_axis_ticks(1:2, 1:2, cex.axis = 1.2)
}
source("~/Documents/Coding/r-useful-funs/pca_importance.R")
source("~/Documents/Coding/r-useful-funs/inRange.R")
source("~/Documents/Coding/r-useful-funs/strip_non_numbers.R")
source("~/Documents/Coding/r-useful-funs/nearest.R")
source("~/Documents/Coding/Base plot quick functions.R")

subset_myb(c("2099", "4899", "COL", "9999SR"))
apca <- prcomp(spectra_subset[, 1:440], center = T, scale. = F)
plot_pca_scores(apca, "plot_pch_alt", "plot_col_exp", x = "bottomleft", plot_cex = 4)
plot_pca_loadings(apca)
mfrow(1, 2)
plot_wn_box(1513)
plot_wn_box(1605)
mfrow(1, 1)
plot_wn_box(1513, 1605)
mfrow(1, 1)
plot_wn_against(1513, 1605, "plot_pch_alt", "plot_col_exp", plot_cex = 4, x = "bottomright")

subset_myb(c("2099", "4899", "COL", "9999SR", "991", "995", "MLP8"))
apca <- prcomp(spectra_subset[, 1:440], center = T, scale. = F)
plot_pca_scores(apca, "plot_pch_alt", "plot_col_exp", x = "bottomleft", plot_cex = 4)
plot_pca_loadings(apca)
mfrow(1, 2)
plot_wn_box(1513)
plot_wn_box(1605)
mfrow(1, 1)
plot_wn_box(1513, 1605)
mfrow(1, 1)
plot_wn_against(1513, 1605, "plot_pch_alt", "plot_col_exp", plot_cex = 4, x = "bottomright")

subset_myb(c("2099", "4899", "COL", "9999SR", "991", "995"))
apca <- prcomp(spectra_subset[, 1:440], center = T, scale. = F)
plot_pca_scores(apca, "plot_pch_alt", "plot_col_exp", x = "bottomleft", plot_cex = 4)
plot_pca_loadings(apca)
mfrow(1, 1)
plot_wn_against_alt(1513, 1605, plot_col = "plot_col_exp", plot_cex = abs(plotting_data_subset$expression) * 2, x = "bottomright", plot_pch = c("-", "-", "-", "+", "+", "+"))
plot_wn_against_expression(1513, 1605, plot_col = "plot_col_exp", x = "topright", plot_cex = 2)

subset_myb(plotting_data$sample)
apca <- prcomp(spectra_subset[, 1:440], center = T, scale. = F)
plot_pca_scores(apca, "plot_pch_alt", "plot_col_exp", x = "bottomleft", plot_cex = 4)
plot_pca_loadings(apca)
mfrow(1, 2)
plot_wn_box(1513)
plot_wn_box(1605)
mfrow(1, 1)
plot_wn_box(1513, 1605)
mfrow(1, 1)
plot_wn_against(1513, 1605, "plot_pch_alt", "plot_col_exp", plot_cex = 4, x = "bottomright")

source("~/Documents/Coding/r-spectral-analysis-tools/get_mean_spectra.R")
source("~/Documents/Coding/r-spectral-analysis-tools/plot_spectra_base.R")
source("~/Documents/Coding/r-useful-funs/dec_to_hex.R")
source("~/Documents/Coding/quick_legend_ALPHA.R")
source("~/Documents/Coding/find_peaks.R")

plot_pars()
tidy_spectra <- get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 7)

plotting_data_sp <- data.frame(group = unique(tidy_spectra$group))
plotting_data_sp <- plotting_data_sp %>% left_join(metadata[c("sample", "has1605")], by = c("group" = "sample")) %>% unique()
plotting_data_sp$col <- ifelse(plotting_data_sp$has1605 == "Yes_1605", 1, 2)
plotting_data_sp$col <- scico::scico(2, palette = "vik")[plotting_data_sp$col]
plotting_data_sp$lty <- 1
plotting_data_sp$lwd <- 2
plotting_data_sp$pch <- NA


quick_legend(plotting_data_sp, "bottomright", cex = 1.5)
plot_spectra_base(tidy_spectra, wn_hi = 1705, wn_lo = 1505, lwd = 2, col_vector = plotting_data_sp$col, alpha = 0.7)
plot_spectra_base(tidy_spectra, wn_hi = 1470, wn_lo = 1270, lwd = 2, col_vector = plotting_data_sp$col, alpha = 0.7)
plot_spectra_base(tidy_spectra, wn_hi = 1270, wn_lo = 1070, lwd = 2, col_vector = plotting_data_sp$col, alpha = 0.7)

####### finding peak areas in derivatives
mfrow(2, 1)

tidy_spectra_2099 <- get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 7) %>% filter(group == "2099")
wn_hi <- 1800
wn_lo <- 950

my_peaks <- find_peaks(tidy_spectra_2099$mean, tidy_spectra_2099$wavenumber, wn_hi = wn_hi, wn_lo = wn_lo, ndowns = 3, type = "trough")
plot_spectra_base(tidy_spectra_2099, wn_hi = wn_hi, wn_lo = wn_lo, col_vector = "#000000", point_size = 0.5)
#points(my_peaks$peak_wavenumber, my_peaks$abs, col = "red")
#points(my_peaks$start_wavenumber, my_peaks$start_absorbance, col = "green", pch = 3)
#points(my_peaks$end_wavenumber, my_peaks$end_absorbance, col = "blue", pch = 4)
abline(v = c(1605, 1514, 1377, 1169), col = "darkgrey")
text(my_peaks$peak_wavenumber, my_peaks$nudge, labels = round(my_peaks$peak_wavenumber, 0), cex = 0.5)

get_peak_areas <- function(x, wns, peak_data, trancated_areas = F){
  lapply(1:nrow(peak_data), function(n){
    section <- peak_data %>% select(peak_start, peak_end, peak_ind) %>% slice(n) %>% as.numeric()
    peak_dat <- x[section[1]:section[2]]
    return_dat <- seq(x[section[2]], x[section[1]], length.out = length(section[1]:section[2]))
    #### this will only work for derivative rather than raw spectrum because negative values
    #### could possibly fix by converting to absolute values and then seeing if return exceeds orig values
    #### or by refering to my_peaks to check for trough = derivative or peak = raw spectrum
    return_dat <- rev(return_dat)
    check_outside_boundary <- return_dat < peak_dat
    return_dat[check_outside_boundary] <- peak_dat[check_outside_boundary]
    return_dat <- rev(return_dat)
    peak_wns <- wns[section[1]:section[2]]
    return_wns <- rev(wns[section[1]:section[2]])
    poly_dat <- data.frame(x = c(peak_wns, return_wns), y = c(peak_dat, return_dat))
    area <- area::polygon_area(as.matrix(poly_dat))
    list(polygon_area = area, polygon_data = poly_dat)
  })
}

my_polygons <- get_peak_areas(tidy_spectra_2099$mean, tidy_spectra_2099$wavenumber, my_peaks)
plot_polygons <- lapply(my_polygons, "[[", 2)
polygon_areas <- sapply(my_polygons, "[[", 1)

plot_peak_polygons <- function(plot_polygons, col = "#FFFF0077"){
  for(i in 1:length(plot_polygons)){
    poly_dat <- plot_polygons[[i]]
    polygon(poly_dat$x, poly_dat$y, col = col)
  }
}

plot_peak_polygons(plot_polygons)




plot_peak_areas <- function(group, wn_hi = 1800, wn_lo = 950){
  group_spectra <- get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 7)
  group_spectra <- group_spectra[group_spectra$group == group, ]
  my_peaks <- find_peaks(group_spectra$mean, group_spectra$wavenumber, wn_hi = wn_hi, wn_lo = wn_lo, ndowns = 3, type = "trough")
  plot_spectra_base(group_spectra, wn_hi = wn_hi, wn_lo = wn_lo, col_vector = "#000000", point_size = 0.5)
  #points(my_peaks$peak_wavenumber, my_peaks$abs, col = "red")
  #points(my_peaks$start_wavenumber, my_peaks$start_absorbance, col = "green", pch = 3)
  #points(my_peaks$end_wavenumber, my_peaks$end_absorbance, col = "blue", pch = 4)
  abline(v = c(1605, 1514, 1377, 1169), col = "darkgrey")
  text(my_peaks$peak_wavenumber, my_peaks$nudge, labels = round(my_peaks$peak_wavenumber, 0), cex = 0.5)
  
  my_polygons <- get_peak_areas(group_spectra$mean, group_spectra$wavenumber, peak_data = my_peaks)
  plot_polygons <- lapply(my_polygons, "[[", 2)
  polygon_areas <- sapply(my_polygons, "[[", 1)
  
  plot_peak_polygons(plot_polygons)
}

mfrow(1, 1)
plot_peak_areas("2099")
plot_peak_areas("4899")
plot_peak_areas("COL")
plot_peak_areas("991")
plot_peak_areas("995")
plot_peak_areas("9999SR")
plot_peak_areas("MLP5")
plot_peak_areas("MLP8")
plot_peak_areas("4M")
plot_peak_areas("99C1")
plot_peak_areas("99C2")

tidy_spectra <- get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 7)

cor_spectra <- tidy_spectra %>% select(-c(upper, lower)) %>% filter(wavenumber < 1800 & wavenumber > 950) %>% tidyr::spread("wavenumber", "mean") %>% tibble::column_to_rownames("group")

wn_cor <- cor(cor_spectra, method = "spearman")
corrplot::corrplot(wn_cor)

get_wn_correlation <- function(wn_cor, wn){
  wn_cor[, nearest(colnames(wn_cor), wn)$index] %>% 
    data.frame() %>% 
    tibble::rownames_to_column("wavenumber") %>%
    setNames(c("wavenumber", "cor")) %>% 
    mutate(wavenumber = strip_non_numbers(wavenumber))
}

plot_wn_correlation <- function(wn_cor, base_plot = T) {
  if(!base_plot){
    wn_cor %>% 
      ggplot(aes(as.numeric(wavenumber), as.numeric(cor))) + geom_line() +
      scale_x_reverse() +
      scale_y_continuous(limits = c(-1, 1)) +
      labs(x = "Wavenumber cm-1", y = "Pearson's r") +
      theme_bw()
  } else {
    plot(cor ~ wavenumber, wn_cor,
         xlim = rev(range(wavenumber)),
         ylim = c(-1.05, 1.05),
         type = "n")
    abline(h = 0, col = "darkgrey")
    lines(cor ~ wavenumber, wn_cor)
  }
}

cor_1605 <- get_wn_correlation(wn_cor, 1605)
x <- cor_1605$cor; wns <- cor_1605$wavenumber
find_peaks(cor_1605$cor, cor_1605$wavenumber)

plot_wn_correlation()


get_wn_correlation(wn_cor, 1514)
get_wn_correlation(wn_cor, 1377)
get_wn_correlation(wn_cor, 1168)

wn_pca <- prcomp(wn_cor, center = T, scale = T)




tidy_spectra_4899 <- get_mean_spectra(spectra, metadata$sample, smoothed_derivative = T, w = 7) %>% filter(group == "4899")

my_peaks <- find_peaks(tidy_spectra_4899$mean, tidy_spectra_4899$wavenumber, ndowns = 3, type = "trough")
plot_spectra_base(tidy_spectra_4899, wn_hi = 1800, wn_lo = 950, col_vector = "#000000", point_size = 0.5)
points(my_peaks$peak_wavenumber, my_peaks$abs, col = "red")
points(my_peaks$start_wavenumber, my_peaks$start_absorbance, col = "green", pch = 3)
points(my_peaks$end_wavenumber, my_peaks$end_absorbance, col = "blue", pch = 4)
abline(v = c(1605, 1514, 1377, 1169), col = "darkgrey")



spec_diff_MLP8W <- colMeans(spectra[which(metadata$sample == "MLP8W"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff_MLP8 <- colMeans(spectra[which(metadata$sample == "MLP8"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff_9999SR <- colMeans(spectra[which(metadata$sample == "9999SR"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff_995W <- colMeans(spectra[which(metadata$sample == "995W"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff_995 <- colMeans(spectra[which(metadata$sample == "995"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff_991W <- colMeans(spectra[which(metadata$sample == "991W"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff_991 <- colMeans(spectra[which(metadata$sample == "991"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff_4899 <- colMeans(spectra[which(metadata$sample == "4899"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff_2099 <- colMeans(spectra[which(metadata$sample == "2099"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])


par(pty = "m", mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 0.6))
plot(wavenumbers, spec_diff_MLP8, xlim = c(1800, 948), 
     ylim = c(-0.5, 0.5),
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(v = c(1710, 1605, 1514, 1436, 1168), lty = 2, col = "#000000DD")
#lutzke_wavenumbers()
abline(v = c(1605, 1515), lty = 2, lwd = 2)
#lines(wavenumbers, spec_diff_MLP8, lwd = 2, col = plot_cols[1])
abline(h = 0, lty = 1, lwd = 2, col = plot_cols[8])
lines(wavenumbers, spec_diff_MLP8W, lwd = 2, col = plot_cols[2])
lines(wavenumbers, spec_diff_9999SR, lwd = 2, col = plot_cols[3])
#lines(wavenumbers, spec_diff_995W, lwd = 2, col = plot_cols[4])
#lines(wavenumbers, spec_diff_995, lwd = 2, col = plot_cols[5])
lines(wavenumbers, spec_diff_991W, lwd = 2, col = plot_cols[6])
#lines(wavenumbers, spec_diff_991, lwd = 2, col = plot_cols[7])
lines(wavenumbers, spec_diff_4899, lwd = 2, col = plot_cols[9])
lines(wavenumbers, spec_diff_2099, lwd = 2, col = plot_cols[10])
#points(wavenumbers[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)

axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, at = 0, label = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = 0, cex = 1.5, text = "Relative difference (sample - control)")
#mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Relative difference (sample - control)")
box()
# text(1710, 0.15, "***", cex = 2)
# text(1710, -0.19, "***", cex = 2)
# text(1605, 0.125, "***", cex = 2)
# text(1605, -0.23, "**", cex = 2)
# text(1514, 0.18, "***", cex = 2)
# text(1514, -0.14, "*", cex = 2)
# text(1436, 0.125, "NS", cex = 1.5)
# text(1436, -0.23, "***", cex = 2)
# text(1168, 0.09, "***", cex = 2)
# text(1168, -0.11, "NS", cex = 1.5)

# mtext(side = 3, line = 0.5, at = 1710, text = "C=O")
# mtext(side = 3, line = 0.5, at = 1605, text = "C=C")
# mtext(side = 3, line = 0.5, at = 1514, text = "C=C")
# mtext(side = 3, line = 0.5, at = 1436, text = "dC-H")
#axis(2)
#rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
#rect(900, 0, 940, -10, col = plot_cols[3], border = NA)
#text(1050, max(spec_diff), "Plot ref: OQZU")
#phil_wavenumbers()
legend("bottomright", legend = c("MLP8", "9999SR", "COL"), col = plot_cols[c(1, 3, 8)], lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
#legend("bottomright", legend = factor_levels, col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
text(965, -0.22, "Plot ref: R43X")
#abline(v = interesting_wavenumbers, lwd = 2, col = "#00000077")

par(pty = "m", mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 0.6))
plot(wavenumbers, spec_diff, xlim = c(1800, 948), 
     ylim = c(-0.23, 0.29),
     type = "n", axes = F, xlab = "", ylab = "")
rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(v = c(1710, 1605, 1514, 1436, 1168), lty = 2, col = "#000000DD")
lutzke_wavenumbers(cex = 0.8)
lines(wavenumbers, spec_diff, lwd = 2, col = plot_cols[4])
lines(wavenumbers, spec_diff2, lwd = 2, col = plot_cols[1])
lines(wavenumbers, spec_diff3, lwd = 2, col = plot_cols[2])
#points(wavenumbers[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
abline(h = 0, lty = 1, lwd = 2, col = plot_cols[3])
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, at = 0, label = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Relative difference (sample - control)")
box()
legend("topleft", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
text(965, -0.22, "Plot ref: V2BE")


# mtext(side = 3, line = 0.5, at = 1710, text = "C=O")
# mtext(side = 3, line = 0.5, at = 1605, text = "C=C")
# mtext(side = 3, line = 0.5, at = 1514, text = "C=C")
# mtext(side = 3, line = 0.5, at = 1436, text = "dC-H")
#axis(2)
#rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
#rect(900, 0, 940, -10, col = plot_cols[3], border = NA)
#text(1050, max(spec_diff), "Plot ref: OQZU")
#phil_wavenumbers()
legend("topright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
text(965, -0.22, "Plot ref: R43X")
#abline(v = interesting_wavenumbers, lwd = 2, col = "#00000077")

ez <- function(x){
  paste(x, collapse = ", ")
}

diff_4899 <- diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2)
diff_2099 <- diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2)
diff_COL <- diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)
diff_9999SR <- diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2)

sd(diff_4899)

(peaks_4899_full <- find_peaks(diff_4899, type = "t", minpeakheight = 0))
(peaks_2099_full <- find_peaks(diff_2099, type = "t", minpeakheight = 0))
(peaks_COL_full <- find_peaks(diff_COL, type = "t", minpeakheight = 0))
(peaks_9999SR_full <- find_peaks(diff_9999SR, type = "t", minpeakheight = 0))

peaks_4899_full[nearest(peaks_4899_full$peak_wavenumber, 3026)$index, ]
peaks_2099_full[nearest(peaks_4899_full$peak_wavenumber, 3026)$index, ]
peaks_COL_full[nearest(peaks_4899_full$peak_wavenumber, 3026)$index, ]
peaks_9999SR_full[nearest(peaks_4899_full$peak_wavenumber, 3026)$index, ]

#peaks in second derivatives
(peaks_4899 <- find_peaks(diff_4899[1:440], type = "t", minpeakheight = 0))
(peaks_2099 <- find_peaks(diff_2099[1:440], type = "t", minpeakheight = 0))
(peaks_COL <- find_peaks(diff_COL[1:440], type = "t", minpeakheight = 0))
(peaks_9999SR <- find_peaks(diff_9999SR[1:440], type = "t", minpeakheight = 0))
peaks_4899$peak_wavenumber
peaks_2099$peak_wavenumber

all_peaks <- unique(c(peaks_4899$peak_wavenumber, peaks_2099$peak_wavenumber, peaks_COL$peak_wavenumber, peaks_9999SR$peak_wavenumber))
unique_peaks <- sort(all_peaks)
peak_present <- matrix(F, ncol = length(unique_peaks), nrow = 4)
colnames(peak_present) <- unique_peaks
rownames(peak_present) <- c("4899", "2099", "COL", "9999SR")
for(i in seq_along(peaks_4899$peak_wavenumber)){
  wn <- peaks_4899$peak_wavenumber[i]
  peak_present[1, which.min(abs(wn - unique_peaks))] <- T
}
for(i in seq_along(peaks_2099$peak_wavenumber)){
  wn <- peaks_2099$peak_wavenumber[i]
  peak_present[2, which.min(abs(wn - unique_peaks))] <- T
}
for(i in seq_along(peaks_COL$peak_wavenumber)){
  wn <- peaks_COL$peak_wavenumber[i]
  peak_present[3, which.min(abs(wn - unique_peaks))] <- T
}
for(i in seq_along(peaks_9999SR$peak_wavenumber)){
  wn <- peaks_9999SR$peak_wavenumber[i]
  peak_present[4, which.min(abs(wn - unique_peaks))] <- T
}
peak_present <- peak_present[c(2, 1, 3, 4), ]
peak_present2 <- lapply(data.frame(peak_present), as.character) %>% data.frame()
row.names(peak_present2) <- c("2099", "4899", "COL", "9999SR")
peak_present2 <- ifelse(peak_present2 == "TRUE", "*", "")
peak_present2 <- t(peak_present2) %>% data.frame()
row.names(peak_present2) <- stringr::str_replace_all(row.names(peak_present2), "X", "")
#write.csv(peak_present2, "Arabidopsis peaks present in second derivatives.csv")

abline(v = strip_non_numbers(names(shared_peaks)), col = "red")

lwd <- 9
plot(1:length(unique_peaks), seq(0, 4, length.out = length(unique_peaks)), xlim = c(length(unique_peaks), 1), type = "n", axes = F)
for(i in 1:length(unique_peaks)){
  if(peak_present[1, i] == T){
    arrows(i, 0, i, 1, length = 0, col = plot_cols[1], lwd = lwd)
  }
  if(peak_present[2, i] == T){
    arrows(i, 1, i, 2, length = 0, col = plot_cols[2], lwd = lwd)
  }
  if(peak_present[3, i] == T){
    arrows(i, 2, i, 3, length = 0, col = plot_cols[3], lwd = lwd)
  }
  if(peak_present[4, i] == T){
    arrows(i, 3, i, 4, length = 0, col = plot_cols[4], lwd = lwd)
  }
}

(shared_peaks <- which(colSums(peak_present) == 4))
(shared_peaks_MYB <- which(colSums(peak_present[1:3, ]) == 3))
(shared_peaks_MYB[!shared_peaks_MYB %in% shared_peaks])

expand_wavenumbers <- function(vector, n = 2, by = 1){
  expand <- seq(-n, n, by)
  do.call(c, lapply(expand, function(x){vector + x}))
}

fuzzy_wavenumbers <- expand_wavenumbers(round(strip_non_numbers(names(shared_peaks))) + 2)
fuzzy_wavenumbers_MYB <- expand_wavenumbers(round(strip_non_numbers(names(shared_peaks_MYB))) + 2)
lutzke$lutzke_wavenumbers %in% fuzzy_wavenumbers
lutzke$lutzke_wavenumbers %in% fuzzy_wavenumbers_MYB

lutzke$lutzke_wavenumbers[which(lutzke$lutzke_wavenumbers %in% fuzzy_wavenumbers)]

spectra2 <- spectra %>% bind_cols(sample = metadata$sample)
spectra2$sample <- factor(spectra2$sample, levels = rev(c("2099", "4899", "COL", "9999SR")))
for(i in 1:length(shared_peaks_MYB)){
  wn <- strip_non_numbers(names(shared_peaks_MYB[i]))
  wn_ind <- which.min(abs(wavenumbers - wn))
  png(paste0("boxplot_", strip_non_numbers(names(shared_peaks_MYB[i])) + 2, ".png"))
  boxplot(as.formula(spectra2[, wn_ind] ~ sample), spectra2, col = rev(plot_cols), notch = T, main = wn)
  dev.off()
}

aov_shared_peaks <- vector("list", length = length(shared_peaks_MYB))
for(i in 1:length(shared_peaks_MYB)){
  wn <- strip_non_numbers(names(shared_peaks_MYB[i])) + 2
  wn_ind <- which.min(abs(wavenumbers - wn))
  aov_shared_peaks[[i]] <- summary(aov(as.formula(spectra2[, wn_ind] ~ sample), spectra2))
}
names(aov_shared_peaks) <- shared_peaks_MYB %>% names() %>% as.numeric() %>% `+`(2) %>% as.character()
aov_shared_peaks

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 2.1, 0.6))
ylims <- c(-0.35, 0.025)
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 900), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(h = 0, lty = 1, col = "black")
lutzke_wavenumbers(cex = 0.8)
abline(h = -c(0, 0.1, 0.2, 0.3))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[1], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2) - 0.1, col = plot_cols[2], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - 0.2, col = plot_cols[3], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - 0.3, col = plot_cols[4], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
rect(948, 1, 800, -1, col = "white", border = NA)
text(940, 0, "4899", cex = 1.5, adj = c(0, 0.5))
text(940, -0.1, "2099", cex = 1.5, adj = c(0, 0.5))
text(940, -0.2, "COL", cex = 1.5, adj = c(0, 0.5))
text(940, -0.3, "9999SR", cex = 1.5, adj = c(0, 0.5))
#legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
text(930, -0.335, "Plot ref: KZZH")
box()

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 2.1, 0.6))
ylims <- c(-0.65, 0.025)
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 900), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
rect(strip_non_numbers(names(shared_peaks_MYB[!shared_peaks_MYB %in% shared_peaks])) + 4,
     -10,
     strip_non_numbers(names(shared_peaks_MYB[!shared_peaks_MYB %in% shared_peaks])) - 2,
     10,
     border = NA,
     col = "#FF000022")
rect(strip_non_numbers(names(shared_peaks_MYB[shared_peaks_MYB %in% shared_peaks])) + 4,
     -10,
     strip_non_numbers(names(shared_peaks_MYB[shared_peaks_MYB %in% shared_peaks])) - 2,
     10,
     border = NA,
     col = "#00000022")
#abline(v = strip_non_numbers(names(shared_peaks_MYB[!shared_peaks_MYB %in% shared_peaks])) + 2, lty = 3, lwd = 2, col = "darkgrey")
#abline(v = strip_non_numbers(names(shared_peaks_MYB[shared_peaks_MYB %in% shared_peaks])) + 2, lty = 5, lwd = 2, col = "darkgrey")
abline(h = -c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), col = "lightgrey")
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2) - 0, col = plot_cols[7], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2) - 0.1, col = plot_cols[6], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "991"), ]), differences = 2) - 0.2, col = plot_cols[5], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - 0.3, col = plot_cols[4], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "995"), ]), differences = 2) - 0.4, col = plot_cols[3], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - 0.5, col = plot_cols[2], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "MLP8"), ]), differences = 2) - 0.6, col = plot_cols[1], lwd = 3)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^" nd"*" derivative)")))
rect(948, 1, 800, -1, col = "white", border = NA)
text(940, 0, "2099", cex = 1.5, adj = c(0, 0.5))
text(940, -0.1, "4899", cex = 1.5, adj = c(0, 0.5))
text(940, -0.2, "991", cex = 1.5, adj = c(0, 0.5))
text(940, -0.3, "COL", cex = 1.5, adj = c(0, 0.5))
text(940, -0.4, "995", cex = 1.5, adj = c(0, 0.5))
text(940, -0.5, "9999SR", cex = 1.5, adj = c(0, 0.5))
text(940, -0.6, "MLP8", cex = 1.5, adj = c(0, 0.5))
#legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
#text(930, -0.335, "Plot ref: VWLR")
box()

spectra_diff <- t(spectra) |> diff(differences = 2) |> t() |> data.frame()
sum_9999SR <- spectra_diff[metadata$sample == "9999SR", ] |> colSums()
mean_9999SR <-  spectra_diff[metadata$sample == "9999SR", ] |> colMeans()

mean_MYB99_diff <-  spectra_diff[metadata$sample %in% c("2099", "4899", "COL"), ] |> colMeans()

plot(wavenumbers[2:1583], mean_MYB99_diff, type = "l", xlim = c(3100, 3000), ylim = c(-0.01, 0.01))

(mean_MYB99_peaks <- find_peaks(-mean_MYB99_diff, type = "p", minpeakheight = 0, nudge = 0.002))
#(mean_MYB99_peaks <- find_peaks(-mean_MYB99_diff, type = "p", minpeakheight = sd(mean_MYB99_diff)))
mean_MYB99_peaks$sd <- mean_MYB99_peaks$abs / sd(mean_MYB99_diff)
mean_MYB99_peaks$sd <- round(mean_MYB99_peaks$sd, 1)
mean_MYB99_peaks[, c(5, 8)]

common_wavenumbers <- c(954,
                        987,
                        1018,
                        1025,
                        1049,
                        1078,
                        1103,
                        1112,
                        1153,
                        1166,
                        1203,
                        1224,
                        1232,
                        1257,
                        1278,
                        1288,
                        1313,
                        1338,
                        1375,
                        1396,
                        1417,
                        1434,
                        1454,
                        1463,
                        1496,
                        1513,
                        1538,
                        1556,
                        1602,
                        1616,
                        1633,
                        1643,
                        1650,
                        1660,
                        1668,
                        1681,
                        1697,
                        1714,
                        1731,
                        1739,
                        1781,
                        1791
)

mean_MYB99_peaks <- mean_MYB99_peaks[str_detect(mean_MYB99_peaks$peak_wavenumber, rebus::or1(common_wavenumbers)), ]
range(mean_MYB99_peaks$abs)

par(mfrow = c(1, 1), mar = c(5.1, 5.1, 0.6, 0.6), pty = "m")
plot(wavenumbers[2:439], mean_MYB99_diff[1:438], 
     xlim = rev(range(wavenumbers[2:439])), 
     ylim = c(min(mean_MYB99_diff) - 0.005, max(mean_MYB99_diff)),
     xlab = "",
     ylab = "",
     type = "n", axes = F)
abline(h = 0, col = "black")
abline(h = -sd(mean_MYB99_diff), col = "lightgrey", lty = 2)
#lines(wavenumbers[2:439], mean_9999SR[1:438], lwd = 2)
lines(wavenumbers[2:439], mean_MYB99_diff[1:438], lwd = 2, col = "#00000099")
#lines(wavenumbers[2:439], diff_COL[1:438], lwd = 2, col = plot_cols[3])
#lines(wavenumbers[2:439], diff_4899[1:438], lwd = 2, col = plot_cols[2])
text(mean_MYB99_peaks$peak_wavenumber, -mean_MYB99_peaks$nudge, rep(c(letters, LETTERS), length.out = nrow(mean_MYB99_peaks)), cex = 0.8)
axis(1, cex.axis = 1.2, tck = 0.01)
mtext(side = 1, line = 3, at = midpoint(wavenumbers[1:440]), cex = 1.5, text = lab("wn"))
axis(2, cex.axis = 1.2, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = midpoint(mean_MYB99_diff), cex = 1.5, text = "Relative absorbance")
box()

par(mfrow = c(1, 1), mar = c(4.1, 5.1, 1.1, 1.1))
# plot(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1515)$index] ~ 
#        sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), 
plot(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1515)$index] ~ 
       sample, bind_cols(spectra, sample = factor(metadata$sample, levels = unique(metadata$sample))), 
     notch = T, 
     axes = F, 
     col = rev(plot_cols), 
     ylim = c(0.5, 0.9), 
     ylab = expression(paste("1605/1515 cm"^-1*"")), 
     cex.lab = 1.5, 
     cex.axis = 1.3, 
     axes = F, 
     boxwex = 0.6, 
     outline = F, 
     staplelty = 0, 
     whisklty = 1, 
     las = 1, 
     xlab = "")
stripchart(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1515)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = unique(metadata$sample))), method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
text(1:length(unique(metadata$sample)), 0.9, unname(table(metadata$sample)), cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:length(unique(metadata$sample)), text = unique(metadata$sample), cex = 1.5)
text(4, 0.6, "Plot ref: BXCC", cex = 1.2)

plot(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ 
       sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), 
     notch = T, 
     axes = F, 
     col = rev(plot_cols), 
     ylim = c(1.1, 1.6), 
     ylab = expression(paste("1605/1378 cm"^-1*"")), 
     cex.lab = 1.5, 
     cex.axis = 1.3, 
     axes = F, 
     boxwex = 0.6, 
     outline = F, 
     staplelty = 0, 
     whisklty = 1, 
     las = 1, 
     xlab = "")
stripchart(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
text(1:4, 1.6, unname(table(metadata$sample))[c(3, 4, 2, 1)], cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 1.5)
text(4, 1.1, "Plot ref: BXC2", cex = 1.2)

plot(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ 
       sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), 
     notch = T, 
     axes = F, 
     col = rev(plot_cols), 
     ylim = c(1.1, 1.6), 
     ylab = expression(paste("1605/1168 cm"^-1*"")), 
     cex.lab = 1.5, 
     cex.axis = 1.3, 
     axes = F, 
     boxwex = 0.6, 
     outline = F, 
     staplelty = 0, 
     whisklty = 1, 
     las = 1, 
     xlab = "")
stripchart(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
text(1:4, 1.6, unname(table(metadata$sample))[c(3, 4, 2, 1)], cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 1.5)
text(4, 1.1, "Plot ref: BXC1", cex = 1.2)

plot(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ 
       sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), 
     notch = T, 
     axes = F, 
     col = rev(plot_cols), 
     ylim = c(1.5, 2.2), 
     ylab = expression(paste("1515/1378 cm"^-1*"")), 
     cex.lab = 1.5, 
     cex.axis = 1.3, 
     axes = F, 
     boxwex = 0.6, 
     outline = F, 
     staplelty = 0, 
     whisklty = 1, 
     las = 1, 
     xlab = "")
stripchart(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
text(1:4, 2.2, unname(table(metadata$sample))[c(3, 4, 2, 1)], cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 1.5)
text(4, 1.5, "Plot ref: BXC3", cex = 1.2)
summary(aov(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099")))))

plot(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ 
       sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), 
     notch = T, 
     axes = F, 
     col = rev(plot_cols), 
     ylim = c(1.7, 2.3), 
     ylab = expression(paste("1515/1168 cm"^-1*"")), 
     cex.lab = 1.5, 
     cex.axis = 1.3, 
     axes = F, 
     boxwex = 0.6, 
     outline = F, 
     staplelty = 0, 
     whisklty = 1, 
     las = 1, 
     xlab = "")
stripchart(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
text(1:4, 2.3, unname(table(metadata$sample))[c(3, 4, 2, 1)], cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 1.5)
text(4, 1.7, "Plot ref: BXC4", cex = 1.2)
summary(aov(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099")))))

par(mfrow = c(1, 4), mar = c(4.1, 5.1, 1.1, 1.1))
plot(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ 
       sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), 
     notch = T, 
     axes = F, 
     col = rev(plot_cols), 
     ylim = c(1.1, 2.3), 
     ylab = expression(paste("1605/1378 cm"^-1*"")), 
     cex.lab = 1.5, 
     cex.axis = 1.3, 
     axes = F, 
     boxwex = 0.6, 
     outline = F, 
     staplelty = 0, 
     whisklty = 1, 
     las = 1, 
     xlab = "")
stripchart(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
#text(1:4, 1.6, unname(table(metadata$sample))[c(3, 4, 2, 1)], cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)


plot(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ 
       sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), 
     notch = T, 
     axes = F, 
     col = rev(plot_cols), 
     ylim = c(1.1, 2.3), 
     ylab = expression(paste("1605/1168 cm"^-1*"")), 
     cex.lab = 1.5, 
     cex.axis = 1.3, 
     axes = F, 
     boxwex = 0.6, 
     outline = F, 
     staplelty = 0, 
     whisklty = 1, 
     las = 1, 
     xlab = "")
stripchart(spectra[, nearest(wavenumbers, 1605)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
#text(1:4, 1.6, unname(table(metadata$sample))[c(3, 4, 2, 1)], cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)


plot(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ 
       sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), 
     notch = T, 
     axes = F, 
     col = rev(plot_cols), 
     ylim = c(1.1, 2.3), 
     ylab = expression(paste("1515/1378 cm"^-1*"")), 
     cex.lab = 1.5, 
     cex.axis = 1.3, 
     axes = F, 
     boxwex = 0.6, 
     outline = F, 
     staplelty = 0, 
     whisklty = 1, 
     las = 1, 
     xlab = "")
stripchart(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
#text(1:4, 2.2, unname(table(metadata$sample))[c(3, 4, 2, 1)], cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
summary(aov(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1378)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099")))))

plot(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ 
       sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), 
     notch = T, 
     axes = F, 
     col = rev(plot_cols), 
     ylim = c(1.1, 2.3), 
     ylab = expression(paste("1515/1168 cm"^-1*"")), 
     cex.lab = 1.5, 
     cex.axis = 1.3, 
     axes = F, 
     boxwex = 0.6, 
     outline = F, 
     staplelty = 0, 
     whisklty = 1, 
     las = 1, 
     xlab = "")
stripchart(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099"))), method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
#text(1:4, 2.3, unname(table(metadata$sample))[c(3, 4, 2, 1)], cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
text(2.5, 1.2, "Plot ref: BXCA", cex = 1.5)
summary(aov(spectra[, nearest(wavenumbers, 1515)$index] / spectra[, nearest(wavenumbers, 1168)$index] ~ sample, bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099")))))

spectra_cor <- cor(spectra[, 1:440], method = "kendall")
names(spectra)[1:440]
spectra_cor[295, ]

par(mfrow = c(1,1))
corrplot::corrplot(spectra_cor)

a <- c(5, 7, 34 , 4, 7, 1, 2345, 8, 3, 77, 9)
b <- a /2
cor.test(a, b)

bound_spectra <- bind_cols(spectra, sample = factor(metadata$sample, levels = c("9999SR", "COL", "4899", "2099")))

bound_spectra |> ggplot(aes(X1604.48, X1515.78)) + geom_point(size = 3) + geom_smooth(method = "lm")
summary(lm(X1515.78 ~ X1604.48, bound_spectra))

bound_spectra |> ggplot(aes(X1604.48, X1515.78, col = sample)) + geom_point(size = 3) + geom_smooth(method = "lm")
summary(lm(X1515.78 ~ X1604.48, bound_spectra))

plot(X1378.85 ~ X1515.78, bound_spectra, las = 1, type = "n", axes = F)
abline(lm(X1378.85 ~ X1515.78, bound_spectra))
points(X1378.85 ~ X1515.78, bound_spectra)
summary(lm(X1378.85 ~ X1515.78, bound_spectra))
internal_axis_ticks(1:2, 1:2)

bound_spectra |> ggplot(aes(X1515.78, X1168.65)) + geom_point(size = 3) + geom_smooth(method = "lm")
summary(lm(X1168.65 ~ X1515.78, bound_spectra))
bound_spectra |> ggplot(aes(X1515.78, X1168.65, col = sample)) + geom_point(size = 3) + geom_smooth(method = "lm")
summary(lm(X1168.65 ~ X1515.78, bound_spectra))
bound_spectra |> ggplot(aes(X1604.48, X1168.65)) + geom_point(size = 3) + geom_smooth(method = "lm")
summary(lm(X1168.65 ~ X1515.78, bound_spectra))
bound_spectra |> ggplot(aes(X1604.48, X1168.65, col = sample)) + geom_point(size = 3) + geom_smooth(method = "lm")
summary(lm(X1168.65 ~ X1515.78, bound_spectra))

bound_spectra |> ggplot(aes(X1604.48, X1378.85)) + geom_point() + geom_smooth(method = "lm")
summary(lm(X1378.85 ~ X1604.48, bound_spectra))
bound_spectra |> ggplot(aes(X1604.48, X1168.65)) + geom_point() + geom_smooth(method = "lm")
summary(lm(X1168.65 ~ X1604.48, bound_spectra))



par(mfrow = c(1, 1))
ylims <- c(-0.05, 0.025)
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 900), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(h = 0, lty = 1, col = "black")
lutzke_wavenumbers(cex = 0.8)
abline(h = -c(0, 0.1, 0.2, 0.3))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[1], lwd = 3)
#lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2) - 0.1, col = plot_cols[2], lwd = 2)
#lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - 0.2, col = plot_cols[3], lwd = 2)
#lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - 0.3, col = plot_cols[4], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(940, 0, "4899", cex = 1.5, adj = c(0, 0.5))
#text(940, -0.1, "2099", cex = 1.5, adj = c(0, 0.5))
#text(940, -0.2, "COL", cex = 1.5, adj = c(0, 0.5))
#text(940, -0.3, "9999SR", cex = 1.5, adj = c(0, 0.5))
#legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
#text(930, -0.335, "Plot ref: KZZH")
box()

peaks_4899 <- find_peaks(diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), minpeakdistance = 8, deriv = T)
peaks_4899 <- peaks_4899[peaks_4899$peak_type == "trough", ]
#peaks_4899 <- peaks_4899[peaks_4899$abs > 0, ]
#abline(v = peaks_4899$peak_wavenumber, col = "red")
peaks_2099 <- find_peaks(diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2), minpeakdistance = 8, deriv = T)
peaks_2099 <- peaks_2099[peaks_2099$peak_type == "trough", ]
#peaks_4899 <- peaks_4899[peaks_4899$abs > 0, ]
#abline(v = peaks_2099$peak_wavenumber, col = "red")
peaks_COL <- find_peaks(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), minpeakdistance = 8, deriv = T)
peaks_COL <- peaks_COL[peaks_COL$peak_type == "trough", ]

abline(v = wavenumbers[-c(1, length(wavenumbers))][intersect(peaks_4899$peak_ind, peaks_2099$peak_ind)], col = "red")

wavenumbers[-c(1, length(wavenumbers))][intersect(peaks_4899$peak_ind, peaks_COL$peak_ind)]
wavenumbers[-c(1, length(wavenumbers))][intersect(peaks_2099$peak_ind, peaks_COL$peak_ind)]

#peaks_4899 <- peaks_4899[peaks_4899$abs > 0, ]
abline(v = peaks_4899$peak_wavenumber, col = "red")


par(mfrow = c(1, 1))
ylims <- c(-0.075, 0.05)
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 900), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(h = 0, lty = 1, col = "black")
lutzke_wavenumbers(cex = 0.8)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[1], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2), col = plot_cols[2], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), col = plot_cols[3], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), col = plot_cols[4], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(940, 0, "4899", cex = 1.5, adj = c(0, 0.5))
text(940, -0.1, "2099", cex = 1.5, adj = c(0, 0.5))
text(940, -0.2, "COL", cex = 1.5, adj = c(0, 0.5))
text(940, -0.3, "9999SR", cex = 1.5, adj = c(0, 0.5))
#legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
text(930, -0.335, "Plot ref: XXXX")
box()

diffs_wns <- 
  seek_seq(
    which(
      diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) < 0 & (
        !(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2) < 0) &
          !(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2) < 0)
      )
    )
  )


diffs_wns <- 
  seek_seq(
    which(
      (diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) > 0) & 
        (diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2) > 0) &
        (diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2) > 0)
    )
  )

diffs_wns$from <- wavenumbers[diffs_wns$from]
diffs_wns$to <- wavenumbers[diffs_wns$to]

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6))
ylims <- c(-0.04, 0.035)
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 950), ylim = -ylims,
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(h = 0, lty = 1, col = "black")
lutzke_wavenumbers(cex = 0.8)
rect(diffs_wns$from, 1, diffs_wns$to, -1, col = "#00000055", border = NA)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[1], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2), col = plot_cols[2], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), col = plot_cols[3], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), col = plot_cols[4], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(-ylims) + ((max(-ylims) - min(-ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
text(1790, -0.03, "Plot ref: T5Y0")
box()

spec_4899 <- colMeans(spectra[which(metadata$sample == "4899"), ])
spec_2099 <- colMeans(spectra[which(metadata$sample == "2099"), ])
spec_COL <- colMeans(spectra[which(metadata$sample == "COL"), ])
spec_MYB99 <- rbind(spec_4899, spec_2099, spec_COL)
spec_MYB99 <- colMeans(spec_MYB99)
spec_MYB99_plus <- rbind(spec_4899, spec_2099)
spec_MYB99_plus <- colMeans(spec_MYB99_plus)
sd_diff_MYB99 <- sd(diff(colMeans(spectra[which(metadata$sample == "9999SR"), 1:440]), differences = 2) - diff(spec_MYB99[1:440], differences = 2))



par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6))
ylims <- rev(c(-0.014, 0.027))
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 950), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lwd = 2, col = plot_cols[4])
#lutzke_wavenumbers(cex = 0.8, only = c(1605, 1515, 1379, 1168))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(spec_MYB99, differences = 2), col = "black", lwd = 2)
sd_diff_MYB99 <- sd(diff(colMeans(spectra[which(metadata$sample == "9999SR"), 1:440]), differences = 2) - diff(spec_MYB99[1:440], differences = 2))
abline(h = 2 * sd_diff_MYB99, col = "#00000088", lty = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 2, line = 1, at = 0, cex = 1.5, text = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(1800, sd_diff_MYB99 + 0.002, expression(paste("2", sigma)), col = "#00000088", cex = 1.5)
text(1790, -0.01, "Plot ref: EZH4")
text(c(1605, 1515, 1379, 1168), 
     c(0.013, 0.025, 0.01, 0.018), 
     c("v(C=C)", "v(C=C)", expression(paste(delta[sym]*"(CH"[3]*")")), expression(paste(delta, "(CH)"))), cex = 1.5)
legend("bottomright", legend = c("MYB99-", "MYB99/MYB99+ mean"), col = c(plot_cols[4], "black"), lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
box()

(diff_peaks <- find_peaks(diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(spec_MYB99, differences = 2),
                          type = "p", nudge = 0.001, minpeakheight = 2 * sd_diff_MYB99))

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6))
ylims <- rev(c(-0.014, 0.027))
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 950), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lwd = 2, col = plot_cols[4])
#lutzke_wavenumbers(cex = 0.8, only = c(1605, 1515, 1379, 1168))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(spec_MYB99_norm, differences = 2), col = "black", lwd = 2)
sd_diff_MYB99 <- sd(diff(colMeans(spectra[which(metadata$sample == "9999SR"), 1:440]), differences = 2) - diff(spec_MYB99[1:440], differences = 2))
abline(h = 2 * sd_diff_MYB99, col = "#00000088", lty = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 2, line = 1, at = 0, cex = 1.5, text = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(1800, sd_diff_MYB99 + 0.002, expression(paste("2", sigma)), col = "#00000088", cex = 1.5)
text(1790, -0.01, "Plot ref: XXXX")
text(c(1605, 1515, 1379, 1168), 
     c(0.013, 0.025, 0.01, 0.018), 
     c("v(C=C)", "v(C=C)", expression(paste(delta[sym]*"(CH"[3]*")")), expression(paste(delta, "(CH)"))), cex = 1.5)
legend("bottomright", legend = c("MYB99-", "MYB99/MYB99+ mean"), col = c(plot_cols[4], "black"), lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
box()

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6))
ylims <- rev(c(-0.025, 0.03))
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 950), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lwd = 2, col = plot_cols[4])
#lutzke_wavenumbers(cex = 0.8, only = c(1605, 1515, 1379, 1168))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(spec_MYB99_plus, differences = 2), col = "black", lwd = 2)
sd_diff_MYB99 <- sd(diff(colMeans(spectra[which(metadata$sample == "9999SR"), 1:440]), differences = 2) - diff(spec_MYB99[1:440], differences = 2))
abline(h = 2 * sd_diff_MYB99, col = "#00000088", lty = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 2, line = 1, at = 0, cex = 1.5, text = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(1800, sd_diff_MYB99 + 0.002, expression(paste("2", sigma)), col = "#00000088", cex = 1.5)
text(1790, -0.01, "Plot ref: XXXX")
text(c(1605, 1515, 1379, 1168), 
     c(0.013, 0.025, 0.01, 0.018), 
     c("v(C=C)", "v(C=C)", expression(paste(delta[sym]*"(CH"[3]*")")), expression(paste(delta, "(CH)"))), cex = 1.5)
legend("bottomright", legend = c("MYB99-", "MYB99/MYB99+ mean"), col = c(plot_cols[4], "black"), lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
box()


which(diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(spec_MYB99, differences = 2) > (2 * sd_diff_MYB99))

diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(spec_MYB99, differences = 2)

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6))
#ylims <- rev(c(-0.014, 0.027))
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(spec_MYB99_plus, differences = 2), 
     xlim = c(1800, 950), #ylim = range,
     type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lwd = 2, col = plot_cols[3])
#lutzke_wavenumbers(cex = 0.8, only = c(1605, 1515, 1379, 1168))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(spec_MYB99_plus, differences = 2), col = "black", lwd = 2)
sd_diff_MYB99 <- sd(diff(colMeans(spectra[which(metadata$sample == "9999SR"), 1:440]), differences = 2) - diff(spec_MYB99[1:440], differences = 2))
abline(h = 2 * sd(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(spec_MYB99_plus, differences = 2)[1:440]), col = "#00000088", lty = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(1800, sd_diff_MYB99 + 0.002, expression(paste("2", sigma)), col = "#00000088", cex = 1.5)
#text(1790, -0.01, "Plot ref: EZH4")
text(c(1605, 1515, 1379, 1168), 
     c(0.013, 0.025, 0.01, 0.018), 
     c("v(C=C)", "v(C=C)", expression(paste(delta[sym]*"(CH"[3]*")")), expression(paste(delta, "(CH)"))), cex = 1.5)
legend("bottomright", legend = c("MYB99-", "MYB99/MYB99+ mean"), col = c(plot_cols[4], "black"), lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
box()

which((diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(spec_MYB99_plus, differences = 2)[1:440]) > 2 * sd(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(spec_MYB99_plus, differences = 2)[1:440]))

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6))
#ylims <- rev(c(-0.014, 0.027))
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 950), #ylim = range,
     type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lwd = 2, col = plot_cols[3])
#lutzke_wavenumbers(cex = 0.8, only = c(1605, 1515, 1379, 1168))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), col = "black", lwd = 2)
sd_diff_COL_9999SR <- -sd(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2)[1:440])
abline(h = 2 * sd_diff_COL_9999SR, col = "#00000088", lty = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(1800, sd_diff_COL_9999SR + 0.002, expression(paste("2", sigma)), col = "#00000088", cex = 1.5)
#text(1790, -0.01, "Plot ref: EZH4")
text(c(1605, 1515, 1379, 1168), 
     c(0.013, 0.025, 0.01, 0.018), 
     c("v(C=C)", "v(C=C)", expression(paste(delta[sym]*"(CH"[3]*")")), expression(paste(delta, "(CH)"))), cex = 1.5)
#legend("bottomright", legend = c("MYB99-", "MYB99/MYB99+ mean"), col = c(plot_cols[4], "black"), lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
box()

which((diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2)[1:440]) < 2 * sd_diff_COL_9999SR)








rect(diffs_wns$from, 1, diffs_wns$to, -1, col = "#00000033", border = NA)



plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), xlim = c(1800, 948), 
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
abline(h = 0, lty = 1, lwd = 2, col = "black")
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[1], lwd = 2)
lutzke_wavenumbers()
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, at = 0, label = "0", las = 1)
box()


spec_diff <- colMeans(spectra[which(metadata$sample == "2099"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
plot(wavenumbers, spec_diff, xlim = c(1800, 948), type = "n", axes = F, xlab = "", ylab = "")
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers, spec_diff, lwd = 2, col = "lightgrey")
points(wavenumbers[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
abline(h = 0, lty = 2, col = "black")
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, at = 0, label = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Mean difference")
box()
rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
rect(900, 0, 940, -10, col = plot_cols[3], border = NA)
text(1050, max(spec_diff), "Plot ref: OQZU")


plot(wavenumbers, spectra[1, ], xlim = rev(range(wavenumbers)), type = "n", ylim = range(spectra), axes = F, xlab = "", ylab = "")
for(i in 1:9){
  lines(wavenumbers, spectra[i, ], col = plot_cols[1])  
}
for(i in 10:17){
  lines(wavenumbers, spectra[i, ], col = plot_cols[2])  
}
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 2, line = 2, at = min(spectra) + ((max(spectra) - min(spectra)) / 2), cex = 1.5, text = "Relative absorbance")
box()

plot(wavenumbers, spec_diff, type = "n", xlim = rev(range(wavenumbers)), axes = F, xlab = "", ylab = "")
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers, spec_diff, lwd = 2, col = "lightgrey")
points(wavenumbers[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
abline(h = 0, lty = 2, col = "black")
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Mean difference")
box()
rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
rect(900, 0, 940, -10, col = plot_cols[2], border = NA)

metadata$treatment <- NA
metadata$treatment[metadata$sample == "9999SR"] <- "MYB99-"
metadata$treatment[metadata$sample == "COL"] <- "CONTROL"
metadata$treatment[is.na(metadata$treatment)] <- "MYB99+"




############################ Principal component analysis



apca <- prcomp(spectra[, 1:440], center = T, scale. = F)

plot_cex <- 2
pc1 <- 1
pc2 <- 2
par(mfrow = c(1, 2), mar = c(5.6, 4.6, 1.1, 1.1), pty = "s")
plot(apca$x[, pc1], apca$x[, pc2], type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- apca$x[, 1:2] %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids$pch <- ifelse(centroids$sample %in% c("991W", "995W", "MLP8W"), 4, 3)
points(PC2 ~ PC1, centroids, pch = centroids$pch, cex = 4, col = rev(plot_cols), lwd = 2)
points(apca$x[metadata$sample == "MLP8W", pc1], apca$x[metadata$sample == "MLP8W", pc2], col = plot_cols[1], pch = 4, cex = plot_cex, lwd = 2)
points(apca$x[metadata$sample == "MLP8", pc1], apca$x[metadata$sample == "MLP8", pc2], col = plot_cols[2], pch = 18, cex = plot_cex)
points(apca$x[metadata$sample == "9999SR", pc1], apca$x[metadata$sample == "9999SR", pc2], col = plot_cols[3], pch = 17, cex = plot_cex)
points(apca$x[metadata$sample == "995W", pc1], apca$x[metadata$sample == "995W", pc2], col = plot_cols[4], pch = 4, cex = plot_cex, lwd = 2)
points(apca$x[metadata$sample == "995", pc1], apca$x[metadata$sample == "995", pc2], col = plot_cols[5], pch = 18, cex = plot_cex)
points(apca$x[metadata$sample == "991W", pc1], apca$x[metadata$sample == "991W", pc2], col = plot_cols[6], pch = 4, cex = plot_cex, lwd = 2)
points(apca$x[metadata$sample == "991", pc1], apca$x[metadata$sample == "991", pc2], col = plot_cols[7], pch = 18, cex = plot_cex)
points(apca$x[metadata$sample == "COL", pc1], apca$x[metadata$sample == "COL", pc2], col = plot_cols[8], pch = 16, cex = plot_cex)
points(apca$x[metadata$sample == "4899", pc1], apca$x[metadata$sample == "4899", pc2], col = plot_cols[9], pch = 15, cex = plot_cex)
points(apca$x[metadata$sample == "2099", pc1], apca$x[metadata$sample == "2099", pc2], col = plot_cols[10], pch = 15, cex = plot_cex)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = min(apca$x[, pc1]) + ((max(apca$x[, pc1]) - min(apca$x[, pc1])) / 2), cex = 1.5, text = "PC1 (62.4%)")
mtext(side = 2, line = 3, at = min(apca$x[, pc2]) + ((max(apca$x[, pc2]) - min(apca$x[, pc2])) / 2), cex = 1.5, text = "PC2 (14.5%)")
box()
legend("bottomright", legend = factor_levels, pch = plot_pch, col = plot_cols, cex = 1.5)
#text(apca$x[, pc1], apca$x[, pc2], metadata$sample)
#text(3.2, 1.5, "Plot ref: YNMP")

which(apca$x[, 1] < -4)

pc1 <- 2
pc2 <- 3
plot(apca$x[, pc1], apca$x[, pc2], type = "n", axes = F, xlab = "", ylab = "", xlim = c(min(apca$x[, pc1]), 2.2))
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- apca$x[, c(pc1, pc2)] %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids$pch <- ifelse(centroids$sample %in% c("991W", "995W", "MLP8W"), 4, 3)
points(centroids[, 2], centroids[, 3], pch = centroids$pch, cex = 3, col = rev(plot_cols), lwd = 2)
points(apca$x[metadata$sample == "MLP8W", pc1], apca$x[metadata$sample == "MLP8W", pc2], col = plot_cols[1], pch = 4, cex = plot_cex, lwd = 2)
points(apca$x[metadata$sample == "MLP8", pc1], apca$x[metadata$sample == "MLP8", pc2], col = plot_cols[2], pch = 18, cex = plot_cex)
points(apca$x[metadata$sample == "9999SR", pc1], apca$x[metadata$sample == "9999SR", pc2], col = plot_cols[3], pch = 17, cex = plot_cex)
points(apca$x[metadata$sample == "995W", pc1], apca$x[metadata$sample == "995W", pc2], col = plot_cols[4], pch = 4, cex = plot_cex, lwd = 2)
points(apca$x[metadata$sample == "995", pc1], apca$x[metadata$sample == "995", pc2], col = plot_cols[5], pch = 18, cex = plot_cex)
points(apca$x[metadata$sample == "991W", pc1], apca$x[metadata$sample == "991W", pc2], col = plot_cols[6], pch = 4, cex = plot_cex, lwd = 2)
points(apca$x[metadata$sample == "991", pc1], apca$x[metadata$sample == "991", pc2], col = plot_cols[7], pch = 18, cex = plot_cex)
points(apca$x[metadata$sample == "COL", pc1], apca$x[metadata$sample == "COL", pc2], col = plot_cols[8], pch = 16, cex = plot_cex)
points(apca$x[metadata$sample == "4899", pc1], apca$x[metadata$sample == "4899", pc2], col = plot_cols[9], pch = 15, cex = plot_cex)
points(apca$x[metadata$sample == "2099", pc1], apca$x[metadata$sample == "2099", pc2], col = plot_cols[10], pch = 15, cex = plot_cex)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = min(apca$x[, pc1]) + ((max(apca$x[, pc1]) - min(apca$x[, pc1])) / 2), cex = 1.5, text = "PC2 (14.5%)")
mtext(side = 2, line = 3, at = min(apca$x[, pc2]) + ((max(apca$x[, pc2]) - min(apca$x[, pc2])) / 2), cex = 1.5, text = "PC3 (9.1%)")
box()
#legend("bottomright", legend = factor_levels, pch = c(18, 17, 18, 16, 18, 15, 15), col = plot_cols, bty = "n", cex = 1.5)
#text(apca$x[, pc1], apca$x[, pc2], metadata$sample)
#text(3.2, 1.5, "Plot ref: YNMP")

plot(wavenumbers[1:440], apca$rotation[, 1], xlim = rev(range(wavenumbers[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[1:440], apca$rotation[, 1], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
#mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca$rotation[, pc1]) + ((max(apca$rotation[, pc1]) - min(apca$rotation[, pc1])) / 2), cex = 1.5, text = "PC1 loadings")
box()
lutzke_wavenumbers()

#par(mfrow = c(1, 1), mar = c(5.1, 6.1, 4.1, 0.1))
plot(wavenumbers[1:440], apca$rotation[, 2], xlim = rev(range(wavenumbers[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[1:440], apca$rotation[, 2], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca$rotation[, 2]) + ((max(apca$rotation[, 2]) - min(apca$rotation[, 2])) / 2), cex = 1.5, text = "PC2 loadings")
box()
lutzke_wavenumbers()
#mtext(side = 1, line = 3, at = 960, text = "Plot ref: 5512")

plot(wavenumbers[1:440], apca$rotation[, 3], xlim = rev(range(wavenumbers[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[1:440], apca$rotation[, 3], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca$rotation[, 3]) + ((max(apca$rotation[, 2]) - min(apca$rotation[, 3])) / 2), cex = 1.5, text = "PC2 loadings")
box()
lutzke_wavenumbers()


apca_crop <- prcomp(spectra[!metadata$sample %in% c("991", "995", "995W", "MLP8"), 1:440], center = T, scale. = F)
metadata_crop <- metadata[!metadata$sample %in% c("991", "995", "995W", "MLP8"), ]

plot_cex <- 2
pc1 <- 1
pc2 <- 2
par(mfrow = c(1, 2), mar = c(5.6, 4.6, 1.1, 1.1), pty = "s")
plot(apca_crop$x[, pc1], apca_crop$x[, pc2], type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- apca_crop$x[, 1:2] %>% bind_cols(sample = metadata_crop$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
#centroids$pch <- ifelse(centroids$sample %in% c("991W", "995W", "MLP8W"), 4, 3)
points(PC2 ~ PC1, centroids, pch = 4, cex = 4, col = rev(plot_cols[c(1, 3, 7:10)]), lwd = 3)
points(apca_crop$x[metadata_crop$sample == "MLP8W", pc1], apca_crop$x[metadata_crop$sample == "MLP8W", pc2], col = plot_cols[1], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "MLP8", pc1], apca_crop$x[metadata_crop$sample == "MLP8", pc2], col = plot_cols[2], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "9999SR", pc1], apca_crop$x[metadata_crop$sample == "9999SR", pc2], col = plot_cols[3], pch = 17, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "995W", pc1], apca_crop$x[metadata_crop$sample == "995W", pc2], col = plot_cols[4], pch = 4, cex = plot_cex, lwd = 2)
points(apca_crop$x[metadata_crop$sample == "995", pc1], apca_crop$x[metadata_crop$sample == "995", pc2], col = plot_cols[5], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "991W", pc1], apca_crop$x[metadata_crop$sample == "991W", pc2], col = plot_cols[6], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "991", pc1], apca_crop$x[metadata_crop$sample == "991", pc2], col = plot_cols[7], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "COL", pc1], apca_crop$x[metadata_crop$sample == "COL", pc2], col = plot_cols[8], pch = 16, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "4899", pc1], apca_crop$x[metadata_crop$sample == "4899", pc2], col = plot_cols[9], pch = 15, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "2099", pc1], apca_crop$x[metadata_crop$sample == "2099", pc2], col = plot_cols[10], pch = 15, cex = plot_cex)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = min(apca_crop$x[, pc1]) + ((max(apca_crop$x[, pc1]) - min(apca_crop$x[, pc1])) / 2), cex = 1.5, text = "PC1 (52.1%)")
mtext(side = 2, line = 3, at = min(apca_crop$x[, pc2]) + ((max(apca_crop$x[, pc2]) - min(apca_crop$x[, pc2])) / 2), cex = 1.5, text = "PC2 (25.5%)")
box()
legend("topright", legend = c("MLP8", "9999SR", "991", "COL", "4899", "2099"), pch = c(18, 17, 18, 16, 15, 15), col = plot_cols[c(1, 3, 7:10)], cex = 1.8)
#text(apca_crop$x[, pc1], apca_crop$x[, pc2], metadata_crop$sample)
#text(3.2, 1.5, "Plot ref: YNMP")

which(apca_crop$x[, 1] < -4)

pc1 <- 2
pc2 <- 3
plot(apca_crop$x[, pc1], apca_crop$x[, pc2], type = "n", axes = F, xlab = "", ylab = "", xlim = c(min(apca_crop$x[, pc1]), 2.2))
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- apca_crop$x[, c(pc1, pc2)] %>% bind_cols(sample = metadata_crop$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids$pch <- ifelse(centroids$sample %in% c("991W", "995W", "MLP8W"), 4, 3)
points(centroids[, 2], centroids[, 3], pch = 4, cex = 4, col = rev(plot_cols[c(1, 3, 7:10)]), lwd = 3)
points(apca_crop$x[metadata_crop$sample == "MLP8W", pc1], apca_crop$x[metadata_crop$sample == "MLP8W", pc2], col = plot_cols[1], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "MLP8", pc1], apca_crop$x[metadata_crop$sample == "MLP8", pc2], col = plot_cols[2], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "9999SR", pc1], apca_crop$x[metadata_crop$sample == "9999SR", pc2], col = plot_cols[3], pch = 17, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "995W", pc1], apca_crop$x[metadata_crop$sample == "995W", pc2], col = plot_cols[4], pch = 4, cex = plot_cex, lwd = 2)
points(apca_crop$x[metadata_crop$sample == "995", pc1], apca_crop$x[metadata_crop$sample == "995", pc2], col = plot_cols[5], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "991W", pc1], apca_crop$x[metadata_crop$sample == "991W", pc2], col = plot_cols[6], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "991", pc1], apca_crop$x[metadata_crop$sample == "991", pc2], col = plot_cols[7], pch = 18, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "COL", pc1], apca_crop$x[metadata_crop$sample == "COL", pc2], col = plot_cols[8], pch = 16, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "4899", pc1], apca_crop$x[metadata_crop$sample == "4899", pc2], col = plot_cols[9], pch = 15, cex = plot_cex)
points(apca_crop$x[metadata_crop$sample == "2099", pc1], apca_crop$x[metadata_crop$sample == "2099", pc2], col = plot_cols[10], pch = 15, cex = plot_cex)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = min(apca_crop$x[, pc1]) + ((max(apca_crop$x[, pc1]) - min(apca_crop$x[, pc1])) / 2), cex = 1.5, text = "PC2 (25.5%)")
mtext(side = 2, line = 3, at = min(apca_crop$x[, pc2]) + ((max(apca_crop$x[, pc2]) - min(apca_crop$x[, pc2])) / 2), cex = 1.5, text = "PC3 (11.7%)")
box()
#legend("bottomright", legend = factor_levels, pch = c(18, 17, 18, 16, 18, 15, 15), col = plot_cols, bty = "n", cex = 1.5)
#text(apca_crop$x[, pc1], apca_crop$x[, pc2], metadata_crop$sample)
#text(3.2, 1.5, "Plot ref: YNMP")
par(mfrow = c(1, 1), mar = c(6.1, 6.6, 2.1, 1.1), pty = "m")
plot(wavenumbers[1:440], apca_crop$rotation[, 1], xlim = rev(range(wavenumbers[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[1:440], apca_crop$rotation[, 1], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca_crop$rotation[, 1]) + ((max(apca_crop$rotation[, 1]) - min(apca_crop$rotation[, 1])) / 2), cex = 1.5, text = "PC1 loadings")
box()
lutzke_wavenumbers()

#par(mfrow = c(1, 1), mar = c(5.1, 6.1, 4.1, 0.1))
plot(wavenumbers[1:440], apca_crop$rotation[, 2], xlim = rev(range(wavenumbers[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[1:440], apca_crop$rotation[, 2], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca_crop$rotation[, 2]) + ((max(apca_crop$rotation[, 2]) - min(apca_crop$rotation[, 2])) / 2), cex = 1.5, text = "PC2 loadings")
box()
lutzke_wavenumbers()
#mtext(side = 1, line = 3, at = 960, text = "Plot ref: 5512")

plot(wavenumbers[1:440], apca_crop$rotation[, 3], xlim = rev(range(wavenumbers[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[1:440], apca_crop$rotation[, 3], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca_crop$rotation[, 3]) + ((max(apca_crop$rotation[, 2]) - min(apca_crop$rotation[, 3])) / 2), cex = 1.5, text = "PC2 loadings")
box()
lutzke_wavenumbers()

spectra[, nearest(names(spectra), 1605)$index] /
  spectra[, nearest(names(spectra), 1515)$index]



################ derivative

spectra_d <- t(apply(spectra, 1, diff, differences = 2)) %>% as.data.frame()
wavenumbers_d <- wavenumbers[2:(length(wavenumbers) - 1)]


# par(mfrow = c(2, 1), mar = c(5.1, 4.1, 0.1, 0.1))
# plot(wavenumbers_d, spectra_d[1, ], xlim = c(1800, 948), type = "n", ylim = range(spectra_d), axes = F, xlab = "", ylab = "")
input.var.equal <- F
t_test_results <- numeric(ncol(spectra_d))
for(i in 1:ncol(spectra_d)){
  t_test <- t.test(spectra_d[1:9, i], spectra_d[10:17, i], var.equal = input.var.equal)
  t_test_results[i] <- t_test$p.value
}

t_test_results[which(t_test_results > 0.05)] <- NA

seek_table <- seek_seq(which(t_test_results <= 0.05))
seek_table$from <- wavenumbers_d[seek_table$from]
seek_table$to <- wavenumbers_d[seek_table$to]
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)

# for(i in 1:9){
#   lines(wavenumbers_d, spectra_d[i, ], col = plot_cols[1])  
# }
# for(i in 10:17){
#   lines(wavenumbers_d, spectra_d[i, ], col = plot_cols[2])  
# }
# axis(1, cex.axis = 1.5, tck = 0.01)
# 
# mtext(side = 2, line = 2, at = min(spectra_d) + ((max(spectra_d) - min(spectra_d)) / 2), cex = 1.5, text = "Relative absorbance")
# #internal_axis_ticks(1, 1)
# legend("bottomright", legend = unique(metadata$sample), lty = 1, col = plot_cols, bty = "n", lwd = 2, cex = 1.5)
# box()

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 2.1, 0.1))
spec_diff <- colMeans(spectra_d[1:9, ]) - colMeans(spectra_d[10:17, ])
plot(wavenumbers_d, spec_diff, xlim = c(1800, 948), type = "n", axes = F, xlab = "", ylab = "")
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers_d, spec_diff, lwd = 2, col = "lightgrey")
points(wavenumbers_d[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
abline(h = 0, lty = 2, col = "black")
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, at = 0, label = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Mean difference")
box()
rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
rect(900, 0, 940, -10, col = plot_cols[2], border = NA)
phil_wavenumbers()
text(1770, max(spec_diff), "Plot ref: EL7Z")



# plot(wavenumbers_d, spectra_d[1, ], xlim = rev(range(wavenumbers_d)), type = "n", ylim = range(spectra_d), axes = F, xlab = "", ylab = "")
# for(i in 1:9){
#   lines(wavenumbers_d, spectra_d[i, ], col = plot_cols[1])  
# }
# for(i in 10:17){
#   lines(wavenumbers_d, spectra_d[i, ], col = plot_cols[2])  
# }
# rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
# axis(1, cex.axis = 1.5, tck = 0.01)
# mtext(side = 2, line = 2, at = min(spectra_d) + ((max(spectra_d) - min(spectra_d)) / 2), cex = 1.5, text = "Relative absorbance")
# box()
# 
# plot(wavenumbers_d, spec_diff, type = "n", xlim = rev(range(wavenumbers_d)), axes = F, xlab = "", ylab = "")
# rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
# lines(wavenumbers_d, spec_diff, lwd = 2, col = "lightgrey")
# points(wavenumbers_d[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
# abline(h = 0, lty = 2, col = "black")
# axis(1, cex.axis = 1.5, tck = 0.01)
# mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
# mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Mean difference")
# box()
# rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
# rect(900, 0, 940, -10, col = plot_cols[2], border = NA)



pc1 <- 2
pc2 <- 3
apca <- prcomp(spectra_d[, 1:440], center = T, scale. = T)
par(mfrow = c(1, 1))
plot(apca$x[, pc1], apca$x[, pc2], type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
points(apca$x[metadata$sample == "9999SR", pc1], apca$x[metadata$sample == "9999SR", pc2], col = plot_cols[1], pch = 16)
points(apca$x[metadata$sample == "COL", pc1], apca$x[metadata$sample == "COL", pc2], col = plot_cols[2], pch = 16)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = min(apca$x[, pc1]) + ((max(apca$x[, pc1]) - min(apca$x[, pc1])) / 2), cex = 1.5, text = "PC1 (35.0%)")
mtext(side = 2, line = 3, at = min(apca$x[, pc2]) + ((max(apca$x[, pc2]) - min(apca$x[, pc2])) / 2), cex = 1.5, text = "PC2 (20.0%)")
box()
legend("topright", legend = unique(metadata$sample), pch = 16, col = plot_cols, bty = "n", cex = 1.5)
text(10, -8, "Plot ref: 3LHR")

par(mfrow = c(1, 1), mar = c(5.1, 6.1, 2.1, 0.1))
plot(wavenumbers_d[1:440], apca$rotation[, pc2], xlim = rev(range(wavenumbers_d[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 2, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 2, col = "#FF000055", lwd = 2)
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers_d[1:440], apca$rotation[, pc2], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca$rotation[, pc2]) + ((max(apca$rotation[, pc2]) - min(apca$rotation[, pc2])) / 2), cex = 1.5, text = "PC2 loadings")
box()
phil_wavenumbers()
text(1770, 0.13, "Plot ref: 43XY")


(detection_sensitivity_compare <- 3)
peak_table_list <- vector("list", length = nrow(spectra))

for(i in 1:nrow(spectra_d)){
  print(i)
  (my_spectrum1 <- spectra_d[i, ])
  
  (xx <- names(my_spectrum1))
  #my_spectrum1 <- as.numeric(my_spectrum1)
  
  (x <- which(my_spectrum1 < detection_sensitivity_compare * -sd(my_spectrum1))) #### compare against one or another
  
  seek_seq <- function(x){
    skipped_indices <- which(diff(x) != 1) #### only works if [1] is not a skip
    to <- x[skipped_indices]
    to <- c(to, x[length(x)])
    from <- c(1, skipped_indices + 1)
    from <- x[from]
    data.frame(from = from, to = to)
  }
  
  (seqs <- seek_seq(x))
  (peak_sequences <- apply(seqs, 1, function(x){x[1]:x[2]}))
  (peak_wns <- lapply(peak_sequences, function(x){
    names(which.min(my_spectrum1[x]))
  }))
  (peak_values <- lapply(peak_sequences, function(x){
    min(my_spectrum1[x])
  }))
  (peaks <- data.frame(wn = unlist(peak_wns), value = -unlist(peak_values)))
  (peaks$ID<- row.names(my_spectrum1))
  (peak_table_list[[i]] <- peaks)
}
peaks_df <- as.data.frame(do.call(rbind, peak_table_list))

peaks_df$value<- as.numeric(peaks_df$value)
#peaks_df$peaks_id<- peaks_df$ID %>% str_replace_all("[0-9]", "") %>% str_replace_all("_", "")
#peaks_df$peaks_id <- factor(peaks_df$peaks_id)
peaks_df$wn <- peaks_df$wn %>% strip_non_numbers()

peaks1 <- peaks_df[as.numeric(peaks_df$ID) <= 9, ]
(nd <- data.frame(table(peaks1$wn)))
(nd$Var1 <- nd$Var1 %>% as.character() %>% as.numeric())
(nd$Freq_per<- (nd$Freq* 100)/ 9)

(peaks2 <- peaks_df[as.numeric(peaks_df$ID) > 9, ])
(nd2 <- data.frame(table(peaks2$wn)))
(nd2$Var1 <- nd2$Var1 %>% as.character() %>% as.numeric())
(nd2$Freq_per <- (nd2$Freq* 100)/ 8)

#nd$quality <- "good"

plot(nd$Var1, nd$Freq_per, type = "n", xlim = rev(range(wavenumbers)), axes = F, xlab = expression(paste("Wavenumber (cm"^-1*")")), ylab = "% of spectra", cex.lab = 1.5, ylim = c(0, 100))
arrows(nd$Var1, 0, nd$Var1, nd$Freq_per, angle = 0, length = 0, col = plot_cols[1], lwd = 3)
arrows(nd2$Var1, 0, nd2$Var1, nd2$Freq_per, angle = 0, length = 0, col = plot_cols[2], lwd = 3)
axis(1, tck = 0.01, cex.axis = 1.3)
axis(2, las = 1, tck = 0.01, cex.axis = 1.3)
box()




#############################################################################################################


par(mfrow = c(1,2))                            # Splitting the plot in 2.
#biplot(apca)                                    # In-built biplot() R func.
PCA <- apca

choices = 1:2                                  # Selecting first two PC's
scale = 1                                      # Default
scores = PCA$x                                  # The scores
lam = PCA$sdev[choices]                        # Sqrt e-vals (lambda) 2 PC's
n = nrow(scores)                               # no. rows scores
lam = lam * sqrt(n)                            # See below.


x = t(t(scores[,choices]) / lam)         # scaled scores
y = t(t(PCA$rotation[,choices]) * lam)         # scaled eigenvecs (loadings)

n = nrow(x)                                    # Same as dataset (150)
p = nrow(y)                                    # Three var -> 3 rows

xlabs = 1L:n
xlabs = as.character(xlabs)                    # no. from 1 to 150 
dimnames(x) = list(xlabs, dimnames(x)[[2L]])   # no's and PC1 / PC2

ylabs = dimnames(y)[[1L]]                      # Iris species
ylabs = as.character(ylabs)
dimnames(y) <- list(ylabs, dimnames(y)[[2L]])  # Species and PC1/PC2


unsigned.range = function(x) c(-abs(min(x, na.rm = TRUE)), 
                               abs(max(x, na.rm = TRUE)))
rangx1 = unsigned.range(x[, 1L])               # Range first col x
rangx2 = unsigned.range(x[, 2L])               # Range second col x
rangy1 = unsigned.range(y[, 1L])               # Range 1st scaled evec
rangy2 = unsigned.range(y[, 2L])               # Range 2nd scaled evec


(xlim = ylim = rangx1 = rangx2 = range(rangx1, rangx2))

(ratio = max(rangy1/rangx1, rangy2/rangx2)) 

par(mfrow = c(1, 2), pty = "s")                                
pc1 <- 1
pc2 <- 2
plot(x, type = "n", 
     #xlim = xlim, 
     #ylim = ylim, 
     xlim = c(-0.32, 0.32), 
     ylim = c(-0.32, 0.32), 
     axes = F, 
     xlab = "", 
     ylab = ""
)  
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- x %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids <- centroids[c(2, 1, 4, 3), ]
points(PC2 ~ PC1, centroids, pch = 3, cex = 3, col = plot_cols, lwd = 2)
points(x[metadata$sample == "4899", pc1], x[metadata$sample == "4899", pc2], col = plot_cols[1], pch = 17)
points(x[metadata$sample == "2099", pc1], x[metadata$sample == "2099", pc2], col = plot_cols[2], pch = 17)
points(x[metadata$sample == "COL", pc1], x[metadata$sample == "COL", pc2], col = plot_cols[3], pch = 16)
points(x[metadata$sample == "9999SR", pc1], x[metadata$sample == "9999SR", pc2], col = plot_cols[4], pch = 15)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 0, cex = 1.5, text = "PC1 scores (59.0%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 scores (26.6%)")
legend("topleft", legend = c("4899", "2099", "COL", "9999SR"), pch = c(17, 17, 16, 15), col = plot_cols, bty = "n", cex = 1.5)
box()

(new_xlim = xlim * ratio)

(new_ylim = ylim * ratio)

lutzke_indices <- which(wavenumbers %in% lutzke_wavenumbers_only[lutzke_wavenumbers_only < 1800])

axis_lim <- 1.15
plot(y, axes = FALSE, type = "n", 
     xlim = new_xlim, 
     ylim = new_ylim, 
     xlab = "", 
     ylab = "", 
     #xlim = c(-axis_lim, axis_lim),
     #ylim = c(-axis_lim, axis_lim)
)

abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 0, cex = 1.5, text = "PC1 loadings (59.0%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 loadings (26.6%)")

text(y[lutzke_indices, ], labels = round(strip_non_numbers(ylabs[lutzke_indices])), col = 1, cex = 0.8)  # This just prints the species
arrow.len <- 0.1
arrows(0, 0, y[lutzke_indices, 1L] * 0.8, y[lutzke_indices, 2L] * 0.8, length = arrow.len, col = "#00000033")

non_sporopollenin_indices <- sapply(c(1642, 1022, 1541, 1558, 1455, 1750, 1772, 1080), function(x){which.min(abs(x - wavenumbers))})
text(y[non_sporopollenin_indices, ], labels = round(strip_non_numbers(ylabs[non_sporopollenin_indices])), cex = 0.8, col = "red")  # This just prints the species
arrows(0, 0, y[non_sporopollenin_indices, 1L] * 0.8, y[non_sporopollenin_indices, 2L] * 0.8, length = arrow.len, col = "#FF000033")
legend("topleft", lty = 1, col = 1:2, legend = c("Sporopollenin", "Other (amide, etc.)"), cex = 1.5, bty = "n")
text(1.2, 1.5, "Plot ref: QCIV")
box()

####################################################

par(mfrow = c(1, 2), mar = c(5.1, 5.1, 1.1, 1.1), pty = "s")
cex <- 1.5
pc1 <- 1
pc2 <- 2
plot(x, type = "n", 
     #xlim = xlim, 
     #ylim = ylim, 
     xlim = c(-0.32, 0.32), 
     ylim = c(-0.32, 0.32), 
     axes = F, 
     xlab = "", 
     ylab = ""
)  
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- x %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids <- centroids[c(2, 1, 4, 3), ]
points(PC2 ~ PC1, centroids, pch = 3, cex = 3, col = plot_cols, lwd = 2)
points(x[metadata$sample == "4899", pc1], x[metadata$sample == "4899", pc2], col = plot_cols[2], pch = 17, cex = cex)
points(x[metadata$sample == "2099", pc1], x[metadata$sample == "2099", pc2], col = plot_cols[1], pch = 17, cex = cex)
points(x[metadata$sample == "COL", pc1], x[metadata$sample == "COL", pc2], col = plot_cols[3], pch = 16, cex = cex)
points(x[metadata$sample == "9999SR", pc1], x[metadata$sample == "9999SR", pc2], col = plot_cols[4], pch = 15, cex = cex)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 0, cex = 1.5, text = "PC1 scores (59.0%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 scores (26.6%)")
legend("topleft", legend = c("2099", "4899", "COL", "9999SR"), pch = c(17, 17, 16, 15), col = plot_cols, bty = "n", cex = 1.5)
box()

par(mfrow = c(1, 1), pty = "m")
plot(wavenumbers[1:440], apca$rotation[, 2], xlim = rev(range(wavenumbers[1:440])), ylim = range(apca$rotation[, 2]) * 1.05, lwd = 2, 
     #xlab = lab("wn"), 
     #ylab = "PC1 loadings", 
     xlab = "", ylab = "",
     axes = F, type = "n")
abline(h = 0)
# abline(v = 1168)
# abline(v = 1515)
# abline(v = 1605)
# abline(v = 1025)
abline(h = sqrt(1/440), col = "#00000055", lwd = 2, lty = 2)
abline(h = -sqrt(1/440), col = "#00000055", lwd = 2, lty = 2)
lines(wavenumbers[1:440], apca$rotation[, 2], lwd = 2, col = "#00000099")
loadings2_threshold <- which(apca$rotation[, 2] > sqrt(1/440) |
                               apca$rotation[, 2] < -sqrt(1/440))
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = midpoint(wavenumbers[1:440]), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = midpoint(apca$rotation[, 2]), cex = 1.5, text = "PC2 loadings")
loadings2_peaks <- find_peaks(apca$rotation[, 2], minpeakheight = sqrt(1/440), nudge = 0.008)

text(loadings2_peaks$peak_wavenumber, loadings2_peaks$nudge, rev(letters[1:length(loadings2_peaks$peak_wavenumber)]), cex = 1.2)
#text(1050, -0.09, "Plot ref: VED7")
box()

plot(wavenumbers[1:440], apca$rotation[, 1], xlim = rev(range(wavenumbers[1:440])), ylim = range(apca$rotation[, 1]) * 1.05, lwd = 2, 
     #xlab = lab("wn"), 
     #ylab = "PC1 loadings", 
     xlab = "", ylab = "",
     axes = F, type = "n")
abline(h = 0)
# abline(v = 1168)
# abline(v = 1515)
# abline(v = 1605)
# abline(v = 1025)
abline(h = sqrt(1/440), col = "#00000055", lwd = 2, lty = 2)
abline(h = -sqrt(1/440), col = "#00000055", lwd = 2, lty = 2)
lines(wavenumbers[1:440], apca$rotation[, 1], lwd = 2, col = "#00000099")
loadings1_threshold <- which(apca$rotation[, 1] > sqrt(1/440) |
                               apca$rotation[, 1] < -sqrt(1/440))

axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = midpoint(wavenumbers[1:440]), cex = 1.5, text = lab("wn"))
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = midpoint(apca$rotation[, 1]), cex = 1.5, text = "PC1 loadings")
loadings1_peaks <- find_peaks(apca$rotation[, 1], minpeakheight = sqrt(1/440), nudge = 0.008)
text(loadings1_peaks$peak_wavenumber, loadings1_peaks$nudge, rev(letters[1:length(loadings1_peaks$peak_wavenumber)]), cex = 1.2)
text(1050, -0.09, "Plot ref: VED7")
box()


find_peaks(apca$rotation[, 2], minpeakheight = sqrt(1/440), nudge = 0.008)
####################################################

par(mfrow = c(1,2))                            # Splitting the plot in 2.
#biplot(apca)                                    # In-built biplot() R func.
PCA <- apca

choices = 1:2                                  # Selecting first two PC's
scale = 1                                      # Default
scores = PCA$x                                  # The scores
lam = PCA$sdev[choices]                        # Sqrt e-vals (lambda) 2 PC's
n = nrow(scores)                               # no. rows scores
lam = lam * sqrt(n)                            # See below.


x = t(t(scores[,choices]) / lam)         # scaled scores
y = t(t(PCA$rotation[,choices]) * lam)         # scaled eigenvecs (loadings)

n = nrow(x)                                    # Same as dataset (150)
p = nrow(y)                                    # Three var -> 3 rows

xlabs = 1L:n
xlabs = as.character(xlabs)                    # no. from 1 to 150 
dimnames(x) = list(xlabs, dimnames(x)[[2L]])   # no's and PC1 / PC2

ylabs = dimnames(y)[[1L]]                      # Iris species
ylabs = as.character(ylabs)
dimnames(y) <- list(ylabs, dimnames(y)[[2L]])  # Species and PC1/PC2


unsigned.range = function(x) c(-abs(min(x, na.rm = TRUE)), 
                               abs(max(x, na.rm = TRUE)))
rangx1 = unsigned.range(x[, 1L])               # Range first col x
rangx2 = unsigned.range(x[, 2L])               # Range second col x
rangy1 = unsigned.range(y[, 1L])               # Range 1st scaled evec
rangy2 = unsigned.range(y[, 2L])               # Range 2nd scaled evec


(xlim = ylim = rangx1 = rangx2 = range(rangx1, rangx2))

(ratio = max(rangy1/rangx1, rangy2/rangx2)) 

par(mfrow = c(1, 2), pty = "s")                                
pc1 <- 1
pc2 <- 2
plot(x, type = "n", 
     #xlim = xlim, 
     #ylim = ylim, 
     xlim = c(-0.32, 0.32), 
     ylim = c(-0.32, 0.32), 
     axes = F, 
     xlab = "", 
     ylab = ""
)  
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- x %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids <- centroids[c(2, 1, 4, 3), ]
points(PC2 ~ PC1, centroids, pch = 3, cex = 3, col = plot_cols, lwd = 2)
points(x[metadata$sample == "4899", pc1], x[metadata$sample == "4899", pc2], col = plot_cols[1], pch = 17)
points(x[metadata$sample == "2099", pc1], x[metadata$sample == "2099", pc2], col = plot_cols[2], pch = 17)
points(x[metadata$sample == "COL", pc1], x[metadata$sample == "COL", pc2], col = plot_cols[3], pch = 16)
points(x[metadata$sample == "9999SR", pc1], x[metadata$sample == "9999SR", pc2], col = plot_cols[4], pch = 15)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 0, cex = 1.5, text = "PC1 scores (59.0%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 scores (26.6%)")
legend("topleft", legend = c("4899", "2099", "COL", "9999SR"), pch = c(17, 17, 16, 15), col = plot_cols, bty = "n", cex = 1.5)
box()

(new_xlim = xlim * ratio)

(new_ylim = ylim * ratio)

nearest(wavenumbers, 1515)
nearest(wavenumbers, 1604)
nearest(wavenumbers, 1378)
nearest(wavenumbers, 1168)

lutzke_indices <- c(115, 224, 295, 341)

axis_lim <- 1.15
plot(y, axes = FALSE, type = "n", 
     xlim = new_xlim, 
     ylim = new_ylim, 
     xlab = "", 
     ylab = "", 
     #xlim = c(-axis_lim, axis_lim),
     #ylim = c(-axis_lim, axis_lim)
)

abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 0, cex = 1.5, text = "PC1 loadings (59.0%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 loadings (26.6%)")

text(y[lutzke_indices, ], labels = round(strip_non_numbers(ylabs[lutzke_indices])), col = 1, cex = 0.8)  # This just prints the species
arrow.len <- 0.1
arrows(0, 0, y[lutzke_indices, 1L] * 0.8, y[lutzke_indices, 2L] * 0.8, length = arrow.len, col = "#00000033")

#non_sporopollenin_indices <- sapply(c(1642, 1022, 1541, 1558, 1455, 1750, 1772, 1080), function(x){which.min(abs(x - wavenumbers))})
#text(y[non_sporopollenin_indices, ], labels = round(strip_non_numbers(ylabs[non_sporopollenin_indices])), cex = 0.8, col = "red")  # This just prints the species
#arrows(0, 0, y[non_sporopollenin_indices, 1L] * 0.8, y[non_sporopollenin_indices, 2L] * 0.8, length = arrow.len, col = "#FF000033")
#legend("topleft", lty = 1, col = 1:2, legend = c("Sporopollenin", "Other (amide, etc.)"), cex = 1.5, bty = "n")
text(1.2, 1.5, "Plot ref: XXXX")
box()

par(mfrow = c(1, 2), pty = "m")
plot(wavenumbers[1:440], apca$rotation[, 1], xlim = rev(range(wavenumbers[1:440])), lwd = 2, xlab = lab("wn"), ylab = "PC1 loadings", axes = F, type = "n")
abline(h = 0)
# abline(v = 1168)
# abline(v = 1515)
# abline(v = 1605)
# abline(v = 1025)
abline(h = sqrt(1/440), col = "#00009999", lwd = 2, lty = 2)
abline(h = -sqrt(1/440), col = "#00009999", lwd = 2, lty = 2)
lines(wavenumbers[1:440], apca$rotation[, 1], lwd = 2, col = "darkgrey")
loadings1_threshold <- which(apca$rotation[, 1] > sqrt(1/440) |
                               apca$rotation[, 1] < -sqrt(1/440))

internal_axis_ticks(1:2, 1:2)
loadings1_peaks <- find_peaks(apca$rotation[, 1], minpeakheight = sqrt(1/440))
text(loadings1_peaks$peak_wavenumber, loadings1_peaks$abs, round(loadings1_peaks$peak_wavenumber), cex = 0.9)
#abline(v = wavenumbers[loadings1_peaks[, 2]])

plot(wavenumbers[1:440], apca$rotation[, 2], type = "n", xlim = rev(range(wavenumbers[1:440])), lwd = 2, axes = F, ylab = "PC2 loadings", xlab = lab("wn"))
abline(h = 0)
abline(h = sqrt(1/440), col = "#00009999", lwd = 2, lty = 2)
abline(h = -sqrt(1/440), col = "#00009999", lwd = 2, lty = 2)
lines(wavenumbers[1:440], apca$rotation[, 2], lwd = 2, col = "darkgrey")
# abline(v = 1168)
# abline(v = 1515)
# abline(v = 1605)
which(apca$rotation[, 2] > sqrt(1/440) |
        apca$rotation[, 2] < -sqrt(1/440))

loadings2_threshold <- strip_non_numbers(names(which(apca$rotation[, 2] > sqrt(1/440) |
                                                       apca$rotation[, 2] < -sqrt(1/440))))

loadings2_peaks <- find_peaks(apca$rotation[, 2], minpeakheight = sqrt(1/440))
text(loadings2_peaks$peak_wavenumber, loadings2_peaks$abs, round(loadings2_peaks$peak_wavenumber), cex = 0.9)
internal_axis_ticks(1:2, 1:2)
#abline(v = wavenumbers[loadings2_peaks[, 2]])

############################ 2nd der

par(mfrow = c(1,2))                            # Splitting the plot in 2.
#biplot(apca)                                    # In-built biplot() R func.
PCA <- prcomp(-data.frame(t(diff(t(spectra[, 1:440]), differences = 2))), center = T, scale. = F)

choices = 1:2                                  # Selecting first two PC's
scale = 1                                      # Default
scores = PCA$x                                  # The scores
lam = PCA$sdev[choices]                        # Sqrt e-vals (lambda) 2 PC's
n = nrow(scores)                               # no. rows scores
lam = lam * sqrt(n)                            # See below.


x = t(t(scores[,choices]) / lam)         # scaled scores
y = t(t(PCA$rotation[,choices]) * lam)         # scaled eigenvecs (loadings)

n = nrow(x)                                    # Same as dataset (150)
p = nrow(y)                                    # Three var -> 3 rows

xlabs = 1L:n
xlabs = as.character(xlabs)                    # no. from 1 to 150 
dimnames(x) = list(xlabs, dimnames(x)[[2L]])   # no's and PC1 / PC2

ylabs = dimnames(y)[[1L]]                      # Iris species
ylabs = as.character(ylabs)
dimnames(y) <- list(ylabs, dimnames(y)[[2L]])  # Species and PC1/PC2


unsigned.range = function(x) c(-abs(min(x, na.rm = TRUE)), 
                               abs(max(x, na.rm = TRUE)))
rangx1 = unsigned.range(x[, 1L])               # Range first col x
rangx2 = unsigned.range(x[, 2L])               # Range second col x
rangy1 = unsigned.range(y[, 1L])               # Range 1st scaled evec
rangy2 = unsigned.range(y[, 2L])               # Range 2nd scaled evec


(xlim = ylim = rangx1 = rangx2 = range(rangx1, rangx2))

(ratio = max(rangy1/rangx1, rangy2/rangx2)) 

par(mfrow = c(1, 2), pty = "s")                                
pc1 <- 1
pc2 <- 2
plot(x, type = "n", 
     #xlim = xlim, 
     #ylim = ylim, 
     xlim = c(-0.32, 0.32), 
     ylim = c(-0.32, 0.32), 
     axes = F, 
     xlab = "", 
     ylab = ""
)  
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- x %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids <- centroids[c(2, 1, 4, 3), ]
points(PC2 ~ PC1, centroids, pch = 3, cex = 3, col = plot_cols, lwd = 2)
points(x[metadata$sample == "4899", pc1], x[metadata$sample == "4899", pc2], col = plot_cols[1], pch = 17)
points(x[metadata$sample == "2099", pc1], x[metadata$sample == "2099", pc2], col = plot_cols[2], pch = 17)
points(x[metadata$sample == "COL", pc1], x[metadata$sample == "COL", pc2], col = plot_cols[3], pch = 16)
points(x[metadata$sample == "9999SR", pc1], x[metadata$sample == "9999SR", pc2], col = plot_cols[4], pch = 15)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 0, cex = 1.5, text = "PC1 scores (82.4%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 scores (6.1%)")
legend("topleft", legend = c("4899", "2099", "COL", "9999SR"), pch = c(17, 17, 16, 15), col = plot_cols, bty = "n", cex = 1.5)
box()

(new_xlim = xlim * ratio)

(new_ylim = ylim * ratio)

lutzke_indices <- which(wavenumbers %in% lutzke_wavenumbers_only[lutzke_wavenumbers_only < 1800]) - 2

axis_lim <- 1.15
plot(y, axes = FALSE, type = "n", 
     xlim = new_xlim, 
     ylim = new_ylim, 
     #ylim = c(-0.1, 0.1),
     xlab = "", 
     ylab = "", 
     #xlim = c(-axis_lim, axis_lim),
     #ylim = c(-axis_lim, axis_lim)
)

abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 0, cex = 1.5, text = "PC1 loadings (82.4%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 loadings (6.1%)")

text(y[lutzke_indices, ], labels = round(strip_non_numbers(ylabs[lutzke_indices])), col = 1, cex = 0.8)  # This just prints the species
arrow.len <- 0.1
arrows(0, 0, y[lutzke_indices, 1L] * 0.8, y[lutzke_indices, 2L] * 0.8, length = arrow.len, col = "#00000033")

non_sporopollenin_indices <- sapply(c(1642, 1022, 1541, 1558, 1455, 1750, 1772, 1080), function(x){which.min(abs(x - wavenumbers))})
text(y[non_sporopollenin_indices, ], labels = round(strip_non_numbers(ylabs[non_sporopollenin_indices])), cex = 0.8, col = "red")  # This just prints the species
arrows(0, 0, y[non_sporopollenin_indices, 1L] * 0.8, y[non_sporopollenin_indices, 2L] * 0.8, length = arrow.len, col = "#FF000033")
legend("topleft", lty = 1, col = 1:2, legend = c("Sporopollenin", "Other (amide, etc.)"), cex = 1.5, bty = "n")
text(0.2, 0.24, "Plot ref: FSG8")
box()

par(mfrow = c(2, 1), mar = c(5.1, 6.1, 4.1, 0.6), pty = "m")
plot(wavenumbers[2:439], PCA$rotation[, 1], xlim = rev(range(wavenumbers[2:439])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[2:439], PCA$rotation[, 1], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
#mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(PCA$rotation[, pc1]) + ((max(PCA$rotation[, pc1]) - min(PCA$rotation[, pc1])) / 2), cex = 1.5, text = "PC1 loadings")
box()
lutzke_wavenumbers()

#par(mfrow = c(1, 1), mar = c(5.1, 6.1, 4.1, 0.1))
plot(wavenumbers[2:439], PCA$rotation[, 2], xlim = rev(range(wavenumbers[2:439])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[2:439], PCA$rotation[, 2], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(PCA$rotation[, 2]) + ((max(PCA$rotation[, 2]) - min(PCA$rotation[, 2])) / 2), cex = 1.5, text = "PC2 loadings")
box()
lutzke_wavenumbers()
mtext(side = 1, line = 3, at = 960, text = "Plot ref: K7SO")


peaks <- find_peaks(spec_diff)
peaks2 <-  find_peaks(spec_diff2)
peaks3 <-  find_peaks(spec_diff3)

intersecting_wavenumbers1 <- intersect(peaks$peak_wavenumber, peaks2$peak_wavenumber)
intersecting_wavenumbers2 <- intersect(peaks$peak_wavenumber, peaks3$peak_wavenumber)
#intersecting_wavenumbers <- intersect(intersecting_wavenumbers1, intersecting_wavenumbers2) %>% sort()
#intersecting_wavenumbers <- unique(intersecting_wavenumbers1, intersecting_wavenumbers2) %>% sort()
intersecting_wavenumbers <- unique(c(peaks$peak_wavenumber, peaks2$peak_wavenumber, peaks3$peak_wavenumber)) %>% sort()

opposite_wavenumbers <- wavenumbers[which(spec_diff < 0 & spec_diff2 > 0 & spec_diff3 > 0)]
opposite_wavenumbers <- sort(unique(c(opposite_wavenumbers, wavenumbers[which(spec_diff > 0 & spec_diff2 < 0 & spec_diff3 < 0)])))
intersecting_wavenumbers <- intersect(intersecting_wavenumbers, opposite_wavenumbers)

par(mfrow = c(1, 1))
plot(wavenumbers, spec_diff, xlim = c(1800, 948), type = "l")
points(wavenumbers, spec_diff)
#abline(v = peaks[peaks$peak_type == "peak", ]$peak_wavenumber, col = "red")
#abline(v = peaks[peaks$peak_type == "trough", ]$peak_wavenumber, col = "blue")
abline(h = 0)
lines(wavenumbers, spec_diff2)
lines(wavenumbers, spec_diff3)
#abline(v = intersecting_wavenumbers2, col = "red")
#abline(v = intersecting_wavenumbers, col = "red")

interesting_wavenumbers <- c(
  1024.020,
  1078.010,
  1168.650,
  1232.290,
  1261.220,
  1290.140,
  1436.710,
  1513.85,
  1604.480,
  1616.060,
  1710.550
)

lutzke_wavenumbers_only <- c(1605, 1588, 1514, 3026, 1202, 1168, 833, 1680, 1630, 985, 2926, 2855, 1463, 720, 1411, 1710, 1437, 1283, 943, 1125, 1103, 3350, 1323)
lutzke_wavenumbers_only <- lutzke_wavenumbers_only[lutzke_wavenumbers_only > 950]
lutzke_wavenumbers_only <- sapply(lutzke_wavenumbers_only, function(x){wavenumbers[which.min(abs(x - wavenumbers))]}) %>% sort()

abline(v = interesting_wavenumbers)

interesting_indices <- which(wavenumbers %in% interesting_wavenumbers)


###############################################
#####   COMPLETE MESS - A LOAD OF CODE FOUND IN RICH BARCLAY'S SCRIPT



mtext(side = 2, line = 2, at = min(spectra) + ((max(spectra) - min(spectra)) / 2), cex = 2, text = "Relative absorbance")
mtext(side = 1, line = 3.5, at = 948 + ((1800 - 948) / 2), cex = 2, text = expression(paste("Wavenumber (cm"^-1*")")))
#internal_axis_ticks(1, 1)
legend("bottomright", legend = c("9999SR", "COL", "2099", "4899"), lty = 1, col = rev(plot_cols), bg = "white", lwd = 2, cex = 1.5)


par(mfrow = c(2, 1), mar = c(5.1, 4.1, 0.6, 0.6), pty = "m")
plot(wavenumbers, spectra[1, ], xlim = rev(range(wavenumbers)), type = "n", ylim = range(spectra), axes = F, xlab = "", ylab = "")
rect(1800, 10, 950, -10, col = "#00000016", border = NA)
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "9999SR"), ]), col = plot_cols[4], lwd = 2)  
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "COL"), ]), col = plot_cols[3], lwd = 2)  
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "4899"), ]), col = plot_cols[1], lwd = 2)  
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "2099"), ]), col = plot_cols[2], lwd = 2)  
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 2, line = 2, at = min(spectra) + ((max(spectra) - min(spectra)) / 2), cex = 1.5, text = "Relative absorbance")
#internal_axis_ticks(1, 1)
#legend("bottomright", legend = c("9999SR", "COL", "2099", "4899"), lty = 1, col = rev(plot_cols), bty = "n", lwd = 2, cex = 1.5)
box()

plot(wavenumbers, spectra[1, ], xlim = c(1800, 948), type = "n", ylim = range(spectra), axes = F, xlab = "", ylab = "")
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "9999SR"), ]), col = plot_cols[4], lwd = 2)  
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "COL"), ]), col = plot_cols[3], lwd = 2)  
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "4899"), ]), col = plot_cols[1], lwd = 2)  
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "2099"), ]), col = plot_cols[2], lwd = 2)  
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 2, line = 2, at = min(spectra) + ((max(spectra) - min(spectra)) / 2), cex = 1.5, text = "Relative absorbance")
mtext(side = 1, line = 3.5, at = 948 + ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
#internal_axis_ticks(1, 1)
legend("bottomright", legend = c("9999SR", "COL", "2099", "4899"), lty = 1, col = rev(plot_cols), bty = "n", lwd = 2, cex = 1.5)
box()

# colMedians <- function(x){
#   apply(x, 2, median)
# }

ggbiplot::ggbiplot(apca, circle = T, groups = metadata$sample, ellipse = T)

t.table <- function(wn){
  wavenumber <- wn
  wn <- which.min(abs(wn - wavenumbers))
  ln1 <- broom::tidy(t.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "9999SR"), wn]))
  ln2 <- broom::tidy(t.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "4899"), wn]))
  ln3 <- broom::tidy(t.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "2099"), wn]))
  results <- rbind(ln1, ln2, ln3) %>% data.frame()
  results$greater_than <- results$estimate1 < results$estimate2
  results$wn <- wavenumber
  results$prob_degree <- NA
  var1 <- broom::tidy(var.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "9999SR"), wn]))
  var2 <- broom::tidy(var.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "4899"), wn]))
  var3 <- broom::tidy(var.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "2099"), wn]))
  for(i in 1:nrow(results)){
    if(results$p.value[i] < 0.001){
      results$prob_degree[i] <- "***"
    } else if(results$p.value[i] < 0.01){
      results$prob_degree[i] <- "**"
    } else if(results$p.value[i] < 0.05){
      results$prob_degree[i] <- "*"
    } else {
      results$prob_degree[i] <- ""
    }
  }
  results
}

wilcox.table <- function(wn){
  wavenumber <- wn
  wn <- which.min(abs(wn - wavenumbers))
  ln1 <- broom::tidy(wilcox.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "9999SR"), wn]))
  ln2 <- broom::tidy(wilcox.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "4899"), wn]))
  ln3 <- broom::tidy(wilcox.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "2099"), wn]))
  results <- rbind(ln1, ln2, ln3) %>% data.frame()
  #results$greater_than <- results$estimate1 < results$estimate2
  results$wn <- wavenumber
  results$prob_degree <- NA
  var1 <- broom::tidy(var.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "9999SR"), wn]))
  var2 <- broom::tidy(var.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "4899"), wn]))
  var3 <- broom::tidy(var.test(spectra[which(metadata$sample == "COL"), wn], spectra[which(metadata$sample == "2099"), wn]))
  for(i in 1:nrow(results)){
    if(results$p.value[i] < 0.001){
      results$prob_degree[i] <- "***"
    } else if(results$p.value[i] < 0.01){
      results$prob_degree[i] <- "**"
    } else if(results$p.value[i] < 0.05){
      results$prob_degree[i] <- "*"
    } else if(results$p.value[i] < 0.1){
      results$prob_degree[i] <- "."
    } else {
      results$prob_degree[i] <- ""
    }
  }
  results
}


(t_test_results <- lapply(c(3026, 1202, 1168, 1605, 1588, 1515), t.table))
(wilcox_results <- lapply(c(3026, 1202, 1168, 1605, 1588, 1515, 1680, 1630, 985), wilcox.table))

spectra[, which.min(abs(1605 - wavenumbers))] 


wn <- 3026
spectra_box <- data.frame(sample = metadata$sample, 
                          abs = spectra[, which.min(abs(wn - wavenumbers))],
                          wn = wn)
spectra_box$sample <- factor(spectra_box$sample, levels = factor_levels)


wn <- 1202
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))

wn <- 1168
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))

wn <- 1605
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))


wn <- 1588
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))

wn <- 1515
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))

wn <- 1680
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))
wn <- 1630
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))
wn <- 985
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))

wn <- 1585
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))

wn <- 1550
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))

wn <- 1465
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))

wn <- 1378
spectra_box <- rbind(spectra_box, data.frame(sample = metadata$sample, 
                                             abs = spectra[, which.min(abs(wn - wavenumbers))],
                                             wn = wn))

spectra_box$myb <- ""
spectra_box$myb[spectra_box$sample == "COL"] <- "CONTROL"
spectra_box$myb[spectra_box$sample %in% c("995", "9999SR", "MLP8")] <- "MYB99-"
spectra_box$myb[spectra_box$sample %in% c("991", "4899", "2099")] <- "MYB99+"
spectra_box$myb <- factor(spectra_box$myb, levels = c("MYB99-", "CONTROL", "MYB99+"))

#spectra_box2 <- spectra_box[spectra_box$wn %in% c(1168, 1515, 1550, 1585, 1605), ]
#spectra_box2 <- spectra_box[spectra_box$wn %in% c(1168, 1465, 1515, 1605), ]
spectra_box2 <- spectra_box

spectra_box2 %>% ggplot(aes(forcats::fct_reorder(myb, abs), abs, fill = myb)) + geom_boxplot(notch = T) + 
  labs(x = "", y = "Relative absorbance") + 
  facet_wrap(~wn, nrow = 1) +
  #theme_classic() +
  scale_fill_manual(values = rev(plot_cols[c(1, 3, 4)])) +
  theme_classic()

par(mfrow = c(1, 3), mar = c(4.1, 5.1, 1.1, 1.6))
plot(abs ~ myb, spectra_box2[spectra_box2$wn == 1605, ], col = rev(plot_cols[c(1, 3, 4)]), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylim = c(0.8, 2.2), ylab = expression(paste("Relative absorbance @ 1605 cm "^-1*"")), cex.lab = 2, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ myb, spectra_box2[spectra_box2$wn == 1605, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 2)
text(1:3, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$myb))), cex = 1.5)
internal_axis_ticks(2, 2, cex.axis = 1.5)
mtext(side = 1, line = 1.5, at = 1:3, text = c("MYB99-", "CONTROL", "MYB99+"), cex = 1.2)
plot(abs ~ myb, spectra_box2[spectra_box2$wn == 1515, ], col = rev(plot_cols[c(1, 3, 4)]), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylim = c(0.8, 2.2), ylab = expression(paste("Relative absorbance @ 1515 cm "^-1*"")), cex.lab = 2, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ myb, spectra_box2[spectra_box2$wn == 1515, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 2)
text(1:3, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$myb))), cex = 1.5)
internal_axis_ticks(2, 2, cex.axis = 1.5)
mtext(side = 1, line = 1.5, at = 1:3, text = c("MYB99-", "CONTROL", "MYB99+"), cex = 1.2)
plot(abs ~ myb, spectra_box2[spectra_box2$wn == 1168, ], col = rev(plot_cols[c(1, 3, 4)]), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylim = c(0.8, 2.2), ylab = expression(paste("Relative absorbance @ 1168 cm "^-1*"")), cex.lab = 2, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ myb, spectra_box2[spectra_box2$wn == 1168, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 2)
text(1:3, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$myb))), cex = 1.5)
internal_axis_ticks(2, 2, cex.axis = 1.5)
mtext(side = 1, line = 1.5, at = 1:3, text = c("MYB99-", "CONTROL", "MYB99+"), cex = 1.2)
text(3, 0.82, "Plot ref: 8BBZ", cex = 1.5)

stretch_axis <- function(vector, pc, at = "hi"){
  stopifnot(at %in% c("lo", "hi", "both"))
  lo <- range(vector)[1]
  hi <- range(vector)[2]
  if(at == "hi"){
    c(lo, hi + ((hi - lo) / 100 * pc))
  } else if(at == "lo"){
    c(lo - ((hi - lo) / 100 * pc), hi)
  } else if(at == "both"){
    c(lo - ((hi - lo) / 100 * pc), hi + ((hi - lo) / 100 * pc))
  }
}

plot_myb_box <- function(df, wn){
  new_lim <- stretch_axis(df[df$wn == wn, ]$abs, 5)
  par(mfrow = c(1, 1), mar = c(4.1, 5.1, 1.1, 1.6))
  plot(abs ~ sample, df[df$wn == wn, ], ylim = new_lim, col = plot_cols, staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylab = substitute(paste("Relative absorbance @ ", the_wn, " cm "^-1*""), list(the_wn = wn)), cex.lab = 2, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
  stripchart(abs ~ sample, df[df$wn == wn, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 2)
  text(1:length(factor_levels), new_lim[2], unname(table(df[df$wn == wn, ]$sample)), cex = 1.5)
  internal_axis_ticks(2, 2, cex.axis = 1.5)
  mtext(side = 1, line = 1.5, at = 1:length(factor_levels), text = factor_levels, cex = 1.2)
}

spectra_box3 <- spectra_box2[!spectra_box2$sample %in% c("991", "995"), ]
spectra_box3$sample <- droplevels(spectra_box3$sample)
plot_myb_box(spectra_box2, 1605)
plot_myb_box(spectra_box2, 1515)
plot_myb_box(spectra_box2, 1465)
plot_myb_box(spectra_box2, 1378)
plot_myb_box(spectra_box2, 1168)

spectra_box2[spectra_box2$wn == 1605, ] / spectra_box2[spectra_box2$wn == 1515, ]

plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1515, ], col = plot_cols, staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylab = expression(paste("Relative absorbance @ 1515 cm "^-1*"")), cex.lab = 2, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1515, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 2)
text(1:length(factor_levels), 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.5)
internal_axis_ticks(2, 2, cex.axis = 1.5)
mtext(side = 1, line = 1.5, at = 1:length(factor_levels), text = factor_levels, cex = 1.2)


plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1515, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylim = c(0.8, 2.2), ylab = expression(paste("Relative absorbance @ 1515 cm "^-1*"")), cex.lab = 2, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1515, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 2)
text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.5)
internal_axis_ticks(2, 2, cex.axis = 1.5)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 1.2)
plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1168, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylim = c(0.8, 2.2), ylab = expression(paste("Relative absorbance @ 1168 cm "^-1*"")), cex.lab = 2, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1168, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 2)
text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.5)
internal_axis_ticks(2, 2, cex.axis = 1.5)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 1.2)
text(3.9, 0.82, "Plot ref: 1UQS", cex = 1.5)

par(mfrow = c(1, 4), mar = c(4.1, 5.1, 1.1, 1.6))
plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1605, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylim = c(0.8, 2.2), ylab = expression(paste("Relative absorbance @ 1605 cm "^-1*"")), cex.lab = 1.5, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1605, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
text(1:4, 2.2, unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample)), cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1515, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylim = c(0.8, 2.2), ylab = expression(paste("Relative absorbance @ 1515 cm "^-1*"")), cex.lab = 1.5, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1515, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
text(1:4, 2.2, unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample)), cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1465, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylim = c(0.8, 2.2), ylab = expression(paste("Relative absorbance @ 1465 cm "^-1*"")), cex.lab = 1.5, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1465, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
text(1:4, 2.2, unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample)), cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1168, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylim = c(0.8, 2.2), ylab = expression(paste("Relative absorbance @ 1168 cm "^-1*"")), cex.lab = 1.5, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F)
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1168, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
text(1:4, 2.2, unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample)), cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
text(3.9, 0.82, "Plot ref: JF7U", cex = 0.9)

#1378
spectra_box2 <- spectra_box[spectra_box$wn %in% c(1168, 1378, 1515, 1605), ]

par(mfrow = c(1, 4), mar = c(4.1, 5.1, 1.1, 1.6))
plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1605, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylab = expression(paste("Relative absorbance @ 1605 cm "^-1*"")), cex.lab = 1.5, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F, ylim = c(0.8, 2.2))
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1605, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
text(1:4, 2.2, unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample)), cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
mtext(side = 2, line = 2.75, at = 2.2, las = 2, text = "b)", cex = 2)
plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1515, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylab = expression(paste("Relative absorbance @ 1515 cm "^-1*"")), cex.lab = 1.5, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F, ylim = c(0.8, 2.2))
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1515, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
text(1:4, 2.2, unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample)), cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
mtext(side = 2, line = 2.75, at = 2.2, las = 2, text = "c)", cex = 2)
plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1378, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylab = expression(paste("Relative absorbance @ 1378 cm "^-1*"")), cex.lab = 1.5, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F, ylim = c(0.8, 1.15))
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1378, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
text(1:4, 1.15, unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample)), cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
mtext(side = 2, line = 2.75, at = 1.15, las = 2, text = "d)", cex = 2)
plot(abs ~ sample, spectra_box2[spectra_box2$wn == 1168, ], col = rev(plot_cols), staplelty = 0, whisklty = 1, notch = T, las = 1, xlab = "", ylab = expression(paste("Relative absorbance @ 1168 cm "^-1*"")), cex.lab = 1.5, cex.axis = 1.3, axes = F, boxwex = 0.6, outline = F, ylim = c(0.8, 1.15))
stripchart(abs ~ sample, spectra_box2[spectra_box2$wn == 1168, ], method = "jitter", col = "#00000099", add = T, vertical = T, pch = 1, cex = 1.5)
#text(1:4, 2.2, paste0("n = ", unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample))), cex = 1.2)
text(1:4, 1.15, unname(table(spectra_box2[spectra_box2$wn == 1168, ]$sample)), cex = 1.2)
internal_axis_ticks(2, 2, cex.axis = 1.3)
mtext(side = 1, line = 1.5, at = 1:4, text = c("9999SR", "COL", "4899", "2099"), cex = 0.8)
mtext(side = 2, line = 2.75, at = 1.15, las = 2, text = "e)", cex = 2)
text(3.9, 0.82, "Plot ref: JF7U", cex = 0.9)


summary(aov(abs ~ myb, spectra_box2[spectra_box2$wn == 1605, ]))
summary(aov(abs ~ myb, spectra_box2[spectra_box2$wn == 1515, ]))
summary(aov(abs ~ myb, spectra_box2[spectra_box2$wn == 1168, ]))

summary(aov(abs ~ myb, spectra_box[spectra_box$wn == 1550, ]))
summary(aov(abs ~ myb, spectra_box[spectra_box$wn == 1585, ]))
summary(aov(abs ~ myb, spectra_box[spectra_box$wn == 1202, ]))
summary(aov(abs ~ myb, spectra_box[spectra_box$wn == 3026, ]))

summary(aov(abs ~ myb, spectra_box[spectra_box$wn == 1680, ]))
summary(aov(abs ~ myb, spectra_box[spectra_box$wn == 1630, ]))
summary(aov(abs ~ myb, spectra_box[spectra_box$wn == 985, ]))


t_test_results_9999SR <- numeric(ncol(spectra))
for(i in 1:ncol(spectra)){
  t_test <- t.test(spectra[which(metadata$sample == "COL"), i], spectra[which(metadata$sample == "9999SR"), i], var.equal = input.var.equal)
  t_test_results_9999SR[i] <- t_test$p.value
}

t_test_results_4899 <- numeric(ncol(spectra))
for(i in 1:ncol(spectra)){
  t_test <- t.test(spectra[which(metadata$sample == "COL"), i], spectra[which(metadata$sample == "4899"), i], var.equal = input.var.equal)
  t_test_results_4899[i] <- t_test$p.value
}

t_test_results_2099 <- numeric(ncol(spectra))
for(i in 1:ncol(spectra)){
  t_test <- t.test(spectra[which(metadata$sample == "COL"), i], spectra[which(metadata$sample == "2099"), i], var.equal = input.var.equal)
  t_test_results_2099[i] <- t_test$p.value
}

t_test_results_4899_9999SR <- numeric(ncol(spectra))
for(i in 1:ncol(spectra)){
  t_test <- t.test(spectra[which(metadata$sample == "9999SR"), i], spectra[which(metadata$sample == "4899"), i], var.equal = input.var.equal)
  t_test_results_4899_9999SR[i] <- t_test$p.value
}

sum(t_test_results_4899_9999SR < 0.05)

plot(wavenumbers, t_test_results_2099, type = "l", col = plot_cols[2], ylim = c(0.05, 0), xlim = c(1800, 948), lwd = 2)
lines(wavenumbers, t_test_results_4899, type = "l", col = plot_cols[1], lwd = 2)
lines(wavenumbers, t_test_results_9999SR, type = "l", col = plot_cols[4], lwd = 2)
abline(v = c(1710, 1605, 1514))

spec_diff <- colMeans(spectra[which(metadata$sample == "9999SR"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff2 <- colMeans(spectra[which(metadata$sample == "4899"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
spec_diff3 <- colMeans(spectra[which(metadata$sample == "2099"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])


#9999SR < COL < MYB99+
a <- which(colMeans(spectra[which(metadata$sample == "9999SR"), ]) < colMeans(spectra[which(metadata$sample == "COL"), ]))
b <- which(colMeans(spectra[which(metadata$sample == "COL"), ]) < colMeans(rbind(colMeans(spectra[which(metadata$sample == "2099"), ]), colMeans(spectra[which(metadata$sample == "4899"), ]))))
indices <- intersect(a, b)
continuum <- seek_seq(indices)
continuum$from <- wavenumbers[continuum$from]
continuum$to <- wavenumbers[continuum$to]
plot(wavenumbers, spectra[1, ], xlim = rev(range(wavenumbers)), type = "n", ylim = range(spectra), axes = F, xlab = "", ylab = "")
rect(continuum$from, 10, continuum$to, -10, col = "#00000016", border = NA)
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "9999SR"), ]), col = plot_cols[4], lwd = 2)  
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "COL"), ]), col = plot_cols[3], lwd = 2)  
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "4899"), ]), col = plot_cols[1], lwd = 2)  
lines(wavenumbers, colMeans(spectra[which(metadata$sample == "2099"), ]), col = plot_cols[2], lwd = 2)  
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 2, line = 2, at = min(spectra) + ((max(spectra) - min(spectra)) / 2), cex = 1.5, text = "Relative absorbance")
#internal_axis_ticks(1, 1)
#legend("bottomright", legend = c("9999SR", "COL", "2099", "4899"), lty = 1, col = rev(plot_cols), bty = "n", lwd = 2, cex = 1.5)
lutzke_wavenumbers()
box()



wavenumbers[indices]

#spec_diff <- colMedians(spectra[which(metadata$sample == "9999SR"), ]) - colMedians(spectra[which(metadata$sample == "COL"), ])
#spec_diff2 <- colMedians(spectra[which(metadata$sample == "4899"), ]) - colMedians(spectra[which(metadata$sample == "COL"), ])
#spec_diff3 <- colMedians(spectra[which(metadata$sample == "2099"), ]) - colMedians(spectra[which(metadata$sample == "COL"), ])

seek_table2 <- seek_seq(which(spec_diff > 0 & spec_diff2 < 0 & spec_diff3 < 0))
seek_table2$from <- wavenumbers[seek_table2$from]
seek_table2$to <- wavenumbers[seek_table2$to]


seek_table3 <- seek_seq(which(spec_diff < 0 & spec_diff2 > 0 & spec_diff3 > 0))
seek_table3$from <- wavenumbers[seek_table3$from]
seek_table3$to <- wavenumbers[seek_table3$to]

par(pty = "m", mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 0.6))
plot(wavenumbers, spec_diff, xlim = c(1800, 948), 
     ylim = c(-0.23, 0.29),
     type = "n", axes = F, xlab = "", ylab = "")
rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(v = c(1710, 1605, 1514, 1436, 1168), lty = 2, col = "#000000DD")
#lutzke_wavenumbers()
lines(wavenumbers, spec_diff, lwd = 2, col = plot_cols[4])
lines(wavenumbers, spec_diff2, lwd = 2, col = plot_cols[1])
lines(wavenumbers, spec_diff3, lwd = 2, col = plot_cols[2])
#points(wavenumbers[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
abline(h = 0, lty = 1, lwd = 2, col = plot_cols[3])
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, at = 0, label = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Relative difference (sample - control)")
box()
abline(v = c(1515, 1605, 1588))
abline(v = c(1202, 1168), lty = 2)
abline(v = c(1680, 1630, 985), lty = 3)
#text(1710, 0.15, "***", cex = 2)
#text(1710, -0.19, "***", cex = 2)
#text(1605, 0.140, "***", cex = 3)
#text(1605, 0.125, "***", cex = 3)
#text(1605, -0.23, "***", cex = 3)
text(1605, 0.140, "***", cex = 2.5, col = plot_cols[1])
text(1605, 0.125, "***", cex = 2.5, col = plot_cols[2])
text(1605, -0.23, "***", cex = 2.5, col = plot_cols[4])

text(1515, 0.18, "***", cex = 2)
text(1515, -0.14, "*", cex = 2)
#text(1436, 0.125, "NS", cex = 1.5)
#text(1436, -0.23, "***", cex = 2)
text(1168, 0.09, "***", cex = 2)
text(1168, -0.11, "NS", cex = 1.5)

# mtext(side = 3, line = 0.5, at = 1710, text = "C=O")
# mtext(side = 3, line = 0.5, at = 1605, text = "C=C")
# mtext(side = 3, line = 0.5, at = 1514, text = "C=C")
# mtext(side = 3, line = 0.5, at = 1436, text = "dC-H")
#axis(2)
#rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
#rect(900, 0, 940, -10, col = plot_cols[3], border = NA)
#text(1050, max(spec_diff), "Plot ref: OQZU")
#phil_wavenumbers()
legend("topright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
text(965, -0.22, "Plot ref: R43X")
#abline(v = interesting_wavenumbers, lwd = 2, col = "#00000077")


par(mfrow = c(1, 1))
ylims <- c(-0.35, 0.025)
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 900), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(h = 0, lty = 1, col = "black")
lutzke_wavenumbers(cex = 0.8)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[1], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2) - 0.1, col = plot_cols[2], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - 0.2, col = plot_cols[3], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - 0.3, col = plot_cols[4], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(940, 0, "4899", cex = 1.5, adj = c(0, 0.5))
text(940, -0.1, "2099", cex = 1.5, adj = c(0, 0.5))
text(940, -0.2, "COL", cex = 1.5, adj = c(0, 0.5))
text(940, -0.3, "9999SR", cex = 1.5, adj = c(0, 0.5))
#legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
text(930, -0.335, "Plot ref: KZZH")
box()

par(mfrow = c(1, 1))
ylims <- c(-0.075, 0.05)
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 900), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(h = 0, lty = 1, col = "black")
lutzke_wavenumbers(cex = 0.8)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[1], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2), col = plot_cols[2], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), col = plot_cols[3], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), col = plot_cols[4], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(940, 0, "4899", cex = 1.5, adj = c(0, 0.5))
text(940, -0.1, "2099", cex = 1.5, adj = c(0, 0.5))
text(940, -0.2, "COL", cex = 1.5, adj = c(0, 0.5))
text(940, -0.3, "9999SR", cex = 1.5, adj = c(0, 0.5))
#legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
text(930, -0.335, "Plot ref: KZZH")
box()

# diffs_wns <- 
#   seek_seq(
#     which(
#       diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) < 0 & (
#         !(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2) < 0) &
#           !(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2) < 0)
#       )
#     )
#   )


diffs_wns <- 
  seek_seq(
    which(
      (diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) > 0) & 
        (diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2) > 0) &
        (diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2) > 0)
    )
  )

diffs_wns$from <- wavenumbers[diffs_wns$from]
diffs_wns$to <- wavenumbers[diffs_wns$to]


yaxlim = 0.08
par(mfrow = c(1, 3), mar = c(5.6, 4.1, 3.6, 0.6))
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), 
     xlim = c(3050, 2850), ylim = c(-yaxlim, yaxlim - 0.03),
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(h = 0, lty = 1, col = "black")
abline(h = 0)
#lutzke_wavenumbers(cex = 0.8, only = 3026)
rect(diffs_wns$from, 1, diffs_wns$to, -1, col = "#00000055", border = NA)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2), col = plot_cols[1], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[2], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), col = plot_cols[3], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), col = plot_cols[4], lwd = 3)
axis(1, cex.axis = 2, tck = 0.01)
#mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
ylims <- c(-yaxlim, yaxlim - 0.03)
mtext(side = 2, line = 1.5, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
abline(v = 3026, col = "black", lty = 2, lwd = 2)
mtext("v(CH)", side = 3, line = 0.3,  at = 3026, cex = 1.2)
mtext("0", side = 2, line = 0.5,  at = 0, cex = 1.5, las = 1)
#legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
legend("bottomright", legend = c("2099", "4899", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 2, lwd = 3, bg = "white", border = "white")
box()

#par(mfrow = c(1, 1), mar = c(5.1, 4.1, 2.1, 0.6))
#ylims <- c(-0.04, 0.035)
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), 
     xlim = c(1650, 1450), ylim = c(-yaxlim, yaxlim - 0.03),
     type = "n", cex.lab = 1.5, axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(h = 0, lty = 1, col = "black")
abline(h = 0)
abline(v = 1605, col = "black", lty = 1, lwd = 2)
mtext("v(C=C)", side = 3, line = 0.3,  at = 1605, cex = 1.2)
abline(v = 1588, col = "black", lty = 2, lwd = 2)
mtext("v(C=C)", side = 3, line = 1.9,  at = 1588, cex = 1.2)
abline(v = 1515, col = "black", lty = 1, lwd = 2)
mtext("v(C=C)", side = 3, line = 0.3,  at = 1515, cex = 1.2)
#lutzke_wavenumbers(cex = 0.8, only = c(1605, 1588, 1515))
rect(diffs_wns$from, 1, diffs_wns$to, -1, col = "#00000055", border = NA)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2), col = plot_cols[1], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[2], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), col = plot_cols[3], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), col = plot_cols[4], lwd = 3)
axis(1, cex.axis = 2, tck = 0.01)
mtext(side = 1, line = 4, at = 1550, cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
#mtext(side = 2, line = 2, at = min(-ylims) + ((max(-ylims) - min(-ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
#legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), col = plot_cols, lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")

box()


plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), 
     xlim = c(1300, 1100), ylim = c(-yaxlim, yaxlim - 0.03),
     type = "n", cex.lab = 1.5, axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
abline(v = 1202, col = "black", lty = 2, lwd = 2)
mtext("\U03B4(CH)", side = 3, line = 0.3,  at = 1202, cex = 1.2)
abline(v = 1168, col = "black", lty = 1, lwd = 2)
mtext("\U03B4(CH)", side = 3, line = 0.3,  at = 1168, cex = 1.2)
abline(h = 0)
#lutzke_wavenumbers(cex = 0.8, only = c(1202, 1168))
rect(diffs_wns$from, 1, diffs_wns$to, -1, col = "#00000055", border = NA)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2), col = plot_cols[1], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[2], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), col = plot_cols[3], lwd = 3)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), col = plot_cols[4], lwd = 3)
axis(1, cex.axis = 2, tck = 0.01)

#mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
#mtext(side = 2, line = 2, at = min(-ylims) + ((max(-ylims) - min(-ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(1125, 0.05, "Plot ref: 0Y9G", cex = 1.5)
box()

#####

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6), pty = "m")
ylims <- c(-0.04, 0.035)
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 950), ylim = -ylims,
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
#abline(h = 0, lty = 1, col = "black")
lutzke_wavenumbers(cex = 0.8)
rect(diffs_wns$from + 2, 1, diffs_wns$to + 2, -1, col = "#00000055", border = NA)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[1], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "2099"), ]), differences = 2), col = plot_cols[2], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2), col = plot_cols[3], lwd = 2)
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), col = plot_cols[4], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(-ylims) + ((max(-ylims) - min(-ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(1790, -0.03, "Plot ref: T5Y0")
box()


spec_MYB99

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6))
ylims <- rev(c(-0.014, 0.027))
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 950), ylim = ylims,
     type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lwd = 2, col = plot_cols[4])
#lutzke_wavenumbers(cex = 0.8, only = c(1605, 1515, 1379, 1168))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(spec_MYB99, differences = 2), col = "black", lwd = 2)
sd_diff_MYB99 <- sd(diff(colMeans(spectra[which(metadata$sample == "9999SR"), 1:440]), differences = 2) - diff(spec_MYB99[1:440], differences = 2))
abline(h = 2 * sd_diff_MYB99, col = "#00000088", lty = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.4, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(1800, sd_diff_MYB99 + 0.002, expression(paste("2", sigma)), col = "#00000088", cex = 1.5)
text(1790, -0.01, "Plot ref: EZH4")
text(c(1605, 1515, 1379, 1168), 
     c(0.013, 0.025, 0.01, 0.018), 
     c(expression(paste(nu,"(C=C)")), 
       expression(paste(nu,"(C=C)")), 
       expression(paste(delta[sym]*"(CH"[3]*")")), 
       expression(paste(delta, "(CH)"))), 
     cex = 1.5)
legend("bottomright", legend = c("MYB99-", "MYB99/MYB99+ mean"), col = c(plot_cols[4], "black"), lty = 1, cex = 1.2, lwd = 3, bg = "white", border = "white")
box()
mtext(side = 2, line = 2, at = -0.015, las = 2, text = "a)", cex = 2)


which(diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(spec_MYB99, differences = 2) > (2 * sd_diff_MYB99))

diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2) - diff(spec_MYB99, differences = 2)

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6))
#ylims <- rev(c(-0.014, 0.027))
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(spec_MYB99_plus, differences = 2), 
     xlim = c(1800, 950), #ylim = range,
     type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lwd = 2, col = plot_cols[3])
#lutzke_wavenumbers(cex = 0.8, only = c(1605, 1515, 1379, 1168))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(spec_MYB99_plus, differences = 2), col = "black", lwd = 2)
sd_diff_MYB99 <- sd(diff(colMeans(spectra[which(metadata$sample == "9999SR"), 1:440]), differences = 2) - diff(spec_MYB99[1:440], differences = 2))
abline(h = 2 * sd(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(spec_MYB99_plus, differences = 2)[1:440]), col = "#00000088", lty = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
#text(1800, sd_diff_MYB99 + 0.002, expression(paste("2", sigma)), col = "#00000088", cex = 1.5)
#text(1790, -0.01, "Plot ref: EZH4")
# text(c(1605, 1515, 1379, 1168), 
#      c(0.013, 0.025, 0.01, 0.018), 
#      c("v(C=C)", "v(C=C)", expression(paste(delta[sym]*"(CH"[3]*")")), expression(paste(delta, "(CH)"))), cex = 1.5)
legend("bottomright", legend = c("MYB99-", "MYB99/MYB99+ mean"), col = c(plot_cols[4], "black"), lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
box()

find_peaks((diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(spec_MYB99_plus, differences = 2))[1:440],
           type = "t")

myb_minus_col_peaks <- find_peaks((diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(spec_MYB99_plus, differences = 2))[1:440],
                                  type = "t")


abline(v = myb_minus_col_peaks[myb_minus_col_peaks$abs < -0.02, ]$peak_wavenumber, col = "red")

which((diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(spec_MYB99_plus, differences = 2)[1:440]) > 2 * sd(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(spec_MYB99_plus, differences = 2)[1:440]))

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 0.6))
#ylims <- rev(c(-0.014, 0.027))
plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), 
     xlim = c(1800, 950), #ylim = range,
     type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lwd = 2, col = plot_cols[3])
#lutzke_wavenumbers(cex = 0.8, only = c(1605, 1515, 1379, 1168))
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2) - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2), col = "black", lwd = 2)
sd_diff_COL_9999SR <- -sd(diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2)[1:440])
abline(h = 2 * sd_diff_COL_9999SR, col = "#00000088", lty = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(ylims) + ((max(ylims) - min(ylims)) / 2), cex = 1.5, text = expression(paste("Relative absorbance (2"^"nd"*" derivative)")))
text(1800, sd_diff_COL_9999SR + 0.002, expression(paste("2", sigma)), col = "#00000088", cex = 1.5)
#text(1790, -0.01, "Plot ref: EZH4")
text(c(1605, 1515, 1379, 1168), 
     c(0.013, 0.025, 0.01, 0.018), 
     c("v(C=C)", "v(C=C)", expression(paste(delta[sym]*"(CH"[3]*")")), expression(paste(delta, "(CH)"))), cex = 1.5)
#legend("bottomright", legend = c("MYB99-", "MYB99/MYB99+ mean"), col = c(plot_cols[4], "black"), lty = 1, cex = 1.3, lwd = 3, bg = "white", border = "white")
box()

which((diff(colMeans(spectra[which(metadata$sample == "COL"), ]), differences = 2)[1:440] - diff(colMeans(spectra[which(metadata$sample == "9999SR"), ]), differences = 2)[1:440]) < 2 * sd_diff_COL_9999SR)



spec_4899 <- colMeans(spectra[which(metadata$sample == "4899"), ])
spec_2099 <- colMeans(spectra[which(metadata$sample == "2099"), ])
spec_COL <- colMeans(spectra[which(metadata$sample == "COL"), ])
spec_MYB99 <- rbind(spec_4899, spec_2099, spec_COL)
spec_MYB99 <- colMeans(spec_MYB99)
spec_MYB99_plus <- rbind(spec_4899, spec_2099)
spec_MYB99_plus <- colMeans(spec_MYB99_plus)
spec_MYB99_norm <- rbind(spec_MYB99_plus, spec_COL)
spec_MYB99_norm <- colMeans(spec_MYB99_norm)

rect(diffs_wns$from, 1, diffs_wns$to, -1, col = "#00000033", border = NA)



plot(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), xlim = c(1800, 948), 
     type = "n", axes = F, xlab = "", ylab = "")
#rect(seek_table2$from, -10, seek_table2$to, 10, col = "#00000016", border = NA)
#rect(seek_table3$from, -10, seek_table3$to, 10, col = "#00000016", border = NA)
abline(h = 0, lty = 1, lwd = 2, col = "black")
lines(wavenumbers[-c(1, length(wavenumbers))], diff(colMeans(spectra[which(metadata$sample == "4899"), ]), differences = 2), col = plot_cols[1], lwd = 2)
lutzke_wavenumbers()
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, at = 0, label = "0", las = 1)
box()


spec_diff <- colMeans(spectra[which(metadata$sample == "2099"), ]) - colMeans(spectra[which(metadata$sample == "COL"), ])
plot(wavenumbers, spec_diff, xlim = c(1800, 948), type = "n", axes = F, xlab = "", ylab = "")
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers, spec_diff, lwd = 2, col = "lightgrey")
points(wavenumbers[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
abline(h = 0, lty = 2, col = "black")
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, at = 0, label = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Mean difference")
box()
rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
rect(900, 0, 940, -10, col = plot_cols[3], border = NA)
text(1050, max(spec_diff), "Plot ref: OQZU")


plot(wavenumbers, spectra[1, ], xlim = rev(range(wavenumbers)), type = "n", ylim = range(spectra), axes = F, xlab = "", ylab = "")
for(i in 1:9){
  lines(wavenumbers, spectra[i, ], col = plot_cols[1])  
}
for(i in 10:17){
  lines(wavenumbers, spectra[i, ], col = plot_cols[2])  
}
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 2, line = 2, at = min(spectra) + ((max(spectra) - min(spectra)) / 2), cex = 1.5, text = "Relative absorbance")
box()

plot(wavenumbers, spec_diff, type = "n", xlim = rev(range(wavenumbers)), axes = F, xlab = "", ylab = "")
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers, spec_diff, lwd = 2, col = "lightgrey")
points(wavenumbers[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
abline(h = 0, lty = 2, col = "black")
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Mean difference")
box()
rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
rect(900, 0, 940, -10, col = plot_cols[2], border = NA)

metadata$treatment <- NA
metadata$treatment[metadata$sample == "9999SR"] <- "MYB99-"
metadata$treatment[metadata$sample == "COL"] <- "CONTROL"
metadata$treatment[is.na(metadata$treatment)] <- "MYB99+"

pc1 <- 1
pc2 <- 2
apca <- prcomp(spectra[, 1:440], center = T, scale. = F)
par(mfrow = c(1, 1), mar = c(5.6, 4.6, 1.1, 1.1))
plot(apca$x[, pc1], apca$x[, pc2], type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- apca$x[, 1:2] %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids <- centroids[c(2, 1, 4, 3), ]
points(PC2 ~ PC1, centroids, pch = 3, cex = 3, col = plot_cols, lwd = 2)
points(apca$x[metadata$sample == "4899", pc1], apca$x[metadata$sample == "4899", pc2], col = plot_cols[1], pch = 17)
points(apca$x[metadata$sample == "2099", pc1], apca$x[metadata$sample == "2099", pc2], col = plot_cols[2], pch = 17)
points(apca$x[metadata$sample == "COL", pc1], apca$x[metadata$sample == "COL", pc2], col = plot_cols[3], pch = 16)
points(apca$x[metadata$sample == "9999SR", pc1], apca$x[metadata$sample == "9999SR", pc2], col = plot_cols[4], pch = 15)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = min(apca$x[, pc1]) + ((max(apca$x[, pc1]) - min(apca$x[, pc1])) / 2), cex = 1.5, text = "PC1 (59.0%)")
mtext(side = 2, line = 3, at = min(apca$x[, pc2]) + ((max(apca$x[, pc2]) - min(apca$x[, pc2])) / 2), cex = 1.5, text = "PC2 (26.6%)")
box()
legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), pch = c(17, 17, 16, 15), col = plot_cols, bty = "n", cex = 1.5)
#text(apca$x[, pc1], apca$x[, pc2], metadata$sample)
text(3.2, 1.5, "Plot ref: YNMP")

# plot(apca$x[, pc1], apca$x[, pc2], type = "n", axes = F, xlab = "", ylab = "")
# abline(h = 0, lty = 2, col = "lightgrey")
# abline(v = 0, lty = 2, col = "lightgrey")
# arrows(0, 0, apca$rotation[, pc1] * 15, apca$rotation[, pc2] * 15)
# points(apca$x[metadata$sample == "2099", pc1], apca$x[metadata$sample == "2099", pc2], col = plot_cols[2], pch = 17)
# points(apca$x[metadata$sample == "COL", pc1], apca$x[metadata$sample == "COL", pc2], col = plot_cols[3], pch = 16)
# points(apca$x[metadata$sample == "9999SR", pc1], apca$x[metadata$sample == "9999SR", pc2], col = plot_cols[4], pch = 15)
# axis(1, cex.axis = 1.5, tck = 0.01)
# axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
# mtext(side = 1, line = 4, at = min(apca$x[, pc1]) + ((max(apca$x[, pc1]) - min(apca$x[, pc1])) / 2), cex = 1.5, text = "PC1 (59.0%)")
# mtext(side = 2, line = 3, at = min(apca$x[, pc2]) + ((max(apca$x[, pc2]) - min(apca$x[, pc2])) / 2), cex = 1.5, text = "PC2 (26.6%)")
# box()
# legend("bottomright", legend = c("4899", "2099", "COL", "9999SR"), pch = c(17, 17, 16, 15), col = plot_cols, bty = "n", cex = 1.5)
# #text(apca$x[, pc1], apca$x[, pc2], metadata$sample)
# text(3.2, 1.5, "Plot ref: YNMP")


par(mfrow = c(2, 1), mar = c(5.1, 6.1, 4.1, 0.6), pty = "m")
plot(wavenumbers[1:440], apca$rotation[, 1], xlim = rev(range(wavenumbers[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[1:440], apca$rotation[, 1], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
#mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca$rotation[, pc1]) + ((max(apca$rotation[, pc1]) - min(apca$rotation[, pc1])) / 2), cex = 1.5, text = "PC1 loadings")
box()
lutzke_wavenumbers()

#par(mfrow = c(1, 1), mar = c(5.1, 6.1, 4.1, 0.1))
plot(wavenumbers[1:440], apca$rotation[, 2], xlim = rev(range(wavenumbers[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[1:440], apca$rotation[, 2], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca$rotation[, 2]) + ((max(apca$rotation[, 2]) - min(apca$rotation[, 2])) / 2), cex = 1.5, text = "PC2 loadings")
box()
lutzke_wavenumbers()
mtext(side = 1, line = 3, at = 960, text = "Plot ref: 5512")

pc1 <- 1
pc2 <- 2
apca <- prcomp(spectra[, 1:440], center = T, scale. = F)
par(mfrow = c(1, 1), mar = c(5.6, 4.6, 1.1, 1.1))
plot(apca$x[, pc1], apca$x[, pc2], type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- apca$x[, 1:2] %>% bind_cols(treatment = metadata$treatment) %>% group_by(treatment) %>% summarise_all("mean") %>% data.frame()
#centroids <- centroids[c(2, 1, 4, 3), ]
points(PC2 ~ PC1, centroids, pch = 3, cex = 3, col = plot_cols[c(3, 4, 1)], lwd = 2)
points(apca$x[metadata$treatment == "MYB99+", pc1], apca$x[metadata$treatment == "MYB99+", pc2], col = plot_cols[1], pch = 16)
#points(apca$x[metadata$sample == "2099", pc1], apca$x[metadata$sample == "2099", pc2], col = plot_cols[2], pch = 17)
points(apca$x[metadata$treatment == "CONTROL", pc1], apca$x[metadata$treatment == "CONTROL", pc2], col = plot_cols[3], pch = 16)
points(apca$x[metadata$treatment == "MYB99-", pc1], apca$x[metadata$treatment == "MYB99-", pc2], col = plot_cols[4], pch = 16)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = min(apca$x[, pc1]) + ((max(apca$x[, pc1]) - min(apca$x[, pc1])) / 2), cex = 1.5, text = "PC1 (47.1%)")
mtext(side = 2, line = 3, at = min(apca$x[, pc2]) + ((max(apca$x[, pc2]) - min(apca$x[, pc2])) / 2), cex = 1.5, text = "PC2 (30.5%)")
box()
legend("bottomleft", legend = c("MYB99-", "CONTROL", "MYB99+"), pch = 16, col = plot_cols[c(4, 3, 1)], bty = "n", cex = 1.5)
#text(apca$x[, pc1], apca$x[, pc2], metadata$sample)
text(24, -17, "Plot ref: 4JXR")


pc1 <- 3
pc2 <- 4
apca <- prcomp(spectra[, 1:440], center = T, scale. = F)
par(mfrow = c(1, 1), mar = c(5.6, 4.6, 1.1, 1.1))
plot(apca$x[, pc1], apca$x[, pc2], type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- apca$x[, 1:2] %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
#centroids <- centroids[c(2, 1, 4, 3), ]
points(PC2 ~ PC1, centroids, pch = 3, cex = 3, col = plot_cols[c(3, 4, 1)], lwd = 2)
points(apca$x[metadata$sample == "4899", pc1], apca$x[metadata$sample == "4899", pc2], col = plot_cols[4], pch = 16)
points(apca$x[metadata$sample == "2099", pc1], apca$x[metadata$sample == "2099", pc2], col = plot_cols[3], pch = 16)
points(apca$x[metadata$sample == "COL", pc1], apca$x[metadata$sample == "COL", pc2], col = plot_cols[2], pch = 15)
points(apca$x[metadata$sample == "9999SR", pc1], apca$x[metadata$sample == "9999SR", pc2], col = plot_cols[1], pch = 17)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = min(apca$x[, pc1]) + ((max(apca$x[, pc1]) - min(apca$x[, pc1])) / 2), cex = 1.5, text = "PC1 (47.1%)")
mtext(side = 2, line = 3, at = min(apca$x[, pc2]) + ((max(apca$x[, pc2]) - min(apca$x[, pc2])) / 2), cex = 1.5, text = "PC2 (30.5%)")
box()
legend("bottomleft", legend = c("MYB99-", "CONTROL", "MYB99+"), pch = 16, col = plot_cols[c(4, 3, 1)], bty = "n", cex = 1.5)
#text(apca$x[, pc1], apca$x[, pc2], metadata$sample)
text(24, -17, "Plot ref: XXXX")


################ derivative

spectra_d <- t(apply(spectra, 1, diff, differences = 2)) %>% as.data.frame()
wavenumbers_d <- wavenumbers[2:(length(wavenumbers) - 1)]


# par(mfrow = c(2, 1), mar = c(5.1, 4.1, 0.1, 0.1))
# plot(wavenumbers_d, spectra_d[1, ], xlim = c(1800, 948), type = "n", ylim = range(spectra_d), axes = F, xlab = "", ylab = "")
input.var.equal <- F
t_test_results <- numeric(ncol(spectra_d))
for(i in 1:ncol(spectra_d)){
  t_test <- t.test(spectra_d[1:9, i], spectra_d[10:17, i], var.equal = input.var.equal)
  t_test_results[i] <- t_test$p.value
}

t_test_results[which(t_test_results > 0.05)] <- NA

seek_table <- seek_seq(which(t_test_results <= 0.05))
seek_table$from <- wavenumbers_d[seek_table$from]
seek_table$to <- wavenumbers_d[seek_table$to]
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)

# for(i in 1:9){
#   lines(wavenumbers_d, spectra_d[i, ], col = plot_cols[1])  
# }
# for(i in 10:17){
#   lines(wavenumbers_d, spectra_d[i, ], col = plot_cols[2])  
# }
# axis(1, cex.axis = 1.5, tck = 0.01)
# 
# mtext(side = 2, line = 2, at = min(spectra_d) + ((max(spectra_d) - min(spectra_d)) / 2), cex = 1.5, text = "Relative absorbance")
# #internal_axis_ticks(1, 1)
# legend("bottomright", legend = unique(metadata$sample), lty = 1, col = plot_cols, bty = "n", lwd = 2, cex = 1.5)
# box()

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 2.1, 0.1))
spec_diff <- colMeans(spectra_d[1:9, ]) - colMeans(spectra_d[10:17, ])
plot(wavenumbers_d, spec_diff, xlim = c(1800, 948), type = "n", axes = F, xlab = "", ylab = "")
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers_d, spec_diff, lwd = 2, col = "lightgrey")
points(wavenumbers_d[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
abline(h = 0, lty = 2, col = "black")
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, at = 0, label = "0", las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Mean difference")
box()
rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
rect(900, 0, 940, -10, col = plot_cols[2], border = NA)
phil_wavenumbers()
text(1770, max(spec_diff), "Plot ref: EL7Z")



# plot(wavenumbers_d, spectra_d[1, ], xlim = rev(range(wavenumbers_d)), type = "n", ylim = range(spectra_d), axes = F, xlab = "", ylab = "")
# for(i in 1:9){
#   lines(wavenumbers_d, spectra_d[i, ], col = plot_cols[1])  
# }
# for(i in 10:17){
#   lines(wavenumbers_d, spectra_d[i, ], col = plot_cols[2])  
# }
# rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
# axis(1, cex.axis = 1.5, tck = 0.01)
# mtext(side = 2, line = 2, at = min(spectra_d) + ((max(spectra_d) - min(spectra_d)) / 2), cex = 1.5, text = "Relative absorbance")
# box()
# 
# plot(wavenumbers_d, spec_diff, type = "n", xlim = rev(range(wavenumbers_d)), axes = F, xlab = "", ylab = "")
# rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
# lines(wavenumbers_d, spec_diff, lwd = 2, col = "lightgrey")
# points(wavenumbers_d[which(t_test_results < 0.05)], spec_diff[which(t_test_results < 0.05)], lwd = 2, col = "black", pch = ".", cex = 3)
# abline(h = 0, lty = 2, col = "black")
# axis(1, cex.axis = 1.5, tck = 0.01)
# mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
# mtext(side = 2, line = 2, at = min(spec_diff) + ((max(spec_diff) - min(spec_diff)) / 2), cex = 1.5, text = "Mean difference")
# box()
# rect(900, 0, 940, 10, col = plot_cols[1], border = NA)
# rect(900, 0, 940, -10, col = plot_cols[2], border = NA)



pc1 <- 2
pc2 <- 3
apca <- prcomp(spectra_d[, 1:440], center = T, scale. = T)
par(mfrow = c(1, 1))
plot(apca$x[, pc1], apca$x[, pc2], type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
points(apca$x[metadata$sample == "9999SR", pc1], apca$x[metadata$sample == "9999SR", pc2], col = plot_cols[1], pch = 16)
points(apca$x[metadata$sample == "COL", pc1], apca$x[metadata$sample == "COL", pc2], col = plot_cols[2], pch = 16)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = min(apca$x[, pc1]) + ((max(apca$x[, pc1]) - min(apca$x[, pc1])) / 2), cex = 1.5, text = "PC1 (35.0%)")
mtext(side = 2, line = 3, at = min(apca$x[, pc2]) + ((max(apca$x[, pc2]) - min(apca$x[, pc2])) / 2), cex = 1.5, text = "PC2 (20.0%)")
box()
legend("topright", legend = unique(metadata$sample), pch = 16, col = plot_cols, bty = "n", cex = 1.5)
text(10, -8, "Plot ref: 3LHR")

par(mfrow = c(1, 1), mar = c(5.1, 6.1, 2.1, 0.1))
plot(wavenumbers_d[1:440], apca$rotation[, pc2], xlim = rev(range(wavenumbers_d[1:440])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 2, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 2, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 2, col = "#FF000055", lwd = 2)
rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers_d[1:440], apca$rotation[, pc2], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(apca$rotation[, pc2]) + ((max(apca$rotation[, pc2]) - min(apca$rotation[, pc2])) / 2), cex = 1.5, text = "PC2 loadings")
box()
phil_wavenumbers()
text(1770, 0.13, "Plot ref: 43XY")


(detection_sensitivity_compare <- 3)
peak_table_list <- vector("list", length = nrow(spectra))

for(i in 1:nrow(spectra_d)){
  print(i)
  (my_spectrum1 <- spectra_d[i, ])
  
  (xx <- names(my_spectrum1))
  #my_spectrum1 <- as.numeric(my_spectrum1)
  
  (x <- which(my_spectrum1 < detection_sensitivity_compare * -sd(my_spectrum1))) #### compare against one or another
  
  seek_seq <- function(x){
    skipped_indices <- which(diff(x) != 1) #### only works if [1] is not a skip
    to <- x[skipped_indices]
    to <- c(to, x[length(x)])
    from <- c(1, skipped_indices + 1)
    from <- x[from]
    data.frame(from = from, to = to)
  }
  
  (seqs <- seek_seq(x))
  (peak_sequences <- apply(seqs, 1, function(x){x[1]:x[2]}))
  (peak_wns <- lapply(peak_sequences, function(x){
    names(which.min(my_spectrum1[x]))
  }))
  (peak_values <- lapply(peak_sequences, function(x){
    min(my_spectrum1[x])
  }))
  (peaks <- data.frame(wn = unlist(peak_wns), value = -unlist(peak_values)))
  (peaks$ID<- row.names(my_spectrum1))
  (peak_table_list[[i]] <- peaks)
}
peaks_df <- as.data.frame(do.call(rbind, peak_table_list))

peaks_df$value<- as.numeric(peaks_df$value)
#peaks_df$peaks_id<- peaks_df$ID %>% str_replace_all("[0-9]", "") %>% str_replace_all("_", "")
#peaks_df$peaks_id <- factor(peaks_df$peaks_id)
peaks_df$wn <- peaks_df$wn %>% strip_non_numbers()

peaks1 <- peaks_df[as.numeric(peaks_df$ID) <= 9, ]
(nd <- data.frame(table(peaks1$wn)))
(nd$Var1 <- nd$Var1 %>% as.character() %>% as.numeric())
(nd$Freq_per<- (nd$Freq* 100)/ 9)

(peaks2 <- peaks_df[as.numeric(peaks_df$ID) > 9, ])
(nd2 <- data.frame(table(peaks2$wn)))
(nd2$Var1 <- nd2$Var1 %>% as.character() %>% as.numeric())
(nd2$Freq_per <- (nd2$Freq* 100)/ 8)

#nd$quality <- "good"

plot(nd$Var1, nd$Freq_per, type = "n", xlim = rev(range(wavenumbers)), axes = F, xlab = expression(paste("Wavenumber (cm"^-1*")")), ylab = "% of spectra", cex.lab = 1.5, ylim = c(0, 100))
arrows(nd$Var1, 0, nd$Var1, nd$Freq_per, angle = 0, length = 0, col = plot_cols[1], lwd = 3)
arrows(nd2$Var1, 0, nd2$Var1, nd2$Freq_per, angle = 0, length = 0, col = plot_cols[2], lwd = 3)
axis(1, tck = 0.01, cex.axis = 1.3)
axis(2, las = 1, tck = 0.01, cex.axis = 1.3)
box()




#############################################################################################################


par(mfrow = c(1,2))                            # Splitting the plot in 2.
#biplot(apca)                                    # In-built biplot() R func.
PCA <- apca

choices = 1:2                                  # Selecting first two PC's
scale = 1                                      # Default
scores = PCA$x                                  # The scores
lam = PCA$sdev[choices]                        # Sqrt e-vals (lambda) 2 PC's
n = nrow(scores)                               # no. rows scores
lam = lam * sqrt(n)                            # See below.


x = t(t(scores[,choices]) / lam)         # scaled scores
y = t(t(PCA$rotation[,choices]) * lam)         # scaled eigenvecs (loadings)

n = nrow(x)                                    # Same as dataset (150)
p = nrow(y)                                    # Three var -> 3 rows

xlabs = 1L:n
xlabs = as.character(xlabs)                    # no. from 1 to 150 
dimnames(x) = list(xlabs, dimnames(x)[[2L]])   # no's and PC1 / PC2

ylabs = dimnames(y)[[1L]]                      # Iris species
ylabs = as.character(ylabs)
dimnames(y) <- list(ylabs, dimnames(y)[[2L]])  # Species and PC1/PC2


unsigned.range = function(x) c(-abs(min(x, na.rm = TRUE)), 
                               abs(max(x, na.rm = TRUE)))
rangx1 = unsigned.range(x[, 1L])               # Range first col x
rangx2 = unsigned.range(x[, 2L])               # Range second col x
rangy1 = unsigned.range(y[, 1L])               # Range 1st scaled evec
rangy2 = unsigned.range(y[, 2L])               # Range 2nd scaled evec


(xlim = ylim = rangx1 = rangx2 = range(rangx1, rangx2))

(ratio = max(rangy1/rangx1, rangy2/rangx2)) 

par(mfrow = c(1, 2), pty = "s")                                
pc1 <- 1
pc2 <- 2
plot(x, type = "n", 
     #xlim = xlim, 
     #ylim = ylim, 
     xlim = c(-0.32, 0.32), 
     ylim = c(-0.32, 0.32), 
     axes = F, 
     xlab = "", 
     ylab = ""
)  
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- x %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids <- centroids[c(2, 1, 4, 3), ]
points(PC2 ~ PC1, centroids, pch = 3, cex = 3, col = plot_cols, lwd = 3)
points(x[metadata$sample == "4899", pc1], x[metadata$sample == "4899", pc2], col = plot_cols[1], pch = 17, cex = 1.5)
points(x[metadata$sample == "2099", pc1], x[metadata$sample == "2099", pc2], col = plot_cols[2], pch = 17, cex = 1.5)
points(x[metadata$sample == "COL", pc1], x[metadata$sample == "COL", pc2], col = plot_cols[3], pch = 16, cex = 1.5)
points(x[metadata$sample == "9999SR", pc1], x[metadata$sample == "9999SR", pc2], col = plot_cols[4], pch = 15, cex = 1.5)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 3, at = 0, cex = 1.5, text = "PC1 scores (59.0%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 scores (26.6%)")
legend("topleft", legend = c("4899", "2099", "COL", "9999SR"), pch = c(17, 17, 16, 15), col = plot_cols, bty = "n", cex = 1.4)
box()

(new_xlim = xlim * ratio)

(new_ylim = ylim * ratio)
#par(new = T)
#plot(y, axes = F, type = "n", xlim = new_xlim, ylim = new_ylim, xlab = "", ylab = "")
#vegan::ordihull(x * 5, metadata$sample, col = plot_cols[c(2, 1, 4, 3)], lwd = 2)

lutzke_indices <- which(wavenumbers %in% lutzke_wavenumbers_only[lutzke_wavenumbers_only < 1800])

axis_lim <- 1.15
plot(y, axes = FALSE, type = "n", 
     xlim = new_xlim, 
     ylim = new_ylim, 
     xlab = "", 
     ylab = "", 
     #xlim = c(-axis_lim, axis_lim),
     #ylim = c(-axis_lim, axis_lim)
)

abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 3, at = 0, cex = 1.5, text = "PC1 loadings (59.0%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 loadings (26.6%)")

text(y[lutzke_indices, ], labels = round(strip_non_numbers(ylabs[lutzke_indices])), col = 1, cex = 0.8)  # This just prints the species
arrow.len <- 0.1
arrows(0, 0, y[lutzke_indices, 1L] * 0.8, y[lutzke_indices, 2L] * 0.8, length = arrow.len, col = "#00000033")

#non_sporopollenin_indices <- sapply(c(1642, 1022, 1541, 1558, 1455, 1750, 1772, 1080), function(x){which.min(abs(x - wavenumbers))})
#text(y[non_sporopollenin_indices, ], labels = round(strip_non_numbers(ylabs[non_sporopollenin_indices])), cex = 0.8, col = "red")  # This just prints the species
#arrows(0, 0, y[non_sporopollenin_indices, 1L] * 0.8, y[non_sporopollenin_indices, 2L] * 0.8, length = arrow.len, col = "#FF000033")
#legend("topleft", lty = 1, lwd = 2,  col = 1:2, legend = c("Sporopollenin band", "Other band"), cex = 1.4, bty = "n")
text(1.2, 1.5, "Plot ref: QCIV")
box()



############################ 2nd der

par(mfrow = c(1,2))                            # Splitting the plot in 2.
#biplot(apca)                                    # In-built biplot() R func.
PCA <- prcomp(-data.frame(t(diff(t(spectra[, 1:440]), differences = 2))), center = T, scale. = F)

choices = 1:2                                  # Selecting first two PC's
scale = 1                                      # Default
scores = PCA$x                                  # The scores
lam = PCA$sdev[choices]                        # Sqrt e-vals (lambda) 2 PC's
n = nrow(scores)                               # no. rows scores
lam = lam * sqrt(n)                            # See below.


x = t(t(scores[,choices]) / lam)         # scaled scores
y = t(t(PCA$rotation[,choices]) * lam)         # scaled eigenvecs (loadings)

n = nrow(x)                                    # Same as dataset (150)
p = nrow(y)                                    # Three var -> 3 rows

xlabs = 1L:n
xlabs = as.character(xlabs)                    # no. from 1 to 150 
dimnames(x) = list(xlabs, dimnames(x)[[2L]])   # no's and PC1 / PC2

ylabs = dimnames(y)[[1L]]                      # Iris species
ylabs = as.character(ylabs)
dimnames(y) <- list(ylabs, dimnames(y)[[2L]])  # Species and PC1/PC2


unsigned.range = function(x) c(-abs(min(x, na.rm = TRUE)), 
                               abs(max(x, na.rm = TRUE)))
rangx1 = unsigned.range(x[, 1L])               # Range first col x
rangx2 = unsigned.range(x[, 2L])               # Range second col x
rangy1 = unsigned.range(y[, 1L])               # Range 1st scaled evec
rangy2 = unsigned.range(y[, 2L])               # Range 2nd scaled evec


(xlim = ylim = rangx1 = rangx2 = range(rangx1, rangx2))

(ratio = max(rangy1/rangx1, rangy2/rangx2)) 

par(mfrow = c(1, 2), pty = "s")                                
pc1 <- 1
pc2 <- 2
plot(x, type = "n", 
     #xlim = xlim, 
     #ylim = ylim, 
     xlim = c(-0.32, 0.32), 
     ylim = c(-0.32, 0.32), 
     axes = F, 
     xlab = "", 
     ylab = ""
)  
abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
centroids <- x %>% bind_cols(sample = metadata$sample) %>% group_by(sample) %>% summarise_all("mean") %>% data.frame()
centroids <- centroids[c(2, 1, 4, 3), ]
points(PC2 ~ PC1, centroids, pch = 3, cex = 3, col = plot_cols, lwd = 2)
points(x[metadata$sample == "4899", pc1], x[metadata$sample == "4899", pc2], col = plot_cols[1], pch = 17)
points(x[metadata$sample == "2099", pc1], x[metadata$sample == "2099", pc2], col = plot_cols[2], pch = 17)
points(x[metadata$sample == "COL", pc1], x[metadata$sample == "COL", pc2], col = plot_cols[3], pch = 16)
points(x[metadata$sample == "9999SR", pc1], x[metadata$sample == "9999SR", pc2], col = plot_cols[4], pch = 15)
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 0, cex = 1.5, text = "PC1 scores (82.4%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 scores (6.1%)")
legend("topleft", legend = c("4899", "2099", "COL", "9999SR"), pch = c(17, 17, 16, 15), col = plot_cols, bty = "n", cex = 1.5)
box()

(new_xlim = xlim * ratio)

(new_ylim = ylim * ratio)

lutzke_indices <- which(wavenumbers %in% lutzke_wavenumbers_only[lutzke_wavenumbers_only < 1800]) - 2

axis_lim <- 1.15
plot(y, axes = FALSE, type = "n", 
     xlim = new_xlim, 
     ylim = new_ylim, 
     #ylim = c(-0.1, 0.1),
     xlab = "", 
     ylab = "", 
     #xlim = c(-axis_lim, axis_lim),
     #ylim = c(-axis_lim, axis_lim)
)

abline(h = 0, lty = 2, col = "lightgrey")
abline(v = 0, lty = 2, col = "lightgrey")
axis(1, cex.axis = 1.5, tck = 0.01)
mtext(side = 1, line = 4, at = 0, cex = 1.5, text = "PC1 loadings (82.4%)")
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 2, line = 4, at = 0, cex = 1.5, text = "PC2 loadings (6.1%)")

text(y[lutzke_indices, ], labels = round(strip_non_numbers(ylabs[lutzke_indices])), col = 1, cex = 0.8)  # This just prints the species
arrow.len <- 0.1
arrows(0, 0, y[lutzke_indices, 1L] * 0.8, y[lutzke_indices, 2L] * 0.8, length = arrow.len, col = "#00000033")

non_sporopollenin_indices <- sapply(c(1642, 1022, 1541, 1558, 1455, 1750, 1772, 1080), function(x){which.min(abs(x - wavenumbers))})
text(y[non_sporopollenin_indices, ], labels = round(strip_non_numbers(ylabs[non_sporopollenin_indices])), cex = 0.8, col = "red")  # This just prints the species
arrows(0, 0, y[non_sporopollenin_indices, 1L] * 0.8, y[non_sporopollenin_indices, 2L] * 0.8, length = arrow.len, col = "#FF000033")
legend("topleft", lty = 1, col = 1:2, legend = c("Sporopollenin", "Other (amide, etc.)"), cex = 1.5, bty = "n")
text(0.2, 0.24, "Plot ref: FSG8")
box()

par(mfrow = c(2, 1), mar = c(5.1, 6.1, 4.1, 0.6), pty = "m")
plot(wavenumbers[2:439], PCA$rotation[, 1], xlim = rev(range(wavenumbers[2:439])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[2:439], PCA$rotation[, 1], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
#mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(PCA$rotation[, pc1]) + ((max(PCA$rotation[, pc1]) - min(PCA$rotation[, pc1])) / 2), cex = 1.5, text = "PC1 loadings")
box()
lutzke_wavenumbers()

#par(mfrow = c(1, 1), mar = c(5.1, 6.1, 4.1, 0.1))
plot(wavenumbers[2:439], PCA$rotation[, 2], xlim = rev(range(wavenumbers[2:439])), type = "n", axes = F, xlab = "", ylab = "")
abline(h = 0, lty = 1, col = "lightgrey", lwd = 2)
abline(h = sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
abline(h = -sqrt(1/440), lty = 1, col = "#FF000055", lwd = 2)
#rect(seek_table$from, -10, seek_table$to, 10, col = "#00000016", border = NA)
lines(wavenumbers[2:439], PCA$rotation[, 2], lwd = 2)
axis(1, cex.axis = 1.5, tck = 0.01)
axis(2, cex.axis = 1.5, tck = 0.01, las = 1)
mtext(side = 1, line = 4, at = 1800 - ((1800 - 948) / 2), cex = 1.5, text = expression(paste("Wavenumber (cm"^-1*")")))
mtext(side = 2, line = 4.5, at = min(PCA$rotation[, 2]) + ((max(PCA$rotation[, 2]) - min(PCA$rotation[, 2])) / 2), cex = 1.5, text = "PC2 loadings")
box()
lutzke_wavenumbers()
mtext(side = 1, line = 3, at = 960, text = "Plot ref: K7SO")


peaks <- find_peaks(spec_diff)
peaks2 <- find_peaks(spec_diff2)
peaks3 <- find_peaks(spec_diff3)

intersecting_wavenumbers1 <- intersect(peaks$peak_wavenumber, peaks2$peak_wavenumber)
intersecting_wavenumbers2 <- intersect(peaks$peak_wavenumber, peaks3$peak_wavenumber)
#intersecting_wavenumbers <- intersect(intersecting_wavenumbers1, intersecting_wavenumbers2) %>% sort()
#intersecting_wavenumbers <- unique(intersecting_wavenumbers1, intersecting_wavenumbers2) %>% sort()
intersecting_wavenumbers <- unique(c(peaks$peak_wavenumber, peaks2$peak_wavenumber, peaks3$peak_wavenumber)) %>% sort()

opposite_wavenumbers <- wavenumbers[which(spec_diff < 0 & spec_diff2 > 0 & spec_diff3 > 0)]
opposite_wavenumbers <- sort(unique(c(opposite_wavenumbers, wavenumbers[which(spec_diff > 0 & spec_diff2 < 0 & spec_diff3 < 0)])))
intersecting_wavenumbers <- intersect(intersecting_wavenumbers, opposite_wavenumbers)

par(mfrow = c(1, 1))
plot(wavenumbers, spec_diff, xlim = c(1800, 948), type = "l")
points(wavenumbers, spec_diff)
#abline(v = peaks[peaks$peak_type == "peak", ]$peak_wavenumber, col = "red")
#abline(v = peaks[peaks$peak_type == "trough", ]$peak_wavenumber, col = "blue")
abline(h = 0)
lines(wavenumbers, spec_diff2)
lines(wavenumbers, spec_diff3)
#abline(v = intersecting_wavenumbers2, col = "red")
#abline(v = intersecting_wavenumbers, col = "red")

interesting_wavenumbers <- c(
  1024.020,
  1078.010,
  1168.650,
  1232.290,
  1261.220,
  1290.140,
  1436.710,
  1513.85,
  1604.480,
  1616.060,
  1710.550
)

lutzke_wavenumbers_only <- c(1605, 1588, 1514, 3026, 1202, 1168, 833, 1680, 1630, 985, 2926, 2855, 1463, 720, 1411, 1710, 1437, 1283, 943, 1125, 1103, 3350, 1323)
lutzke_wavenumbers_only <- lutzke_wavenumbers_only[lutzke_wavenumbers_only > 950]
lutzke_wavenumbers_only <- sapply(lutzke_wavenumbers_only, function(x){wavenumbers[which.min(abs(x - wavenumbers))]}) %>% sort()

abline(v = interesting_wavenumbers)

interesting_indices <- which(wavenumbers %in% interesting_wavenumbers)


###### using lutzke spectra
ferulic <- read.csv("C:\\Users\\sbzmk2\\Documents\\Post-doc\\Data from other sources\\Lutzke et al. Pinus ponderosa\\ATR-FTIR\\R-5 References\\R-5.11.csv", header = F)
head(ferulic)
names(ferulic) <- c("wavenumber", "abs")
#plot(abs ~ wavenumber, ferulic, xlim = rev(range(wavenumber)), type = "l")
plot(abs ~ wavenumber, ferulic, xlim = c(1800, 700), ylim = c(50, 100), type = "l")
abline(v = 1411, col = "blue")
abline(v = 1379, col = "red")

coumaric <- read.csv("C:\\Users\\sbzmk2\\Documents\\Post-doc\\Data from other sources\\Lutzke et al. Pinus ponderosa\\ATR-FTIR\\R-5 References\\R-5.05.csv", header = F)
names(coumaric) <- c("wavenumber", "abs")
plot(abs ~ wavenumber, coumaric, xlim = c(1800, 700), ylim = c(0, 100), type = "l")
abline(v = 1411, col = "blue")
abline(v = 1379, col = "red")

ext <- 5
ferulic <- ferulic[!ferulic$wavenumber < 900, ]
coumaric <- coumaric[!coumaric$wavenumber < 900, ]
plot(abs ~ wavenumber, ferulic, xlim = c(1800, 950), ylim = c(10, 100), type = "n", las = 1, axes = F, cex.lab = 1.5, xlab = lab("wn"), ylab = "")
#abline(v = 1411, col = "blue")
#abline(v = 1379, col = "red")
rect(1411 - ext, -100, 1411 + ext, 1000, border = NA, col = "#FF000022")
rect(1379 - ext, -100, 1379 + ext, 1000, border = NA, col = "#00000022")
lines(abs ~ wavenumber, ferulic, lwd = 2)
lines(abs - 25 ~ wavenumber, coumaric, xlim = c(1800, 950), ylim = c(0, 100), type = "l", lwd = 2)
internal_axis_ticks(1, 1, cex.axis = 1.5)
text(1350, 97, "trans-4-hydroxy-3-methoxycinnamic acid", adj = 0, cex = 1.5)
text(1350, 15, "trans-4-hydroxycinnamic acid", adj = 0, cex = 1.5)
mtext(side = 2, line = 2, at = midpoint(c(10, 100)), text = "Transmittance", cex = 1.5)
text(1750, 15, "Plot ref: KS91")

methyl_hc <- read.csv("C:\\Users\\sbzmk2\\Documents\\Post-doc\\Data from other sources\\Lutzke et al. Pinus ponderosa\\ATR-FTIR\\R-5 References\\R-5.16.csv", header = F)
names(methyl_hc) <- c("wavenumber", "abs")
plot(abs ~ wavenumber, methyl_hc, xlim = c(1800, 700), ylim = c(0, 100), type = "l")
abline(v = 1411, col = "blue")
abline(v = 1379, col = "red")

methoxy <- read.csv("C:\\Users\\sbzmk2\\Documents\\Post-doc\\Data from other sources\\Lutzke et al. Pinus ponderosa\\ATR-FTIR\\R-5 References\\R-5.15.csv", header = F)
names(methoxy) <- c("wavenumber", "abs")
plot(abs ~ wavenumber, methoxy, xlim = c(1800, 700), ylim = c(0, 100), type = "l")
abline(v = 1411, col = "blue")
abline(v = 1379, col = "red")

ext <- 5
ferulic <- ferulic[!ferulic$wavenumber < 900, ]
coumaric <- coumaric[!coumaric$wavenumber < 900, ]
methyl_hc <- methyl_hc[!methyl_hc$wavenumber < 900, ]
methoxy <- methoxy[!methoxy$wavenumber < 900, ]
plot(abs ~ wavenumber, ferulic, 
     xlim = c(1800, 950), 
     ylim = c(-100, 110),
     type = "n", las = 1, axes = F, cex.lab = 1.5, xlab = lab("wn"), ylab = "")
#abline(v = 1411, col = "blue")
#abline(v = 1379, col = "red")
rect(1411 - ext, -1000, 1411 + ext, 1000, border = NA, col = "#FF000022")
rect(1379 - ext, -1000, 1379 + ext, 1000, border = NA, col = "#00000022")
lines(abs + 2 ~ wavenumber, ferulic, lwd = 2)
lines(abs - 38 ~ wavenumber, coumaric, xlim = c(1800, 950), ylim = c(0, 100), type = "l", lwd = 2)
lines(abs - 93 ~ wavenumber, methyl_hc, xlim = c(1800, 950), ylim = c(0, 100), type = "l", lwd = 2)
lines(abs - 140 ~ wavenumber, methoxy, xlim = c(1800, 950), ylim = c(0, 100), type = "l", lwd = 2)
internal_axis_ticks(1, 1, cex.axis = 1.5)
text(1350, 100, "trans-4-hydroxy-3-methoxycinnamic acid", adj = 0, cex = 1.5)
text(1350, 54, "trans-4-hydroxycinnamic acid", adj = 0, cex = 1.5)
text(1350, -3, "trans-4-methoxycinnamic acid", adj = 0, cex = 1.5)
text(1350, -50, "methyl trans-4-hydroxycinnamate", adj = 0, cex = 1.5)
mtext(side = 2, line = 2, at = midpoint(c(-75, 100)), text = "Transmittance", cex = 1.5)
text(1750, -95, "Plot ref: KS91")

