
# let's start with spectra pre-processed (scatter corrected and scaled) by my app "specsplorer"
# which is currently set as a private repository on GitHub whilst it's in development
md <- readRDS("Pinaceae_Cupressaceae_metadata.RDS")
sp <- readRDS("Pinaceae_Cupressaceae_spectra.RDS")
md <- md$metadata
sp <- sp$spectra

# this is the head of the metadata. there's lots of scan information, but what is relevant to
# what we'll be plotting is the taxonomic information found in columns such as "family" and
# "genus"
head(md)

# and these are the spectra. the spectra are currently in a wide format with a sample identifier
# in the "X" column and wavenumbers as other columns. observations are in rows
spectra[1:10, 1:10]

# there are many wavenumbers at which measurements were made
ncol(spectra) - 1

# my first handy function is get_mean_spectra, which calculates the mean spectrum according to
# a grouping vector which is currently in the metadata. this information gets incorporated
# into the output tidy spectra file which is immediately ready for plotting, containing also
# data for plotting standard deviation (sd) upper and lower bounds

# there are two taxonomic families, so we are going to plot the data that way, and so we will
# also choose two colours for the plotting
sp_tidy <- get_mean_spectra(sp, md$family)
col_vec <- scico::scico(length(unique(md$family)), palette = "roma")

# and to plot those data...
sp_tidy %>% ggplot(aes(wavenumber, mean)) + 
  geom_ribbon(aes(wavenumber, ymin = lower, ymax = upper, fill = group), alpha = 0.3) +
  geom_line(aes(col = group)) + 
  scale_x_reverse() + 
  theme_bw() +
  theme(legend.position = "top",
        axis.title = element_text(size = 20)) +
  scale_color_manual(values = col_vec,
                     name = "Taxon") +
  scale_fill_manual(values = col_vec,
                    name = "Taxon") +
  labs(x = expression(paste("Wavenumber (cm"^"-1"*")")),
       y = "Relative absorbance")

# get_mean_spectra also can quickly get a smoothed derivative. it defaults to second derivative,
# but this can be controlled. Now plotting without the sd ribbons and a smoothing window size
# of 7:
sp_tidy <- get_mean_spectra(sp, md$family, smoothed_derivative = T, w = 7)

sp_tidy %>% ggplot(aes(wavenumber, mean)) + 
  #geom_ribbon(aes(wavenumber, ymin = lower, ymax = upper, fill = group), alpha = 0.3) +
  geom_line(aes(col = group)) + 
  scale_x_reverse() + 
  theme_bw() +
  theme(legend.position = "top",
        axis.title = element_text(size = 20)) +
  scale_color_manual(values = col_vec,
                     name = "Taxon") +
  scale_fill_manual(values = col_vec,
                    name = "Taxon") +
  labs(x = expression(paste("Wavenumber (cm"^"-1"*")")),
       y = "Relative absorbance")

# and now plotting the first derivative (same w) without the sd ribbons:
sp_tidy <- get_mean_spectra(sp, md$family, smoothed_derivative = T, w = 7, m = 1)

sp_tidy %>% ggplot(aes(wavenumber, mean)) + 
  #geom_ribbon(aes(wavenumber, ymin = lower, ymax = upper, fill = group), alpha = 0.3) +
  geom_line(aes(col = group)) + 
  scale_x_reverse() + 
  theme_bw() +
  theme(legend.position = "top",
        axis.title = element_text(size = 20)) +
  scale_color_manual(values = col_vec,
                     name = "Taxon") +
  scale_fill_manual(values = col_vec,
                    name = "Taxon") +
  labs(x = expression(paste("Wavenumber (cm"^"-1"*")")),
       y = "Relative absorbance")

sp_tidy <- get_mean_spectra(sp, md$sample)
col_vec <- scico::scico(length(unique(md$sample)), palette = "batlow")
plot_spectra_base(sp_tidy, col_vector = col_vec, lwd = 2)
ggplot


sp_tidy <- get_mean_spectra(sp, md$family, smoothed_derivative = T, w = 9)
col_vec <- scico::scico(length(unique(md$family)), palette = "batlow")
plot_spectra_base(sp_tidy, col_vector = col_vec, lwd = 2)







md <- read.csv("MarkerPen_metadata.csv")
sp <- read.csv("MarkerPen_spectra.csv")


names_sp <- names(sp)
sp[-1] <- EMSC::EMSC(sp[-1], degree = 2)$corrected %>% data.frame()

sp[-1] <- baseline::baseline(as.matrix(sp[-1]), method = "modpolyfit", degree = 3)@corrected %>% t %>% scale %>% t %>% data.frame()
names(sp) <- names_sp

md$adj <- NA
md$base_col <- NA
adj <- c("navy", "pale", "royal", "dark", "light")
for(a in adj){
  is_a <- which(str_detect(md$sample, toupper(a)))
  md$adj[is_a] <- a
  md$base_col[is_a] <- md$sample[is_a] %>% str_replace_all(toupper(a), "") %>% tolower()
}
md$base_col <- ifelse(is.na(md$base_col), tolower(md$sample), md$base_col)
md$r_col <- ifelse(is.na(md$adj), md$base_col, paste0(md$adj, md$base_col))
md$plot_col <- sapply(md$r_col, col_str_to_hex)
unique_cols <- unique(md$r_col)
not_in_colors <- unique_cols[!unique_cols %in% colors()]

col_str_to_hex <- function(col_str){
  if(!col_str %in% colors()){
    warning("Color string not in R base colors()")
    NA
  } else {
    rgb <- col2rgb(col_str) %>% as.integer
    paste0("#", sapply(1:3, function(i){str_pad(as.hexmode(rgb[i]), pad = "0", width = 2)}) %>% paste(collapse = ""))
  }
}

# just had to quickly look up these others
md <- md %>% mutate(plot_col = case_when(r_col == "paleblue" ~ "#9bc5f8",
                                         r_col == "darkpurple" ~ "#30062e",
                                         r_col == "lightpurple" ~ "#bf17b9",
                                         r_col == "teal" ~ "#008080",
                                         .default = as.character(plot_col)))

qp_plotting_data <- md %>% select(sample, brand, brand_x_colour, plot_col) %>% unique()

sp_tidy <- get_mean_spectra(sp, grouping_vector = md$brand_x_colour)

sp_tidy$group %>% unique()

plot_spectra_base(sp_tidy, col_vector = qp_plotting_data$plot_col)

