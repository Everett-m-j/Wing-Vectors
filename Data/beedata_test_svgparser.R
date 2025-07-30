# 07/08/2025
# Everett Mathiason
# EFA for wing similarity

rm(list = ls())
library(tidyverse)
library(imager)
library(magick)
library(svgparser)
library(dplyr)
library(Momocs)
#install.packages("devtools")
#devtools::install_github("MomX/Momocs")



svg_L <- ("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_L_1.svg")
svg_df <- svgparser::read_svg(svg_L, obj_type = 'data.frame')

nrow(svg_df)
ggplot(svg_df) + 
  geom_path(aes(x, y, colour = name, group = interaction(elem_idx, path_idx))) +
  scale_y_reverse()

# Threshold for black shape
binary <- img_cimg > .2
plot(binary)

# Label and isolate largest objects
lab <- label(binary)
label_vals <- table(lab)

# Get the largest object (excluding background)
main_obj <- as.integer(names(sort(label_vals[names(label_vals) != "0"], decreasing = TRUE)[1]))

# Create masks for outline
shape_mask1 <- lab == main_obj[1]

# Extract contours - using higher n for better resolution
contours_list1 <- contours(shape_mask1, n = 1000)

# Check the structure of the contours
str(contours_list1)

df <- data.frame(
  x = contours_list1[[1]]$x,
  y = contours_list1[[1]]$y,
  outline = "outline"
)
str(df)

# Convert to matrix format first
coords <- as.matrix(df[, c("x", "y")])  # Essential: x must be first column
coords <- na.omit(coords)
coords <- apply(coords, 2, as.numeric) %>% na.omit()
coords_clean <- coords[is.finite(rowSums(coords)), ]
coords_clean <- apply(coords_clean, 2, as.numeric)
coords_clean <- coords_clean[complete.cases(coords_clean), ]  # Remove NA
coords_clean <- coords_clean[is.finite(rowSums(coords_clean)), ]  # Remove Inf
coords_clean <- coords_clean[!duplicated(coords_clean), ]  # Remove duplicates
any(is.na(coords))  # Should return FALSE
str(coords_clean)
head(coords_clean)
tail(coords_clean)
nrow(coords_clean)  # Should be >= 100

# for second svg
svg_path2 <- "/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_R_1.svg"
img2 <- image_read_svg(svg_path2, width = 500, height = 500) %>%
  image_background("white") %>%  # remove transparency
  image_flatten() %>%            # blend alpha over white
  image_convert(colorspace = "gray")
img_cimg2 <- magick2cimg(img2)

# Threshold for black shape
binary2 <- img_cimg2 > .2
plot(binary2)

# Label and isolate largest objects
lab2 <- label(binary2)
label_vals2 <- table(lab2)

# Get the largest object (excluding background)
main_obj2 <- as.integer(names(sort(label_vals2[names(label_vals2) != "0"], decreasing = TRUE)[1]))

# Create masks for outline
shape_mask2 <- lab2 == main_obj2[1]

# Extract contours - using higher n for better resolution
contours_list2 <- contours(shape_mask2, n = 1000)

# Check the structure of the contours
str(contours_list2)

df2 <- data.frame(
  x = contours_list2[[1]]$x,
  y = contours_list2[[1]]$y,
  outline = "outline"
)
str(df2)

# Convert to matrix format first
coords2 <- as.matrix(df2[, c("x", "y")])  # Essential: x must be first column
any(is.na(coords2))  # Should return FALSE


# Remove NA/NaN rows if they exist
coords2 <- na.omit(coords2)
coords2 <- apply(coords2, 2, as.numeric) %>% na.omit()
clean_coords <- function(coords) {
  coords <- as.matrix(coords)
  coords <- coords[complete.cases(coords), ]                     # Remove NA
  coords <- coords[is.finite(rowSums(coords)), ]                 # Remove Inf/NaN
  coords <- coords[!duplicated(coords), ]                        # Remove duplicates
  coords
}

coords_clean <- clean_coords(coords)
coords_clean2 <- clean_coords(coords2)
coords_clean_flipped <- coords_clean
coords_clean[, 2] <- -coords_clean_flipped[, 2]  # Flip y-axis


str(coords_clean2)
head(coords_clean2)
tail(coords_clean2)
nrow(coords_clean2)  # Should be >= 100

# Check coords_L
str(coords_clean)
summary(coords_clean)
any(is.na(coords_clean))
nrow(coords_clean)

# Check coords_R
str(coords_clean2)
summary(coords_clean2)
any(is.na(coords_clean2))
nrow(coords_clean2)

plot(coords_clean, type = "l", asp = 1, main = "Left Wing")
plot(coords_clean2, type = "l", asp = 1, main = "Right Wing")

head(coords_clean)
tail(coords_clean)
coo_plot(coords_clean)
coo_plot(coords_clean2)
class(coords_clean2)


##########trying to make comparable list to ef with ef_list
# 1. Wrap coords in Out class (Momocs expects this)
n_points <- 300
wing_L <- Out(list(coords_clean)) %>%
  coo_close() %>%
  coo_slide(id = which.min(coords_clean[, 1])) %>%
  coo_smooth(3) %>%
  coo_interpolate(300)

wing_R <- Out(list(coords_clean2)) %>%
  coo_close() %>%
  coo_slide(id = which.min(coords_clean2[, 1])) %>%
  coo_smooth(3) %>%
  coo_interpolate(300)

# Then interpolate both to have same number of coordinates
wing_L <- coo_interpolate(wing_L, n_points)
wing_R <- coo_interpolate(wing_R, n_points)
str(wing_L)

#combine them for fgProcrustes
wing_cells <- combine(wing_L, wing_R)
wings_aligned <- fgProcrustes(wing_cells)
stack(wings_aligned, title = "Stacked Procrustes-Aligned Wings")
aligned1 <- wings_aligned$coo[[1]]
aligned2 <- wings_aligned$coo[[2]]
ef1 <- efourier(aligned1, nb.h = 12)
ef2 <- efourier(aligned2, nb.h = 12)

efi <- efourier_i(ef1)
coo_draw(efi, border='red', col=NA)

efi2 <- efourier_i(ef2)