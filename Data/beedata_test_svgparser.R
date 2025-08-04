# 07/31/2025
# Everett Mathiason
# EFA for wing similarity

rm(list = ls())
library(tidyverse)
library(imager)
library(magick)
library(svgparser)
library(dplyr)
library(Momocs)
library(sf)
#install.packages("devtools")
#devtools::install_github("MomX/Momocs")



svg_L <- ("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round/UCSB-IZC00065438_L_1.svg")
svg_df <- svgparser::read_svg(svg_L, obj_type = 'data.frame')

# Convert to matrix format first
coords <- as.matrix(svg_df[, c("x", "y")])  # Essential: x must be first column
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
svg_R <- "/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round/UCSB-IZC00065438_R_1.svg"
svg2_df <- svgparser::read_svg(svg_R, obj_type = 'data.frame')

coords2 <- as.matrix(svg2_df[, c("x", "y")])  # Essential: x must be first column
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
coo_draw(efi2, border='blue', col=NA)

polygon_area <- function(coords) {
  x <- coords[, 1]
  y <- coords[, 2]
  n <- length(x)
  
  # Close the polygon if not already closed
  if (!all(coords[1, ] == coords[n, ])) {
    x <- c(x, x[1])
    y <- c(y, y[1])
  }
  
  return(abs(sum(x[-1] * y[-length(y)] - x[-length(x)] * y[-1])) / 2)
}

area1 <- polygon_area(efi)  # Your first reconstruction (matrix)
area2 <- polygon_area(efi2)  # Your second reconstruction

abs_diff <- abs(area1 - area2)
percent_area <- (abs_diff / ((area1 + area2) / 2)) * 100

##curvature command
compute_curvature <- function(coords) {
  # Ensure coords is a closed curve
  if (!all(coords[1, ] == coords[nrow(coords), ])) {
    coords <- rbind(coords, coords[1, ])
  }
  
  x <- coords[, 1]
  y <- coords[, 2]
  n <- length(x)
  
  # Compute first derivatives (dx, dy)
  dx <- c(x[2:n], x[2]) - c(x[n], x[1:(n-1)])
  dy <- c(y[2:n], y[2]) - c(y[n], y[1:(n-1)])
  
  # Compute second derivatives (ddx, ddy)
  ddx <- c(dx[2:n], dx[2]) - c(dx[n], dx[1:(n-1)])
  ddy <- c(dy[2:n], dy[2]) - c(dy[n], dy[1:(n-1)])
  
  # Compute curvature: k = (dx * ddy - dy * ddx) / (dx^2 + dy^2)^(3/2)
  denom <- (dx^2 + dy^2)^(3/2)
  # Avoid division by zero
  denom[denom == 0] <- .Machine$double.eps
  curvature <- (dx * ddy - dy * ddx) / denom
  
  return(curvature)
}

curv1 <- compute_curvature(efi)
curv2 <- compute_curvature(efi2)

# Ensure same length by interpolating (if needed)
min_len <- min(length(curv1), length(curv2))
curv1 <- curv1[1:min_len]
curv2 <- curv2[1:min_len]
curv2 <- curv2*(-1)

# Compare curvature profiles
curv_diff <- abs(curv1 - curv2)
mean_curv_diff <- mean(curv_diff)
max_curv_diff <- max(curv_diff)

plot(curv1, type = "l", col = "blue", ylim = range(c(curv1, curv2)), ylab = "Curvature", main = "Curvature Profiles")
lines(curv2, col = "red")
legend("topright", legend = c("efi1", "efi2"), col = c("blue", "red"), lty = 1)

#finding polygon overlap
# Ensure polygons are closed
close_polygon <- function(mat) {
  if (!all(mat[1, ] == mat[nrow(mat), ])) {
    mat <- rbind(mat, mat[1, ])
  }
  mat
}

efi_closed <- close_polygon(efi)
efi2_closed <- close_polygon(efi2)

# Convert to sf POLYGON objects
poly1 <- st_polygon(list(efi_closed)) %>% st_sfc()
poly2 <- st_polygon(list(efi2_closed)) %>% st_sfc()

#compute area, intersection, and union of two polygons
area_1 <- st_area(poly1)
area_2 <- st_area(poly2)
intersection <- st_intersection(poly1, poly2)
union <- st_union(poly1, poly2)

area_intersection <- st_area(intersection)
area_union <- st_area(union)

# Calculate Overlap relative to union (Jaccard Index × 100)
percent_overlap <- as.numeric(area_intersection / area_union) * 100

# Or relative to shape 1 or 2
percent_overlap_wing1 <- as.numeric(area_intersection / area_1) * 100
percent_overlap_wing2 <- as.numeric(area_intersection / area_2) * 100

plot(st_geometry(poly1), border = "blue", col = adjustcolor("blue", alpha.f = 0.3))
plot(st_geometry(poly2), border = "red", col = adjustcolor("red", alpha.f = 0.3), add = TRUE)


# Function to get length and width from coordinates
get_length_width <- function(coords) {
  x_range <- range(coords[, 1], na.rm = TRUE)
  y_range <- range(coords[, 2], na.rm = TRUE)
  
  length_x <- diff(x_range)  # Horizontal extent (length)
  width_y <- diff(y_range)   # Vertical extent (width)
  
  return(list(length = length_x, width = width_y))
}

# Apply to both reconstructions
dims_efi1 <- get_length_width(efi)
dims_efi2 <- get_length_width(efi2)

# Print results
cat("Wing 1 — Length:", dims_efi1$length, " Width:", dims_efi1$width, "\n")
cat("Wing 2 — Length:", dims_efi2$length, " Width:", dims_efi2$width, "\n")

length1 <- dims_efi1$length
width1  <- dims_efi1$width
length2 <- dims_efi2$length
width2  <- dims_efi2$width

# Function to compute perimeter
compute_perimeter <- function(coords) {
  # Ensure shape is closed
  if (!all(coords[1, ] == coords[nrow(coords), ])) {
    coords <- rbind(coords, coords[1, ])
  }
  dists <- sqrt(rowSums((coords[-1, ] - coords[-nrow(coords), ])^2))
  sum(dists)
}

# Use your existing polygon_area() function for area
# Already defined in your script:
# polygon_area(coords)

# Compute for efi (wing 1)
area1 <- polygon_area(efi)
perim1 <- compute_perimeter(efi)
circ1 <- (4 * pi * area1) / (perim1^2)

# Compute for efi2 (wing 2)
area2 <- polygon_area(efi2)
perim2 <- compute_perimeter(efi2)
circ2 <- (4 * pi * area2) / (perim2^2)

# Display results
cat("Circularity of Wing 1:", round(circ1, 4), "\n")
cat("Circularity of Wing 2:", round(circ2, 4), "\n")


data_existing <- read.csv("occurrences.csv", stringsAsFactors = FALSE)
data_existing$percent_asymmetry <- NA
data_existing$area <- NA
data_existing$circularity <- NA
data_existing$length <- NA
data_existing$width <- NA
row_index <- which(data_existing$ImageName == "UCSB-IZC00055288_L_1.svg")

# Update left wing row
data_existing$percent_asymmetry[row_index] <- percent_asym
data_existing$area[row_index] <- area1
data_existing$circularity[row_index] <- circ1
data_existing$length[row_index] <- length1
data_existing$width[row_index] <- width1

# Right wing index (adjust if it has a different barcode)
row_index_R <- which(data_existing$ImageName == "UCSB-IZC00055288_R_1.svg")
data_existing$percent_asymmetry[row_index_R] <- percent_asym
data_existing$area[row_index_R] <- area2
data_existing$circularity[row_index_R] <- circ2
data_existing$length[row_index_R] <- length2
data_existing$width[row_index_R] <- width2

write.csv(data_existing, "occurrences.csv", row.names = FALSE)

