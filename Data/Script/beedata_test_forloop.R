#For Loop of Wing Imaging
rm(list = ls())
library(tidyverse)
library(imager)
library(magick)
library(svgparser)
library(dplyr)
library(Momocs)
library(sf)


# Define the folder where your SVGs are
svg_dir <- "/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/"
svg_files <- list.files(svg_dir, pattern = "\\.svg$", full.names = TRUE)

# Load your existing CSV
data_existing <- read.csv("occurrences.csv", stringsAsFactors = FALSE)

# Loop over each left wing SVG
for (svg_L in svg_files) {
  # Derive the right-wing filename
  svg_R <- sub("_L_1\\.svg$", "_R_1.svg", svg_L)
  if (!file.exists(svg_R)) next  # Skip if right-wing file is missing
 
# Load and clean coordinates
  coords_L <- svgparser::read_svg(svg_L, obj_type = 'data.frame')[, c("x", "y")] %>%
    as.matrix() %>% apply(2, as.numeric) %>% na.omit()
  coords_R <- svgparser::read_svg(svg_R, obj_type = 'data.frame')[, c("x", "y")] %>%
    as.matrix() %>% apply(2, as.numeric) %>% na.omit()
  
  clean_coords <- function(coords) {
    coords <- coords[complete.cases(coords) & is.finite(rowSums(coords)), ]
    coords[!duplicated(coords), ]
  }
  
  coords_L <- clean_coords(coords_L)
  coords_R <- clean_coords(coords_R)
  
  
  # Wrap in Momocs Out class and preprocess
  n_points <- 300
  wing_L <- Out(list(coords_L)) %>% coo_close() %>%
    coo_slide(id = which.min(coords_L[, 1])) %>% coo_smooth(3) %>%
    coo_interpolate(300)
  wing_R <- Out(list(coords_R)) %>% coo_close() %>%
    coo_slide(id = which.min(coords_R[, 1])) %>% coo_smooth(3) %>%
    coo_interpolate(300)
  
  # Procrustes alignment
  wings <- combine(wing_L, wing_R)
  aligned <- fgProcrustes(wings)
  aligned1 <- aligned$coo[[1]]
  aligned2 <- aligned$coo[[2]]
  
  # Fourier and reconstruction
  ef1 <- efourier(aligned1, nb.h = 12)
  ef2 <- efourier(aligned2, nb.h = 12)
  efi1 <- efourier_i(ef1)
  efi2 <- efourier_i(ef2)
  
  # Area, Perimeter, Circularity
  polygon_area <- function(xy) {
    if (!all(xy[1, ] == xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
    x <- xy[, 1]; y <- xy[, 2]
    abs(sum(x[-1] * y[-length(y)] - x[-length(x)] * y[-1])) / 2
  }
  
  compute_perimeter <- function(xy) {
    if (!all(xy[1, ] == xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
    sum(sqrt(rowSums((xy[-1, ] - xy[-nrow(xy), ])^2)))
  }
  
  area1 <- polygon_area(efi1)
  area2 <- polygon_area(efi2)
  perim1 <- compute_perimeter(efi1)
  perim2 <- compute_perimeter(efi2)
  circ1 <- (4 * pi * area1) / (perim1^2)
  circ2 <- (4 * pi * area2) / (perim2^2)
  percent_asym <- abs(area1 - area2) / ((area1 + area2)/2) * 100
  
  # Length and width
  get_length_width <- function(coords) {
    list(length = diff(range(coords[, 1])), width = diff(range(coords[, 2])))
  }
  dims1 <- get_length_width(efi1)
  dims2 <- get_length_width(efi2)
  
  # Overlap via sf
  close_poly <- function(mat) if (!all(mat[1, ] == mat[nrow(mat), ])) rbind(mat, mat[1, ]) else mat
  poly1 <- st_polygon(list(close_poly(efi1))) %>% st_sfc()
  poly2 <- st_polygon(list(close_poly(efi2))) %>% st_sfc()
  inter <- st_intersection(poly1, poly2)
  union <- st_union(poly1, poly2)
  percent_overlap <- as.numeric(st_area(inter) / st_area(union)) * 100
  
  # Update CSV
  left_name <- basename(svg_L)
  right_name <- basename(svg_R)
  idx_L <- which(data_existing$ImageName == left_name)
  idx_R <- which(data_existing$ImageName == right_name)
  
  data_existing[idx_L, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <- 
    list(percent_asym, area1, circ1, dims1$length, dims1$width, percent_overlap)
  data_existing[idx_R, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <- 
    list(percent_asym, area2, circ2, dims2$length, dims2$width, percent_overlap)
}


# Example filenames
left_name <- "UCSB-IZC00058724_L_1.svg"
right_name <- "UCSB-IZC00058724_R_1.svg"

# Find the rows
row_L <- which(data_existing$ImageName == left_name)
row_R <- which(data_existing$ImageName == right_name)

plot(efi1, type = "l", col = "red", asp = 1)
lines(efi2, col = "blue")

# Only process "_L_4" files (each with a matching "_R_4")
svg_L4_files <- svg_files[grepl("_L_4\\.svg$", svg_files)]
for (svg_L in svg_L4_files) {
  # Find the matching _R_4 file
  svg_R <- sub("_L_4\\.svg$", "_R_4.svg", svg_L)
  
  # Skip if R wing doesn’t exist
  if (!file.exists(svg_R)) {
    warning(paste("Right wing not found for:", svg_L))
    next
  }
  
  # Extract just the filenames for row-matching later
  svg_L_name <- basename(svg_L)
  svg_R_name <- basename(svg_R)
  
  if (grepl("_R_4", svg_L_name)) {
    coords_R[, 2] <- -coords_R[, 2]
  }
  
  # Then run your shape processing and feature extraction here...
  # Load and clean coordinates
  coords_L <- svgparser::read_svg(svg_L, obj_type = 'data.frame')[, c("x", "y")] %>%
    as.matrix() %>% apply(2, as.numeric) %>% na.omit()
  coords_R <- svgparser::read_svg(svg_R, obj_type = 'data.frame')[, c("x", "y")] %>%
    as.matrix() %>% apply(2, as.numeric) %>% na.omit()
  
  clean_coords <- function(coords) {
    coords <- coords[complete.cases(coords) & is.finite(rowSums(coords)), ]
    coords[!duplicated(coords), ]
  }
  
  coords_L <- clean_coords(coords_L)
  coords_R <- clean_coords(coords_R)
  coords_L[, 2] <- -coords_L[, 2]  # Flip Y-axis for alignment
  
  
  # Wrap in Momocs Out class and preprocess
  n_points <- 300
  wing_L <- Out(list(coords_L)) %>% coo_close() %>%
    coo_slide(id = which.min(coords_L[, 1])) %>% coo_smooth(3) %>%
    coo_interpolate(300)
  wing_R <- Out(list(coords_R)) %>% coo_close() %>%
    coo_slide(id = which.min(coords_R[, 1])) %>% coo_smooth(3) %>%
    coo_interpolate(300)
  
  # Procrustes alignment
  wings <- combine(wing_L, wing_R)
  aligned <- fgProcrustes(wings)
  aligned1 <- aligned$coo[[1]]
  aligned2 <- aligned$coo[[2]]
  
  # Fourier and reconstruction
  ef1 <- efourier(aligned1, nb.h = 12)
  ef2 <- efourier(aligned2, nb.h = 12)
  efi1 <- efourier_i(ef1)
  efi2 <- efourier_i(ef2)
  
  # Area, Perimeter, Circularity
  polygon_area <- function(xy) {
    if (!all(xy[1, ] == xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
    x <- xy[, 1]; y <- xy[, 2]
    abs(sum(x[-1] * y[-length(y)] - x[-length(x)] * y[-1])) / 2
  }
  
  compute_perimeter <- function(xy) {
    if (!all(xy[1, ] == xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
    sum(sqrt(rowSums((xy[-1, ] - xy[-nrow(xy), ])^2)))
  }
  
  area1 <- polygon_area(efi1)
  area2 <- polygon_area(efi2)
  perim1 <- compute_perimeter(efi1)
  perim2 <- compute_perimeter(efi2)
  circ1 <- (4 * pi * area1) / (perim1^2)
  circ2 <- (4 * pi * area2) / (perim2^2)
  percent_asym <- abs(area1 - area2) / ((area1 + area2)/2) * 100
  
  # Length and width
  get_length_width <- function(coords) {
    list(length = diff(range(coords[, 1])), width = diff(range(coords[, 2])))
  }
  dims1 <- get_length_width(efi1)
  dims2 <- get_length_width(efi2)
  
  # Overlap via sf
  close_poly <- function(mat) if (!all(mat[1, ] == mat[nrow(mat), ])) rbind(mat, mat[1, ]) else mat
  poly1 <- st_polygon(list(close_poly(efi1))) %>% st_sfc()
  poly2 <- st_polygon(list(close_poly(efi2))) %>% st_sfc()
  inter <- st_intersection(poly1, poly2)
  union <- st_union(poly1, poly2)
  percent_overlap <- as.numeric(st_area(inter) / st_area(union)) * 100
  
  # Update CSV
  left_name <- basename(svg_L)
  right_name <- basename(svg_R)
  idx_L <- which(data_existing$ImageName == left_name)
  idx_R <- which(data_existing$ImageName == right_name)
  
  data_existing[idx_L, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <- 
    list(percent_asym, area1, circ1, dims1$length, dims1$width, percent_overlap)
  data_existing[idx_R, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <- 
    list(percent_asym, area2, circ2, dims2$length, dims2$width, percent_overlap)
  
  # Find rows in dataset
  row_L <- which(data_existing$ImageName == svg_L_name)
  row_R <- which(data_existing$ImageName == svg_R_name)
  
  # Insert your extracted values:
  data_existing$area[row_L] <- area1
  data_existing$length[row_L] <- length1
  data_existing$width[row_L] <- width1
  data_existing$circularity[row_L] <- circ1
  data_existing$percent_asymmetry[row_L] <- percent_asym
  
  data_existing$area[row_R] <- area2
  data_existing$length[row_R] <- length2
  data_existing$width[row_R] <- width2
  data_existing$circularity[row_R] <- circ2
  data_existing$percent_asymmetry[row_R] <- percent_asym
}
# Save updated file
write.csv(data_existing, "occurrences.csv", row.names = FALSE)

# Example filenames
left_name <- "UCSB-IZC00058724_L_4.svg"
right_name <- "UCSB-IZC00058724_R_4.svg"

# Find the rows
row_L <- which(data_existing$ImageName == left_name)
row_R <- which(data_existing$ImageName == right_name)

plot(efi1, type = "l", col = "red", asp = 1)
lines(efi2, col = "blue")
