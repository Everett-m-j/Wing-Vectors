rm(list = ls())
library(tidyverse)
library(imager)
library(magick)
library(svgparser)
library(dplyr)
library(Momocs)
library(sf)


# Define the folder where your SVGs are
svg_dir <- "/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round"
svg_files <- list.files(svg_dir, pattern = "\\.svg$", full.names = TRUE)

# Load your existing CSV
setwd("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/")
data_existing <- read.csv("occurrences_complete.csv", stringsAsFactors = FALSE)
svg_L_files <- svg_files[grepl("_L_1\\.svg$", svg_files)]

###only use for testing for loop
svg_L <- ("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round/UCSB-IZC00059182_L_1.svg")

# Loop over each left wing SVG
for (svg_L in svg_L_files) {
  # Derive the right-wing filename
  svg_R <- sub("_L_1\\.svg$", "_R_1.svg", svg_L)
  if (!file.exists(svg_R)) next  # Skip if right-wing file is missing
  cat("Left file:", svg_L, "\n")
  cat("Right file:", svg_R, "\n")
  
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
  # Make flipped and unflipped versions of the left wing
  coords_L_unflipped <- coords_L
  coords_L_flipped <- coords_L
  coords_L_flipped[, 2] <- -coords_L_flipped[, 2]  # flip vertically
  
  # Function to create Momocs Out, align, and compute overlap/asymmetry
  analyze_pair <- function(coords_L, coords_R) {
    wing_L <- Out(list(coords_L)) %>% coo_close() %>%
      coo_slide(id = which.min(coords_L[, 1])) %>% coo_smooth(3) %>%
      coo_interpolate(300)
    wing_R <- Out(list(coords_R)) %>% coo_close() %>%
      coo_slide(id = which.min(coords_R[, 1])) %>% coo_smooth(3) %>%
      coo_interpolate(300)

    wings <- combine(wing_L, wing_R)
    aligned <- fgProcrustes(wings)
    aligned1 <- aligned$coo[[1]]
    aligned2 <- aligned$coo[[2]]
    
    ef1 <- efourier(aligned1, nb.h = 12)
    ef2 <- efourier(aligned2, nb.h = 12)
    efi1 <- efourier_i(ef1)
    efi2 <- efourier_i(ef2)
    
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
    
    
    # Overlap via sf
    close_poly <- function(mat) if (!all(mat[1, ] == mat[nrow(mat), ])) rbind(mat, mat[1, ]) else mat
    poly1 <- st_polygon(list(close_poly(efi1))) %>% st_sfc()
    poly2 <- st_polygon(list(close_poly(efi2))) %>% st_sfc()
    inter <- st_intersection(poly1, poly2)
    union <- st_union(poly1, poly2)
    percent_overlap <- as.numeric(st_area(inter) / st_area(union)) * 100
    
    # Length and width
    get_length_width <- function(coords) {
      list(length = diff(range(coords[, 1])), width = diff(range(coords[, 2])))
    }
    dims1 <- get_length_width(efi1)
    dims2 <- get_length_width(efi2)
    
    return(list(
      overlap = percent_overlap,
      asym = percent_asym,
      area1 = area1,
      area2 = area2,
      circ1 = circ1,
      circ2 = circ2,
      dims1 = dims1,
      dims2 = dims2
    ))
  }
  
  # Run analysis for both flipped and unflipped
  result_unflipped <- analyze_pair(coords_L_unflipped, coords_R)
  result_flipped <- analyze_pair(coords_L_flipped, coords_R)
  
  # Choose the version with higher overlap
  if (result_flipped$overlap > result_unflipped$overlap) {
    final_result <- result_flipped
  } else {
    final_result <- result_unflipped
  }
  # Update CSV
  left_name <- basename(svg_L)
  right_name <- basename(svg_R)
  idx_L <- which(data_existing$ImageName == left_name)
  idx_R <- which(data_existing$ImageName == right_name)
  
  # Now use `final_result` to write to your CSV:
  data_existing[idx_L, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <-
    list(final_result$asym, final_result$area1, final_result$circ1, final_result$dims1$length, final_result$dims1$width, final_result$overlap)
  
  data_existing[idx_R, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <-
    list(final_result$asym, final_result$area2, final_result$circ2, final_result$dims2$length, final_result$dims2$width, final_result$overlap)
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
svg_L <- ("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round/UCSB-IZC00058724_L_4.svg")

svg_L4_files <- svg_files[grepl("_L_4\\.svg$", svg_files)]
for (svg_L in svg_L4_files) {
  # Derive the right-wing filename
  svg_R <- sub("_L_4\\.svg$", "_R_4.svg", svg_L)
  if (!file.exists(svg_R)) next  # Skip if right-wing file is missing
  cat("Left file:", svg_L, "\n")
  cat("Right file:", svg_R, "\n")
  
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
  # Make flipped and unflipped versions of the left wing
  coords_L_unflipped <- coords_L
  coords_L_vflip <- coords_L; coords_L_vflip[, 2] <- -coords_L_vflip[, 2]
  coords_L_hflip <- coords_L; coords_L_hflip[, 1] <- -coords_L_hflip[, 1]
  coords_L_bothflip <- coords_L; coords_L_bothflip[, 1:2] <- -coords_L_bothflip[, 1:2]
  
  # Function to create Momocs Out, align, and compute overlap/asymmetry
  analyze_pair <- function(coords_L, coords_R) {
    wing_L <- Out(list(coords_L)) %>% coo_close() %>%
      coo_slide(id = which.min(coords_L[, 1])) %>% coo_smooth(3) %>%
      coo_interpolate(300)
    wing_R <- Out(list(coords_R)) %>% coo_close() %>%
      coo_slide(id = which.min(coords_R[, 1])) %>% coo_smooth(3) %>%
      coo_interpolate(300)
    print(wing_L)
    wing_L$ldk <- list()
    wing_R$ldk <- list()
    wings <- combine(wing_L, wing_R)
    aligned <- fgProcrustes(wings)
    aligned1 <- aligned$coo[[1]]
    aligned2 <- aligned$coo[[2]]
    
    ef1 <- efourier(aligned1, nb.h = 12)
    ef2 <- efourier(aligned2, nb.h = 12)
    efi1 <- efourier_i(ef1)
    efi2 <- efourier_i(ef2)
    best_overlap <- -Inf
    best_angle <- 0
    best_efi1 <- efi1
    
    for (angle in seq(0, 340, by = 20)) {
      theta <- angle * pi / 180
      rot_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
      rotated_efi1 <- efi1 %*% rot_matrix
      
      # Close polygon
      close_poly <- function(mat) if (!all(mat[1, ] == mat[nrow(mat), ])) rbind(mat, mat[1, ]) else mat
      poly1 <- st_polygon(list(close_poly(rotated_efi1))) %>% st_sfc()
      poly2 <- st_polygon(list(close_poly(efi2))) %>% st_sfc()
      
      inter <- st_intersection(poly1, poly2)
      union <- st_union(poly1, poly2)
      if (length(inter) > 0 && length(union) > 0) {
        percent_overlap <- as.numeric(st_area(inter) / st_area(union)) * 100
        if (!is.na(percent_overlap) && percent_overlap > best_overlap) {
          best_overlap <- percent_overlap
          best_angle <- angle
          best_efi1 <- rotated_efi1
        }
      }
    }
    
    # Area & perimeter
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
    
    # Overlap via sf
    close_poly <- function(mat) if (!all(mat[1, ] == mat[nrow(mat), ])) rbind(mat, mat[1, ]) else mat
    poly1 <- st_polygon(list(close_poly(efi1))) %>% st_sfc()
    poly2 <- st_polygon(list(close_poly(efi2))) %>% st_sfc()
    inter <- st_intersection(poly1, poly2)
    union <- st_union(poly1, poly2)
    percent_overlap <- as.numeric(st_area(inter) / st_area(union)) * 100
    
    # Length & width
    get_length_width <- function(coords) {
      list(length = diff(range(coords[, 1])), width = diff(range(coords[, 2])))
    }
    dims1 <- get_length_width(efi1)
    dims2 <- get_length_width(efi2)
    
    return(list(
      overlap = percent_overlap,
      asym = percent_asym,
      area1 = area1,
      area2 = area2,
      circ1 = circ1,
      circ2 = circ2,
      dims1 = dims1,
      dims2 = dims2,
      efi1 = efi1,
      efi2 = efi2  # <–– RETURN THESE
    ))
  }
  
  # Run analysis for both flipped and unflipped
  result_unflipped <- analyze_pair(coords_L_unflipped, coords_R)
  result_vflip     <- analyze_pair(coords_L_vflip, coords_R)
  result_hflip     <- analyze_pair(coords_L_hflip, coords_R)
  result_bothflip  <- analyze_pair(coords_L_bothflip, coords_R)
  
  all_results <- list(
    unflipped = result_unflipped,
    vflip = result_vflip,
    hflip = result_hflip,
    bothflip = result_bothflip
  )
  
  best_version <- names(which.max(sapply(all_results, function(res) res$overlap)))
  final_result <- all_results[[best_version]]
  # Update CSV
  left_name <- basename(svg_L)
  right_name <- basename(svg_R)
  idx_L <- which(data_existing$ImageName == left_name)
  idx_R <- which(data_existing$ImageName == right_name)
  
  # Now use `final_result` to write to your CSV:
  data_existing[idx_L, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <-
    list(final_result$asym, final_result$area1, final_result$circ1, final_result$dims1$length, final_result$dims1$width, final_result$overlap)
  
  data_existing[idx_R, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <-
    list(final_result$asym, final_result$area2, final_result$circ2, final_result$dims2$length, final_result$dims2$width, final_result$overlap)
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
lines(efi2, col = "black")

plot(best_efi1, type = "l", col = "black", asp = 1, main = "Unflipped")
lines(efi2, col = "red")

plot(result_vflip$efi1, type = "l", col = "black", asp = 1, main = "Vertically flipped")
lines(result_vflip$efi2, col = "red")

plot(result_hflip$efi1, type = "l", col = "black", asp = 1, main = "Horizontally flipped")
lines(result_hflip$efi2, col = "red")

plot(result_bothflip$efi1, type = "l", col = "black", asp = 1, main = "Both flipped")
lines(result_bothflip$efi2, col = "red")

view(data_existing)

##############################
#non flipping code
rm(list = ls())
library(tidyverse)
library(imager)
library(magick)
library(svgparser)
library(dplyr)
library(Momocs)
library(sf)

# testing orientation of wings
svg_L <- ("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round/UCSB-IZC00058724_L_4.svg")
svg_R <- ("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round/UCSB-IZC00058724_R_4.svg")

coords_L <- svgparser::read_svg(svg_L, obj_type = 'data.frame')[, c("x", "y")] %>%
  as.matrix() %>% apply(2, as.numeric) %>% na.omit()
coords_R <- svgparser::read_svg(svg_R, obj_type = 'data.frame')[, c("x", "y")] %>%
  as.matrix() %>% apply(2, as.numeric) %>% na.omit()
wing_L <- Out(list(coords_L)) %>% coo_close() %>%
  coo_slide(id = which.min(coords_L[, 1])) %>% coo_smooth(3) %>%
  coo_interpolate(300)
wing_R <- Out(list(coords_R)) %>% coo_close() %>%
  coo_slide(id = which.min(coords_R[, 1])) %>% coo_smooth(3) %>%
  coo_interpolate(300)
wings <- combine(wing_L, wing_R)
aligned <- fgProcrustes(wings)
aligned1 <- aligned$coo[[1]]
aligned2 <- aligned$coo[[2]]

plot(aligned1, type = "l", col = "red", asp = 1)
lines(aligned2, col = "blue")

# Clockwise rotation by 90 degrees
rotate_shape <- function(coords, angle_degrees) {
  theta <- angle_degrees * pi / 180  # Convert to radians
  rotation_matrix <- matrix(
    c(cos(theta), sin(theta), -sin(theta), cos(theta)),
    ncol = 2, byrow = TRUE
  )
  rotated_coords <- coords %*% rotation_matrix
  return(rotated_coords)
}

aligned1 <- rotate_shape(aligned1, angle_degrees = 180)

cat("Best overlap at angle:", aligned1, "degrees\n")

plot(aligned1, type = "l", col = "red", asp = 1)
lines(aligned2, col = "blue")
###fixed vrsion


# Define the folder where your SVGs are
svg_dir <- "/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round"
svg_files <- list.files(svg_dir, pattern = "\\.svg$", full.names = TRUE)

# Load your existing CSV
setwd("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data")
data_existing <- read.csv("occurrences_complete.csv", stringsAsFactors = FALSE)
data_existing <- data_existing %>% select(-Total.Complete.Agr, -Total.Complete.Dev, -Total.Complete.For, -Notes., -Asymmetry.,)
svg_L_files <- svg_files[grepl("_L_1\\.svg$", svg_files)]

###only use for testing for loop
svg_L <- ("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round/UCSB-IZC00065438_L_1.svg")

# Loop over each left wing SVG
for (svg_L in svg_L_files) {
  # Derive the right-wing filename
  svg_R <- sub("_L_1\\.svg$", "_R_1.svg", svg_L)
  if (!file.exists(svg_R)) next  # Skip if right-wing file is missing
  cat("Left file:", svg_L, "\n")
  cat("Right file:", svg_R, "\n")
  
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
  
  # Function to create Momocs Out, align, and compute overlap/asymmetry
    wing_L <- Out(list(coords_L)) %>% coo_close() %>%
      coo_slide(id = which.min(coords_L[, 1])) %>% coo_smooth(3) %>%
      coo_interpolate(300)
    wing_R <- Out(list(coords_R)) %>% coo_close() %>%
      coo_slide(id = which.min(coords_R[, 1])) %>% coo_smooth(3) %>%
      coo_interpolate(300)
    
    wings <- combine(wing_L, wing_R)
    aligned <- fgProcrustes(wings)
    aligned1 <- aligned$coo[[1]]
    aligned2 <- aligned$coo[[2]]
    
    ef1 <- efourier(aligned1, nb.h = 12)
    ef2 <- efourier(aligned2, nb.h = 12)
    efi1 <- efourier_i(ef1)
    efi2 <- efourier_i(ef2)
    #best_overlap <- -Inf
    #best_angle <- 0
   # best_efi1 <- NULL
    
   # angles <- seq(-180, 180, by = 10)  # Or Try every 5 degrees
    
 #   for (angle in angles) {
      #theta <- angle * pi / 180
      #rotation_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol = 2)
      #rotated <- efi1 %*% rotation_matrix
      
     # poly1 <- st_polygon(list(rbind(rotated, rotated[1, ]))) %>% st_sfc()
     # poly2 <- st_polygon(list(rbind(efi2, efi2[1, ]))) %>% st_sfc()
      
     # inter <- st_intersection(poly1, poly2)
    #  union <- st_union(poly1, poly2)
      
    #  if (length(inter) == 0 || st_area(union) == 0) next
      
    #  overlap <- as.numeric(st_area(inter) / st_area(union)) * 100
      
    #  if (overlap > best_overlap) {
      #  best_overlap <- overlap
      #  best_angle <- angle
      #  best_efi1 <- rotated
      #}
    #}
    
    #cat("Best overlap at angle:", best_angle, "degrees\n")
    
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
    
    
    # Overlap via sf
    close_poly <- function(mat) if (!all(mat[1, ] == mat[nrow(mat), ])) rbind(mat, mat[1, ]) else mat
    poly1 <- st_polygon(list(close_poly(efi1))) %>% st_sfc()
    poly2 <- st_polygon(list(close_poly(efi2))) %>% st_sfc()
    inter <- st_intersection(poly1, poly2)
    union <- st_union(poly1, poly2)
    percent_overlap <- as.numeric(st_area(inter) / st_area(union)) * 100
    
    # Length and width
    get_length_width <- function(coords) {
      list(length = diff(range(coords[, 1])), width = diff(range(coords[, 2])))
    }
    dims1 <- get_length_width(efi1)
    dims2 <- get_length_width(efi2)
  
  # Update CSV
  left_name <- basename(svg_L)
  right_name <- basename(svg_R)
  idx_L <- which(data_existing$ImageName == left_name)
  idx_R <- which(data_existing$ImageName == right_name)
  
  # Now use `final_result` to write to your CSV:
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
svg_L <- ("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/Tracings_1st_Round/UCSB-IZC00058724_L_4.svg")

svg_L4_files <- svg_files[grepl("_L_4\\.svg$", svg_files)]
for (svg_L in svg_L4_files) {
  # Derive the right-wing filename
  svg_R <- sub("_L_4\\.svg$", "_R_4.svg", svg_L)
  if (!file.exists(svg_R)) next  # Skip if right-wing file is missing
  cat("Left file:", svg_L, "\n")
  cat("Right file:", svg_R, "\n")
  
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
  
  # Function to create Momocs Out, align, and compute overlap/asymmetry
  wing_L <- Out(list(coords_L)) %>% coo_close() %>%
    coo_slide(id = which.min(coords_L[, 1])) %>% coo_smooth(3) %>%
    coo_interpolate(300)
  wing_R <- Out(list(coords_R)) %>% coo_close() %>%
    coo_slide(id = which.min(coords_R[, 1])) %>% coo_smooth(3) %>%
    coo_interpolate(300)
  
  wings <- combine(wing_L, wing_R)
  aligned <- fgProcrustes(wings)
  aligned1 <- aligned$coo[[1]]
  aligned2 <- aligned$coo[[2]]
  
  ef1 <- efourier(aligned1, nb.h = 12)
  ef2 <- efourier(aligned2, nb.h = 12)
  efi1 <- efourier_i(ef1)
  efi2 <- efourier_i(ef2)
  

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
  
  
  # Overlap via sf
  close_poly <- function(mat) if (!all(mat[1, ] == mat[nrow(mat), ])) rbind(mat, mat[1, ]) else mat
  poly1 <- st_polygon(list(close_poly(efi1))) %>% st_sfc()
  poly2 <- st_polygon(list(close_poly(efi2))) %>% st_sfc()
  inter <- st_intersection(poly1, poly2)
  union <- st_union(poly1, poly2)
  percent_overlap <- as.numeric(st_area(inter) / st_area(union)) * 100
  
  # Length and width
  get_length_width <- function(coords) {
    list(length = diff(range(coords[, 1])), width = diff(range(coords[, 2])))
  }
  dims1 <- get_length_width(efi1)
  dims2 <- get_length_width(efi2)
  
  # Update CSV
  left_name <- basename(svg_L)
  right_name <- basename(svg_R)
  idx_L <- which(data_existing$ImageName == left_name)
  idx_R <- which(data_existing$ImageName == right_name)
  
  # Now use `final_result` to write to your CSV:
  data_existing[idx_L, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <-
    list(percent_asym, area1, circ1, dims1$length, dims1$width, percent_overlap)
  
  data_existing[idx_R, c("percent_asymmetry", "area", "circularity", "length", "width", "percent_overlap")] <-
    list(percent_asym, area2, circ2, dims2$length, dims2$width, percent_overlap)
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
lines(efi2, col = "black")
