# 07/07/2025
# Everett Mathiason

library(tidyverse)
library(rsvg)
library(grImport2)
library(rsvg)
library(xml2)
library(magick)
library(dplyr)
library(tibble)

svg_folder <- "/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus"
svg_files <- list.files(path = svg_folder, pattern = "\\.svg$", full.names = TRUE)
svg_df <- tibble(
  file_name = basename(svg_files),
  file_path = svg_files,
  xml_content = lapply(svg_files, read_xml),
  image_obj = lapply(svg_files, magick::image_read_svg)
)


svg_55288_L_1 <- read_xml("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_L_1.svg"
)
svg_55288_R_1 <- read_xml("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_R_1.svg"
)
xml_ns(svg_55288_L_1)
xml_ns(svg_55288_R_1)

poly_L_1 <- xml_find_all(svg_55288_L_1, ".//d1:polygon", ns = xml_ns(svg_55288_L_1))
poly_R_1 <- xml_find_all(svg_55288_R_1, ".//d1:polygon", ns = xml_ns(svg_55288_R_1))

get_coords <- function(xml_poly) {
  points_str <- xml_attr(xml_poly, "points")
  points <- strsplit(points_str, " ")[[1]]
  coords <- do.call(rbind, strsplit(points, ","))
  coords <- as.data.frame(matrix(as.numeric(unlist(coords)), ncol = 2, byrow = TRUE))
  colnames(coords) <- c("x", "y")
  return(coords)
}

coords1 <- get_coords(poly_L_1[[1]])
coords2 <- get_coords(poly_R_1[[1]])
svg_clean <- xml_ns_strip(svg_55288_L_1)
poly_L_1 <- xml_find_all(svg_clean, ".//polygon")
length(poly_L_1)

paths <- xml_find_all(svg_clean, ".//path")
length(paths)

circles <- xml_find_all(svg_clean, ".//circle")
length(circles)

svg_clean <- xml_ns_strip(svg_55288_L_1)
paths <- xml_find_all(svg_clean, ".//path")
length(paths)

# Strip namespace
svg_clean1 <- xml_ns_strip(read_xml("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_L_1.svg"))
svg_clean2 <- xml_ns_strip(read_xml("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_R_1.svg"))

# Get path strings
paths1 <- xml_find_all(svg_clean1, ".//path")
paths2 <- xml_find_all(svg_clean2, ".//path")

d1 <- xml_attr(paths1, "d")
d2 <- xml_attr(paths2, "d")

parse_path_coords <- function(d_string) {
  # Remove command letters (M, L, Z) and split into numbers
  nums <- as.numeric(unlist(strsplit(gsub("[MLZmlz]", "", d_string), "[ ,]+")))
  coords <- matrix(nums, ncol = 2, byrow = TRUE)
  colnames(coords) <- c("x", "y")
  as.data.frame(coords)
}

# Apply to first path in each file (you could loop through all)
coords1 <- parse_path_coords(d1[1])
coords2 <- parse_path_coords(d2[1])

compare_svg_shapes <- function(svg1_path, svg2_path, width = 500, height = 500, show_diff = TRUE) {
  # Step 1: Read and render SVGs
  img1 <- image_read_svg(svg1_path, width = width, height = height)
  img2 <- image_read_svg(svg2_path, width = width, height = height)
  
  mask1 <- image_resize(mask1, "500x500!")
  mask2 <- image_resize(mask2, "500x500!")
  
  # Step 2: Convert to grayscale and then threshold to binary (black & white)
  mask1 <- image_threshold(image_convert(img1, colorspace = "gray"), "white")
  mask2 <- image_threshold(image_convert(img2, colorspace = "gray"), "white")
  
  # Step 3: Extract raw pixel data and count "black" pixels (shape area)
  data1 <- as.integer(image_data(mask1)[1,,])
  data2 <- as.integer(image_data(mask2)[1,,])
  
  area1 <- sum(data1 < 255)
  area2 <- sum(data2 < 255)
  
  # Step 4: Calculate difference
  difference <- abs(area1 - area2)
  ratio <- area2 / area1
  
  cat("Shape 1 Area (pixels):", area1, "\n")
  cat("Shape 2 Area (pixels):", area2, "\n")
  cat("Absolute Area Difference:", difference, "pixels\n")
  cat("Ratio (SVG2 / SVG1):", round(ratio, 3), "\n")
  
  # Step 5: Show visual difference image (optional)
  if (show_diff) {
    diff_img <- image_compare(mask1, mask2, metric = "AE")
    print(diff_img)
  }
}

library(magick)

compare_svg_shapes <- function(svg1_path, svg2_path, width = 500, height = 500, show_diff = TRUE) {
  # Step 1: Read SVGs
  cat("🔍 Reading SVGs...\n")
  img1 <- tryCatch(image_read_svg(svg1_path, width = width, height = height),
                   error = function(e) stop("Failed to read first SVG: ", e$message))
  img2 <- tryCatch(image_read_svg(svg2_path, width = width, height = height),
                   error = function(e) stop("Failed to read second SVG: ", e$message))
  
  # Step 2: Convert to grayscale
  cat("🎨 Converting to grayscale...\n")
  gray1 <- image_convert(img1, colorspace = "gray")
  gray2 <- image_convert(img2, colorspace = "gray")
  
  # Step 3: Threshold to binary (black-and-white mask)
  cat("🖼️ Thresholding images...\n")
  mask1 <- tryCatch(image_threshold(gray1, "white"),
                    error = function(e) stop("Failed to threshold first image: ", e$message))
  mask2 <- tryCatch(image_threshold(gray2, "white"),
                    error = function(e) stop("Failed to threshold second image: ", e$message))
  
  # Step 4: Force same size
  cat("📐 Resizing for comparison...\n")
  mask1 <- image_resize(mask1, paste0(width, "x", height, "!"))
  mask2 <- image_resize(mask2, paste0(width, "x", height, "!"))
  
  # Step 5: Count non-white pixels (black = shape area)
  cat("🔢 Counting shape pixels...\n")
  data1 <- as.integer(image_data(mask1)[1,,])
  data2 <- as.integer(image_data(mask2)[1,,])
  area1 <- sum(data1 < 255)
  area2 <- sum(data2 < 255)
  
  # Step 6: Output results
  cat("🟦 Shape 1 Area (pixels):", area1, "\n")
  cat("🟨 Shape 2 Area (pixels):", area2, "\n")
  cat("📏 Absolute Area Difference:", abs(area1 - area2), "\n")
  cat("📊 Ratio (SVG2 / SVG1):", round(area2 / area1, 3), "\n")
  
  # Step 7: Show difference image
  if (show_diff) {
    cat("🔬 Comparing images visually...\n")
    diff_img <- image_compare(mask1, mask2, metric = "AE")
    print(diff_img)
  }
}


file.exists("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_L_1.svg")
img1 <- image_read_svg("/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_L_1.svg", width = 500, height = 500)
image_info(img1)

gray1 <- image_convert(img1, colorspace = "gray")         # Convert to grayscale
mask1 <- image_threshold(gray1, "white")                  # Threshold to binary mask
image_info(mask1)
image_info(mask2)
mask1 <- image_resize(mask1, "500x500!")
mask2 <- image_resize(mask2, "500x500!")
diff_img <- image_compare(mask1, mask2, metric = "AE")
print(diff_img)
exists("diff_img")
------------------------------------------------------------------
  library(magick)

compare_svg_shapes <- function(svg1_path, svg2_path, width = 500, height = 500, show_diff = TRUE) {
  img1 <- image_read_svg(svg1_path, width = width, height = height)
  img2 <- image_read_svg(svg2_path, width = width, height = height)
  
  gray1 <- image_convert(img1, colorspace = "gray")
  gray2 <- image_convert(img2, colorspace = "gray")
  
  # Flatten and threshold to avoid transparency issues
  mask1 <- image_flatten(gray1) %>% image_threshold("white")
  mask2 <- image_flatten(gray2) %>% image_threshold("white")
  
  # Force exact same size
  size_str <- paste0(width, "x", height, "!")
  mask1 <- image_resize(mask1, size_str)
  mask2 <- image_resize(mask2, size_str)
  
  # Compare areas via pixel count
  area1 <- sum(as.integer(image_data(mask1)[1,,]) < 255)
  area2 <- sum(as.integer(image_data(mask2)[1,,]) < 255)
  
  cat("🟦 Shape 1 Area (pixels):", area1, "\n")
  cat("🟨 Shape 2 Area (pixels):", area2, "\n")
  cat("📏 Absolute Area Difference:", abs(area1 - area2), "\n")
  cat("📊 Ratio (SVG2 / SVG1):", round(area2 / area1, 3), "\n")
  
  # Now create diff image
  if (show_diff) {
    cat("🔬 Creating visual diff image...\n")
    diff_img <- image_compare(mask1, mask2, metric = "AE")
    print(diff_img)  # This should now work
  }
}
--------------------------------------
  #working code past this point: got an odd image in viewer tab
rm(list = ls())
library(magick)
library(imager)
svg_path <- "/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_L_1.svg"

file.exists(svg_path)
img1 <- image_read_svg(svg_path, width = 500, height = 500)
print(img1)
gray1 <- image_convert(img1, colorspace = "gray")
mask1 <- image_flatten(gray1) %>% image_threshold("black")
print(mask1)

svg2_path <- "/Users/everettmathiason/Documents/CCBER_Bee_Data/Data/Panurginus_polytrichus/UCSB-IZC00055288_R_1.svg"
img2 <- image_read_svg(svg2_path, width = 500, height = 500)
gray2 <- image_convert(img2, colorspace = "gray")
gray2rotate <- image_rotate(gray2, 180)
mirrored_gray2 <- image_flop(gray2rotate)
mask2 <- image_flatten(mirrored_gray2) %>% image_threshold("black")
print(mask2)

# Resize to exact match
mask1 <- image_resize(mask1, "500x500!")
mask2 <- image_resize(mask2, "500x500!")
print(mask1)

#align images
align_by_centroid <- function(mask1, mask2) {
  get_centroid <- function(mask_img) {
    mat <- as.integer(image_data(mask_img)[1,,])
    idx <- which(mat > 0, arr.ind = TRUE)  # white pixels = shape pixels
    
    # Check if idx is a matrix with rows
    if (is.null(idx) || (is.matrix(idx) && nrow(idx) == 0) || length(idx) == 0) {
      stop("⚠️ No shape pixels found — check your threshold or image content.")
    }
    
    # If idx is a vector (length 1), convert to matrix for colMeans
    if (is.vector(idx)) {
      idx <- matrix(idx, nrow = 1)
    }
    
    colMeans(idx)
  }
  
  # Get centroids
  c1 <- get_centroid(mask1)
  c2 <- get_centroid(mask2)
  
  # Calculate shift needed
  dx <- round(c1[2] - c2[2])  # x shift
  dy <- round(c1[1] - c2[1])  # y shift
  
  # Shift img2 by difference
  mask2_aligned <- image_extent(mask2, geometry_area(500, 500, dx, dy), gravity = "NorthWest")
  
  list(mask1 = mask1, mask2_aligned = mask2_aligned)
}

aligned <- align_by_centroid(mask1, mask2)
image_append(c(aligned$mask1, aligned$mask2_aligned)) %>% 
  print()

# Compare areas
area1 <- sum(as.integer(image_data(mask1)[1,,]) < 255)
area2 <- sum(as.integer(image_data(mask2)[1,,]) < 255)
cat("Area 1:", area1, "\nArea 2:", area2, "\n")

# Visual diff
diff_img <- image_compare(mask1, mask2, metric = "AE")
print(diff_img)