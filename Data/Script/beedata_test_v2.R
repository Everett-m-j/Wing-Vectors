#07/08/2025
#Everett Mathiason
#working code, just need to align images directly
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
align_by_best_overlap <- function(mask1, mask2) {

  im1 <- magick2cimg(mask1)
  im2 <- magick2cimg(mask2)
  
  im1 <- resize(im1, 500, 500)
  im2 <- resize(im2, 500, 500)
  
  im1_gray <- grayscale(im1)
  im2_gray <- grayscale(im2)
  
  match <- imager::correlate(im1_gray, im2_gray)
  
  max_point <- which.max(match)
  dim_match <- dim(match)
  y_offset <- (max_point - 1) %% dim_match[1]
  x_offset <- ((max_point - 1) %/% dim_match[1]) %% dim_match[2]
  
  dx <- round(x_offset - 250)
  dy <- round(y_offset - 250)
  
  img2_aligned <- image_extent(mask2, geometry_area(500, 500, dx, dy), gravity = "NorthWest")
  
  return(list(aligned_img1 = mask1, aligned_img2 = img2_aligned))
}
print(img2_aligned)

diff_img <- image_compare(mask1, img2_aligned, metric = "AE")
print(diff_img)