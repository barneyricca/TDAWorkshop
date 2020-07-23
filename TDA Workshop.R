# ########################################### #
## Setup                                   ####
# ########################################### #
{         # Put the cursor on this line, and run it; because
  #  of the curly braces, all the code within
  #  the braces will be run.
  c("bmp",          # Read .BMP files
    "data.table",   # fread() for fast data input
    "dplyr",        # Data wrangling
    "dtplyr",       # For faster tidyverse data wrangling
    "dtw",          # Dynamic time warping
    "gridExtra",    # For better table printing
    "here",         # To assist with folder structure
    "igraph",       # Needed for some raster functions
    "raster",       # For image processing and coverings
    "R.matlab",     # To read matlab data from Zhang et al. (2020)
    "Rcpp",         # For faster processing of Betti numbers
    "TDA") ->       # Topological Data Analysis package
    package_names

  for (package_name in package_names) {
    if (!is.element(package_name, installed.packages()[, 1])) {
      install.packages(package_name, repos = "http://lib.stat.cmu.edu/R/CRAN")
    }
    library(package_name,
            character.only = TRUE,
            quietly = TRUE,
            verbose = FALSE
    )
  }
  rm(list = c("package_names", "package_name"))

  # Reproducibility:
  set_here()      # This helps with folder structure when shifting computers

  # Preferences:
  options(show.signif.stars = FALSE) # Do not show significance stars! Showing
  #  significance stars encourages the
  #  conflation of significant and effect
  #  size.
  options(digits = 5)                # Round to 5 digits
}

# ########################################### #
## Helper Functions                        ####
# ########################################### #
mat_cover <- function(bwmat, rad) {
  # This function accepts a square matrix of 0s and 1s, and covers
  #  it with circles of radius rad.
  # The units of both of these will be in pixels.
  if(ncol(bwmat) != nrow(bwmat)) {
    return(NA)
  }
  bwmat -> new_bw_mat

  # Probably it is fastest to create a mask of changes, and then
  #  apply it. Here's one quarter of it.
  # Note that the axes are at row & column, and so are not
  #  in the matrix.
  #
  # Might go faster by doing only 1/8 of the circle and
  #  then using that 8 times, plus axes and 45 degree lines.
  matrix(1,
         nrow = rad,
         ncol = rad) -> mask
  rad * rad -> rad2
  for(i in 1:rad) {
    for(j in 1:rad) {
      if(i*i+j*j <= rad2) {
        0 -> mask[i,j]
      }
    }
  }

  # The application may call for Rcpp. For now, let's
  #  try to avoid that. See how long this takes
  nrow(new_bw_mat) -> nr
  for(r in rad:(nr - rad)) {
    for(c in rad:(nr - rad)) {
      if(bwmat[r,c] == 0) {
        # 4 quadrants, then do the axes
        for(rz in 1:rad) {
          for(cz in 1:rad) {
            if(mask[rz, cz] == 0) {
              0 -> new_bw_mat[r - rz, c - cz]
              0 -> new_bw_mat[r - rz, c + cz]
              0 -> new_bw_mat[r + rz, c - cz]
              0 -> new_bw_mat[r + rz, c + cz]
            }
          }
        }
        # Now for the axes
        for(ax in 1:rad) {
          0 -> new_bw_mat[r - ax, c]
          0 -> new_bw_mat[r + ax, c]
          0 -> new_bw_mat[r, c + ax]
          0 -> new_bw_mat[r, c - ax]
        }
      }
    }
  }

  return(new_bw_mat)
}

cover_figure_out <- function(bwmat, file_name) {
  png(here(paste0("/Images/",file_name)),
      width = 5,
      height = 5,
      res = 300,
      units = "in")
  plot(c(1, nrow(bwmat)), c(1, nrow(bwmat)),
       type = 'n',
       xlab = "",
       ylab = "",
       asp = 1)
  as.raster(bwmat) -> image
  rasterImage(bwmat, 1, 1, nrow(bwmat), nrow(bwmat),
              interpolate = FALSE)
  dev.off()
  return()
}


# ########################################### #
## A-B coverings                           ####
# ########################################### #
# Here's an example of the covering stuff
# This is just a quick image that I made for this purpose
if(is.bmp(here("AB16.bmp"))) {
  read.bmp(here("AB16.bmp")) -> AB_arr
}

# Do some cleaning of the image; convert it to a black-and-white
#  image. You don't care about this except to know that the next
#  command converts an almost black-and-white image to an
#  actual black-and-white image.
floor(AB_arr[,,] / 128.0)[,,1] -> bw_mat

# Now, create a bunch of figures with coverings of different
#  radii balls.
#
#
# Nota bene:
#  In order to run the next command, you must have a folder "Images"
#  in the current working folder/directory!
#
#
cover_figure_out(bw_mat, "AB.png")

# The following loop takes about 5 minutes on my laptop:
for(index in 1:70) {
  mat_cover(bw_mat, index) -> mat1
  cover_figure_out(mat1, paste0("AB-",index,".png"))
}


# ########################################### #
## Persistence diagram of A-B covering     ####
# ########################################### #
if(is.bmp(here("AB16.bmp"))) {
  read.bmp(here("AB16.bmp")) -> AB_arr
}

# Do some cleaning of the image; convert it to a black-and-white
#  image.
floor(AB_arr[,,] / 128.0)[,,1] -> bw_mat

# Convert the image to a point cloud (i.e., a list of points
#  that are in the image)
which(as.matrix(bw_mat == 0)) -> zeroes
as.numeric(zeroes %% nrow(bw_mat)) -> row_num
as.numeric(zeroes %/% nrow(bw_mat)) -> col_num
cbind(row_num, col_num) -> X1

# Set the limits to include all the points in X1. This can be done
#  using min() and max(), but I already know the size of bw_mat.
c(1, 450) -> Xlim
c(1, 450) -> Ylim
cbind(Xlim, Ylim) -> lim

1 -> by    # Equivalent to the step size for radius change
# Using too small a "by" creates problems. For example,
#  by = 0.5 in this analysis created more than 9000 features and
#  took about 12 minutes, as opposed to creating only 133
#  features and taking about 8 seconds with b = 1. This sample should
#  only have 34 features; in the latter case, 99 of the features were
#  holes that existed only momentarily. The other 34 features were
#  the correct ones.

gridDiag(X = X1,distFct, lim = lim, by = by,
         printProgress = TRUE) -> bw_grid

plot(bw_grid[["diagram"]])

# confidence set
B <- 10       ## the number of bootstrap iterations should be higher!
## this is just an example
alpha <- 0.05
0.3 -> h

# With B = 10, started at 1:18. Finished 1:20
set.seed(42)
cc <- bootstrapDiagram(X1, kde, lim = lim, by = by, sublevel = FALSE, B = B,
                       alpha = alpha, dimension = 1, printProgress = TRUE, h = h)

# The band is too small to see on the diagram, so I increased it's size.
plot(bw_grid[["diagram"]], band = 10000 * cc)


# ########################################### #
## Time series data: Walking data example  ####
# ########################################### #
fread(file = here("Walking sample.csv"),
      header = TRUE) -> walk1_df

plot(walk1_df$Time, walk1_df$ankleAngle,
     type = 'l')

walk1_df$ankleAngle -> funval

# Now, work with this
walk.diag <- gridDiag(FUNvalues = funval, sublevel = TRUE)
walk.diag$diagram

min(funval)
# -43.4
max(funval)
# 10.3
# Limits must include the entire thing, so the next works. You need to
#  be careful checking this because of positive/negative values.
c(1.1 * min(funval), 1.1*max(funval)) -> Ylim
by <- 0.1

1 -> maxdimension

10 -> B
set.seed(42)
cc <- bootstrapDiagram(X = funval, kde,
                       lim = Ylim, # Only one dimension, not two
                       by = by,
                       sublevel = FALSE,
                       B = B,
                       alpha = 0.05,
                       dimension = 0,  # Only dimension 0 features exit
                       printProgress = TRUE,
                       h = 0.3)
# The 200 is arbitrary to make things a bit easier. This needs some
#  further investigation in the context of smoothing. With that
#  value, however, only 11 of the 37 features are significant.
plot(walk.diag[['diagram']], band = 200 * cc)
