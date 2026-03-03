# Least cost analysis - converting from raster to terra
# Jeff Oliver
# jcoliver@arizona.edu
# 2026-03-02

################################################################################
################# MODERNIZED SCRIPT FOR FIELD ET AL. 2021 ######################
################################################################################

library(terra)
# library(sf)
library(movecost)
library(leastcostpath)
library(tidyverse)
# library(spatialEco)
# library(magrittr)
# library(sp)
# library(devtools)

# Set global projection
master_crs <- "EPSG:26913" # UTM Zone 13N NAD83

###########################################################
################## STEP 1: IMPORT DATA ####################
###########################################################

# Ensure these files exist in your directory
origins <- terra::vect("raw_data/origin_locations.shp") 
sites <- terra::vect("raw_data/Sites_Test.shp")

destination1 <- sites[sites$Name == "Corridor 2 destination", ]
destination2 <- sites[sites$Name == "Corridor 3 destination", ]

# Import DEM using terra
dem <- terra::rast("raw_data/DEM_corridor1.tif")
# Convert to a RasterLayer object for use with...?
dem_raster <- raster::raster(dem)

# Reproject to master CRS if necessary
origins      <- terra::project(origins, master_crs)
sites        <- terra::project(sites,   master_crs)
destination1 <- terra::project(destination1, master_crs)
dem          <- terra::project(dem, master_crs)

###########################################################
####### STEP 2: CREATE CONDUCTANCE SURFACES  ##############
###########################################################

# 2.1 Terrain Ruggedness (using terra)
tri_dem <- terra::terrain(dem, v = "TRI")

# 2.2 Create Slope-based Cost Surfaces (leastcostpath)
# neigh <- 16
# Since default neighbor directions (16) is used, we don't need to declare it
tobler_cs   <- leastcostpath::create_slope_cs(dem, cost_function = "tobler")
ic_off_m_cs <- leastcostpath::create_slope_cs(dem, cost_function = "irmischer-clarke offpath male")
campbell_cs <- leastcostpath::create_slope_cs(dem, cost_function = "campbell 2019 50")
herzog_cs   <- leastcostpath::create_slope_cs(dem, cost_function = "herzog")
ls_cs       <- leastcostpath::create_slope_cs(dem, cost_function = "llobera-sluckin")

# 2.3 Pandolf Surface using movecost
# movecost handles the complex Pandolf/Santee slope-switching logic automatically
# Note: movecost creates its own transition layers
W <- 63.5      # body mass
L <- 20        # load mass (kg)
V <- 1.26      # velocity in km/h 
N <- 1.2       # terrain coefficient

# We use 'side_t' to explicitly define the terrain coefficient 
# and 'v' for velocity to avoid any abbreviation confusion.

# Compute slope raster (rise/run, dimensionless) with terra
# terrain() returns slope in degrees by default; convert to tangent (rise/run)
slope_deg <- terra::terrain(dem, v = "slope", unit = "degrees", neighbors = 8)
slope_tan <- tan(slope_deg * pi / 180)   # dimensionless slope (rise/run)

##### Define cost functions
# Pandolf cost function (positive slopes)
# Energy expenditure (W/kg) following Pandolf et al. 1977
# @param s numeric slope
pandolf_energy <- function(s) {
  # s = slope as decimal (rise/run), converted to % inside formula (* 100)
  1.5 * W + 2.0 * (W + L) * (L / W)^2 +
    N * (W + L) * (1.5 * V^2 + 0.35 * V * abs(s) * 100)
}

# Santee et al. (2001) cost function (negative slopes); requires previously 
# defined pandolf_energy() function
# @param s numeric slope
santee_energy <- function(s) {
  pandolf_energy(s) -
    N * (V * abs(s) * 100 * (W + L) / 3.5 -
           (W + L) * (abs(s) * 100 + 6)^2 / W +
           (25 - V^2))
}

# Apply functions to slope raster; conductance = 1 / metabolic cost
pandolf_pos_r <- terra::app(slope_tan, function(s) {
  ifelse(s >= 0, 1 / pandolf_energy(s), NA)
})
# TODO: this is returning a raster of all NaN. Maybe that's OK, because the 
# ifelse only uses the Santee cost function if the slope is negative
pandolf_neg_r <- terra::app(slope_tan, function(s) {
  ifelse(s < 0,  1 / santee_energy(s),  NA)
})

# Merge positive- and negative-slope surfaces
pandolf_cs_r <- terra::cover(pandolf_pos_r, pandolf_neg_r)

# Remove variables no longer needed
rm(pandolf_pos_r, pandolf_neg_r, slope_deg, slope_tan, W, L, N, V)

###########################################################
########### STEP #3: CALCULATE LEAST COST PATHS ###########
########### FROM SERIES OF ORIGIN POINTS TO  ##############
########### SEPARATE DESTINATIONS AND COMPUTE #############
########### DETERMINING VARIABLES (TERRAIN, DISTANCE) #####
###########################################################

num_origins <- nrow(origins)   # terra::nrow() or length() both work on SpatVector
# Build results table
result_destination1 <- as.data.frame(matrix(NA, nrow = num_origins, ncol = 9))
colnames(result_destination1) <- c(
  "time (sq.m)", "energy (sq.m)", "overlap (sq.m)",
  "mean time and energy (sq.m)", "overlap as percent of average",
  "mean TRI within 1km of time lcps",    "median TRI within 1km of time lcps",
  "mean TRI within 1km of energy lcps",  "median TRI within 1km of energy lcps"
)

# Iterate over the different origins and...do some stuff
for (i in seq_len(num_origins)) {
  # Extract the ith origin from SpatVector
  origin_i <- origins[i, ]

  # A little message to assuage our paranoia
  message("Calculating Tobler LCPs for origin ", i)
    
  # Calculate LCPs
  # leastcostpath 2.x: create_lcp() signature is largely unchanged but
  # accepts SpatVector (terra) instead of Spatial* objects.
  
  t_lcp  <- leastcostpath::create_lcp(tobler_cs,   origin_i, destination1)
  # ic_lcp <- leastcostpath::create_lcp(ic_off_m_cs, origin_i, destination1)
  # c_lcp  <- leastcostpath::create_lcp(campbell_cs, origin_i, destination1)
  # h_lcp  <- leastcostpath::create_lcp(herzog_cs,   origin_i, destination1)
  # ls_lcp <- leastcostpath::create_lcp(ls_cs,       origin_i, destination1)
  
  # Let's try the movecost::movecost approach, too
  p_lcp <- movecost::movecost(
    dtm    = dem_raster,
    origin = sf::st_as_sf(origin_i),
    destin = sf::st_as_sf(destination1),
    funct  = "p", # Was "Pa"
    W      = 63.5,
    L      = 20,
    V      = 0.35,
    N      = 1.2
  )$LCPs
  
  # The above call to movecost::movecost errors out. Most likely at line 1126 
  # of the movecost function, in a call to gdistance::accCost
  # accum_final <- gdistance::accCost(Conductance, sp::coordinates(origin))
  #
  # If funct = "Pa" it gives this error:
  # error in evaluating the argument 'x' in selecting a method for function 'accCost': 
  #    object 'Conductance' not found
  # I don't think "Pa" is a valid value for funct. There's no check on funct
  # values, so a user can pass it and end up with a weird error message, like 
  # the one above, unrelated to the actual problem. 
  #
  # If funct = "p" it gives a different error:
  # error in evaluating the argument 'fromCoords' in selecting a method for function 'accCost': 
  #    unable to find an inherited method for function ‘coordinates’ for signature ‘"sf"’
  # The first version of the error (with funct = "Pa"), probably means that 
  # the Conductance variable never got initialized. The second seems to 
  # indicate Conductance does exist, but it doesn't like whatever is getting 
  # passed as sp::coordinates(origin). The 'origin' object is the unmodified
  # thing that is passed to the movecost 'origin' parameter. In this case, it 
  # is sf::st_as_sf(origin_i), which is a sf data frame. From gdistance::accCost
  # documentation, it needs to be a SpatialPoints, matrix or numeric. As near 
  # as I can tell, origin_i is a single point. Maybe we can just pull it out as 
  # a matrix or two-element numeric vector?
  # origin_coords <- geom(origin_i)[, c("x", "y")]
  # dest1_coords <- geom(destination1)[, c("x", "y")]
  # No love there.
  # error in evaluating the argument 'fromCoords' in selecting a method for function 'accCost': 
  #    unable to find an inherited method for function ‘coordinates’ for signature ‘"numeric"’
  # Maybe a matrix?
  origin_coords <- matrix(geom(origin_i)[, c("x", "y")], nrow = 1, ncol = 2)
  dest1_coords <- matrix(geom(destination1)[, c("x", "y")], nrow = 1, ncol = 2)
  
  p_lcp <- movecost::movecost(
    dtm    = dem_raster,
    origin = origin_coords,
    destin = dest1_coords,
    funct  = "p", # Was "Pa"
    W      = 63.5,
    L      = 20,
    V      = 0.35,
    N      = 1.2
  )
  # )$LCPs
  # Wow.
  # Error in seq.default(min(accum_final[][is.finite(accum_final[])]), max(accum_final[]
  #    [is.finite(accum_final[])]),  : 
  #      'from' must be a finite number
  # Sending it a matrix gives an even MOAR esoteric error message. Probably 
  # good to skip movecost approach for now.

}
  

