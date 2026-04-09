################################################################################
################# MODERNIZED SCRIPT FOR FIELD ET AL. 2021 ######################
################################################################################

library(terra)
library(sf)
#library(movecost)
library(gdistance)
library(leastcostpath)
library(tidyverse)
library(spatialEco)
#library(magrittr)
# library(sp)
# library(devtools)

# set working directory and set master projection if re-projecting
setwd("C:/Users/Colin Omilanowski/Downloads/Field_Test/Arch_test")

# Set global projection
master_crs <- "EPSG:26913" # UTM Zone 13N NAD83

# Whether or not to re-calculate surfaces. If TRUE, will re-calculate all 
# conductance surfaces, regardless of whether or not they are already on disk.
# recalc_surfaces should only really need to be set to TRUE if the input DEM 
# has changed or if the number of neighbors changed.
recalc_surfaces <- FALSE

###########################################################
################## STEP 1: IMPORT DATA ####################
###########################################################

# ── FIX #1: Create output directory if it doesn't exist ──────────────────────
# Added here at the top so the directory is guaranteed to exist before
# any st_write() calls fire inside the loop.
if (!dir.exists("./OUTPUT")) dir.create("./OUTPUT", recursive = TRUE)

# Create folder to hold conductance surfaces (if not already on computer)
if (!dir.exists("./OUTPUT/surfaces")) {
  dir.create("./OUTPUT/surfaces")
}

# Import vectors using sf
# Ensure these files exist in your directory
origins <- terra::vect("./raw_data/origin_locations.shp") 
sites <- terra::vect("./raw_data/Sites_Test.shp")

destination1 <- sites[sites$Name == "Corridor 2 destination", ]
destination2 <- sites[sites$Name == "Corridor 3 destination", ]

# Import DEM using terra
dem <- terra::rast('./raw_data/DEM_corridor1.tif')
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
neigh <- 16

# For each of the conductance surfaces, we do not necessarily need to calculate 
# them again, if they are on disk. We check to see if a saved version exists, 
# and if so, we load it in. If it does not exist on disk, the surface is 
# calculated and saved to disk.

# Tobler surface
if (recalc_surfaces || !file.exists("./OUTPUT/surfaces/tobler_cs.rds")) {
  message("Calculating Tobler surface")
  tobler_cs   <- leastcostpath::create_slope_cs(dem, cost_function = "tobler", neighbours = neigh)
  saveRDS(object = tobler_cs, file = "./OUTPUT/surfaces/tobler_cs.rds")
} else {
  tobler_cs <- readRDS(file = "./OUTPUT/surfaces/tobler_cs.rds")
  message("Tobler surface loaded from disk")
}

# Irmischer-Clarke offpath male surface
if (recalc_surfaces || !file.exists("./OUTPUT/surfaces/ic_off_m_cs.rds")) {
  message("Calculating Irmischer-Clarke offpath male surface")
  ic_off_m_cs <- leastcostpath::create_slope_cs(dem, cost_function = "irmischer-clarke offpath male", neighbours = neigh)
  saveRDS(object = ic_off_m_cs, file = "./OUTPUT/surfaces/ic_off_m_cs.rds")
} else {
  ic_off_m_cs <- readRDS(file = "./OUTPUT/surfaces/ic_off_m_cs.rds")
  message("Irmischer-Clarke offpath male surface loaded from disk")
}

# Campbell 2019 surface
if (recalc_surfaces || !file.exists("./OUTPUT/surfaces/campbell_cs.rds")) {
  message("Calculating Campbell 2019 surface")
  campbell_cs <- leastcostpath::create_slope_cs(dem, cost_function = "campbell 2019 50", neighbours = neigh)
  saveRDS(object = campbell_cs, file = "./OUTPUT/surfaces/campbell_cs.rds")
} else {
  campbell_cs <- readRDS(file = "./OUTPUT/surfaces/campbell_cs.rds")
  message("Campbell 2019 surface loaded from disk")
}

# Herzog surface
if (recalc_surfaces || !file.exists("./OUTPUT/surfaces/herzog_cs.rds")) {
  message("Calculating Herzog surface")
  herzog_cs <- leastcostpath::create_slope_cs(dem, cost_function = "campbell 2019 50", neighbours = neigh)
  saveRDS(object = herzog_cs, file = "./OUTPUT/surfaces/herzog_cs.rds")
} else {
  herzog_cs <- readRDS(file = "./OUTPUT/surfaces/herzog_cs.rds")
  message("Herzog surface loaded from disk")
}

# Llobera-Sluckin surface
if (recalc_surfaces || !file.exists("./OUTPUT/surfaces/llobera-sluckin_cs.rds")) {
  message("Calculating Llobera-Sluckin surface")
  ls_cs <- leastcostpath::create_slope_cs(dem, cost_function = "llobera-sluckin", neighbours = neigh)
  saveRDS(object = ls_cs, file = "./OUTPUT/surfaces/llobera-sluckin_cs.rds")
} else {
  ls_cs <- readRDS(file = "./OUTPUT/surfaces/llobera-sluckin_cs.rds")
  message("Llobera-Sluckin surface loaded from disk")
}

################################################################################
# Pandolf cost surface requires a lot more gymnastics
if (recalc_surfaces || !file.exists("./OUTPUT/surfaces/pandolf_cs.rds")) {
  message("Calculating Pandolf surface")
  
  # Adapting from https://github.com/sfield2/Energetics-in-Least-Cost-Analysis
  # See the file ANALYSIS.R
  
  # Define inputs for Pandolf et al. (1977) function
  W <- 63.5 # body mass
  L <- 20   # load mass (kg)
  V <- 0.35 # velocity (m/s)
  N <- 1.2  # terrain
  
  # Create conductance surface using Van Etten (2017) methodology
  # First, create function to calculate the attitudinal difference between 
  # adjacent cells
  altDiff <- function(x) {
    return(x[2] - x[1])
  }
  
  # We need a RasterLayer version of the DEM
  dem_raster <- raster::raster(dem)
  
  # Apply the altidutinal distance function to the (RasterLayer) DEM
  hd <- gdistance::transition(x = dem_raster, 
                              transitionFunction = altDiff, 
                              directions = neigh, 
                              symm = FALSE)
  # gdistance::transition above returns:
  # In .TfromR_old(x, transitionFunction, directions, symm) :
  #   transition function gives negative values
  # But we are going to just hum through the warning for now :-\
    
  # Use the geoCorrection function to divide the attitudinal difference by the 
  # distance between cells (i.e., calculating slope as rise over run)
  slope <- gdistance::geoCorrection(x = hd)
  
  # Values between non-adjacent cells is 0, but the slope between these cells is 
  # not 0! so therefore, we need to restrict the calculation to adjacent cells. 
  # This is done by creating an index for adjacent cells (`adj`) with the 
  # function `adjacent` (Van Etten 2017)
  adj <- raster::adjacent(x = dem_raster, 
                          cells = 1:ncell(dem_raster), 
                          pairs = TRUE, 
                          directions = neigh)
  
  # Replicate slope for duplicated use in following steps
  cost <- slope
  
  # Create Pandolf et al. (1977) function and apply this function to the entire 
  # surface  to create cost function 
  cost_function_pan <- function(x) { 
    return(1 / (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * (abs(x[adj])*100)))) 
  }
  cost[adj] <- cost_function_pan(x = slope)
  Conductance_1 <- gdistance::geoCorrection(x = cost)
  
  # NOTE: Importantly, Pandolf et al. (1977) function creates negative values 
  # on negative slopes, so the function needs to be applied differently 
  # depending on the slope (Herzog 2010, 2014; White 2012; White and Barber 
  # 2012).Therefore two conductance surfaces are going to be created, delimited 
  # depending on the slope they represent, and then combined 
  
  # Create a second surface using Santee et al. (2001) function, which fixes 
  # issue of using Pandolf et al. (1977) over negative slopes 
  cost_function_santee <- function(x) { 
    return(1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * (abs(x[adj])*100)))- N * ((V * (abs(x[adj]) * 100)* (W + L) / 3.5) - ((W + L) * ((abs(x[adj])*100) + 6)^2 / W) + (25-(V^2)))))
  }
  cost[adj] <- cost_function_santee(x = slope)
  Conductance_2 <- gdistance::geoCorrection(x = cost)
  
  # Determine which portions of DEM are positive and negative slope by 
  # classifying rasters. This will produce two rasters where cells equal 1 
  # (positive slopes) and 0 (negative slopes), OR vice versa.
  # NOTE: ensure "0" of classified rasters equal N/A
  slope_r <- raster::raster(x = slope)
  slope_r[is.na(slope_r)] <- raster::maxValue(x = slope_r)
  positiveslope <- slope_r >=0
  negativeslope <- slope_r < 0
  negativeslope[negativeslope == 0] <- NA
  positiveslope[positiveslope==0] <- NA
  
  # Create barrier out of positive and negative raster using function from 
  # earlier version of leastcostpath package (Lewis 2021). 
  
  # Below is a stripped-down version of create_barrier_cs() function, which was 
  # removed from the leastcostpath package around 2023. It can be found in the 
  # repo's history at 
  # https://github.com/josephlewis/leastcostpath/blob/d3f322b7c328a23caa1e7e6905425998704d7664/R/create_barrier_cs.R
  
  #' @param raster \code{RasterLayer} (raster package). The Resolution, Extent, 
  #' and Spatial Reference System of the provided RasterLayer is used when 
  #' creating the resultant Barrier Cost Surface
  #' @param barrier \code{RasterLayer} (raster package). Area within the 
  #' landscape that movement is inhibited
  #' @param neighbours \code{numeric} value. Number of directions used in the 
  #' Least Cost Path calculation. See Huber and Church (1985) for methodological 
  #' considerations when choosing number of neighbours. Expected numeric values 
  #' are 4, 8, 16, 32, 48 or a matrix object. Default is numeric value 16
  #' @param field \code{numeric} value. Value assigned to cells that coincide 
  #' with the barrier. Default is \code{numeric value 0}
  #' @param background \code{numeric} value. Value assigned to cells that do not 
  #' coincide with the barrier. Default is \code{numeric value 1}
  
  create_barrier_cs <- function(raster, barrier, neighbours = 16, field = 0, 
                                background = 1) {
    # A little defensive programming making sure the value passed to the 
    # neighbors argument is OK.
    if (any(!neighbours %in% c(4, 8, 16, 32, 48)) & (!inherits(neighbours, "matrix"))) {
      stop("neighbours argument is invalid. Expecting 4, 8, 16, 32, 48, or matrix object")
    }
    
    
    # create TransitionLayer of zeroes based on raster dimensions
    barrier_cs <- new("TransitionLayer", nrows = as.integer(nrow(raster)), 
                      ncols = as.integer(ncol(raster)), 
                      extent = raster::extent(raster), 
                      crs = raster::projection(raster, asText = FALSE), 
                      transitionMatrix = Matrix::Matrix(0, raster::ncell(raster), 
                                                        raster::ncell(raster)), 
                      transitionCells = 1:raster::ncell(raster))
    
    barrier_cells <- which(!is.na(raster::getValues(barrier)))
    
    # get adjacent cells to limit value change
    adj <- raster::adjacent(raster, cells = 1:raster::ncell(raster), 
                            pairs = TRUE, directions = neighbours)
    
    # change values that coincide with barrier
    barrier_cs[adj[adj[, 2] %in% barrier_cells, ]] <- field
    
    # change values that don't coincide with barrier
    barrier_cs[adj[!adj[, 2] %in% barrier_cells, ]] <- background
    
    # drop zeroes from matrix
    barrier_cs@transitionMatrix <- Matrix::drop0(barrier_cs@transitionMatrix)
    
    return(barrier_cs)
    
  }
  
  positiveslopesonly <- create_barrier_cs(raster = dem_raster, 
                                          barrier = negativeslope, 
                                          neighbours = neigh, 
                                          field = 0, 
                                          background = 1)
  
  negativeslopeonly <- create_barrier_cs(raster = dem_raster, 
                                         barrier = positiveslope, 
                                         neighbours = neigh, 
                                         field = 0, 
                                         background = 1)
  
  # Multiply Pandolf et al. (1977) derived conductance surface by barrier 
  # raster and multiply Santee et al. (2001) derived conductance surface by 
  # other barrier raster. 
  pandolfpositive_conductance <- Conductance_1*positiveslopesonly
  pandolfnegative_conductance <- Conductance_2*negativeslopeonly
  
  # Merge the two delimited conductance surfaces 
  pandolf_cs <- pandolfpositive_conductance+pandolfnegative_conductance
  
  # We aren't done yet! At this point, we have a TransitionLayer (from 
  # gdistance), but we are going to need a conductanceMatrix (like in 
  # leastcostpath). So here is some ugly
  pandolf_cs <- list(conductanceMatrix = pandolf_cs@transitionMatrix, 
                         costFunction = NA, 
                         "maxSlope" = NA, 
                         exaggeration = FALSE, 
                         criticalSlope = NA,
                         neighbours = sum(neigh, na.rm = TRUE), 
                         nrow = terra::nrow(dem), 
                         ncol = terra::ncol(dem), 
                         "resolution" = terra::res(dem), 
                         "extent" = as.vector(terra::ext(dem)),  
                         crs = terra::crs(dem, proj = TRUE))
  class(pandolf_cs) <- "conductanceMatrix"
  
    # There are some large-ish temporary objects in memory, so we will remove them
  rm(dem_raster, hd, slope, cost, adj, Conductance_1, Conductance_2, slope_r,
     positiveslope, negativeslope, positiveslopesonly, negativeslopeonly,
     pandolfpositive_conductance, pandolfnegative_conductance)
  # And explicitly clean up the garbage
  invisible(gc())

  # After all that work, save the cost surface to disk!  
  saveRDS(object = pandolf_cs, file = "./OUTPUT/surfaces/pandolf_cs.rds")
} else {
  pandolf_cs <- readRDS(file = "./OUTPUT/surfaces/pandolf_cs.rds")
  message("Pandolf surface loaded from disk")
}

###########################################################
########### STEP #3: CALCULATE LEAST COST PATHS ###########
########### FROM SERIES OF ORIGIN POINTS TO  ##############
########### SEPARATE DESTINATIONS AND COMPUTE #############
########### DETERMINING VARIABLES (TERRAIN, DISTANCE) #####
###########################################################

n <- nrow(origins)   # terra::nrow() or length() both work on SpatVector

# Build results table
result_destination1 <- as.data.frame(matrix(NA, nrow = n, ncol = 9))
colnames(result_destination1) <- c(
  "time (sq.m)", "energy (sq.m)", "overlap (sq.m)",
  "mean time and energy (sq.m)", "overlap as percent of average",
  "mean TRI within 1km of time lcps",    "median TRI within 1km of time lcps",
  "mean TRI within 1km of energy lcps",  "median TRI within 1km of energy lcps"
)

for (i in seq_len(n)){ 
  origin_i <- origins[i, ]
  message("calculating least cost path for origin " , i)
  # ── Calculate LCPs ──────────────────────────────────────────────────────
  # leastcostpath 2.x: create_lcp() signature is largely unchanged but
  # accepts SpatVector (terra) instead of Spatial* objects.
  
  # All the invisible(gc()) calls are attempts to help with garbage collection 
  # so RAM doesn't cook
  t_lcp  <- leastcostpath::create_lcp(tobler_cs,   origin_i, destination1)
  invisible(gc)
  ic_lcp <- leastcostpath::create_lcp(ic_off_m_cs, origin_i, destination1)
  invisible(gc)
  c_lcp  <- leastcostpath::create_lcp(campbell_cs, origin_i, destination1)
  invisible(gc)
  h_lcp  <- leastcostpath::create_lcp(herzog_cs,   origin_i, destination1)
  invisible(gc)
  ls_lcp <- leastcostpath::create_lcp(ls_cs,       origin_i, destination1)
  invisible(gc)
  p_lcp <-  leastcostpath::create_lcp(pandolf_cs,  origin_i, destination1)
  invisible(gc)
  
  # terra::crs(tobler_cs, describe = TRUE)
  # terra::crs(dem, describe = TRUE)
  # terra::crs(origins, describe = TRUE)
  # terra::crs(destination1, describe = TRUE)
  # 
  # terra::ext(dem)
  # terra::geom(origin_i)
  # terra::geom(destination1)
  # 
  # class(tobler_cs)
  # print(tobler_cs)
  
  
  # Replace the pandolf_cs create_lcp line with this:
  # Then use dem_raster instead of dem in the movecost call
  # p_lcp <- movecost::movecost(
  #   dtm    = dem_raster,
  #   origin = sf::st_as_sf(origin_i),
  #   destin = sf::st_as_sf(destination1),
  #   funct  = "Pa",
  #   W      = 63.5,
  #   L      = 20,
  #   V      = 0.35,
  #   N      = 1.2
  # )$LCPs
  
  
  # Convert lcp SpatVector lines to sf for geometry operations
  # sf is more feature-rich for line/polygon operations than terra
  t_lcp_sf  <- sf::st_as_sf(t_lcp)
  ic_lcp_sf <- sf::st_as_sf(ic_lcp)
  c_lcp_sf  <- sf::st_as_sf(c_lcp)
  h_lcp_sf  <- sf::st_as_sf(h_lcp)
  ls_lcp_sf <- sf::st_as_sf(ls_lcp)
  # p_lcp_sf  <- sf::st_as_sf(p_lcp)
  dest1_sf  <- sf::st_as_sf(destination1)
  
  # ── Sample points along each time-based LCP at ~100 m intervals ─────────
  # sf::st_line_sample() replaces sp::spsample()
  sample_pts <- function(lcp_sf, interval = 100) {
    len_m <- as.numeric(sf::st_length(lcp_sf))
    n_pts <- max(2L, round(len_m / interval))
    sf::st_line_sample(lcp_sf, n = n_pts, type = "regular") |>
      sf::st_cast("POINT") |>
      sf::st_as_sf()
  }
  
  t_p  <- sample_pts(t_lcp_sf)
  ic_p <- sample_pts(ic_lcp_sf)
  c_p  <- sample_pts(c_lcp_sf)
 
  
  # Harmonise to minimum point count across the three sets
  l <- min(nrow(t_p), nrow(ic_p), nrow(c_p))
  t_p  <- t_p[seq_len(l), ]
  ic_p <- ic_p[seq_len(l), ]
  c_p  <- c_p[seq_len(l), ]
  
  # Helper: nearest distance from a point to a line (sf)
  dist_pt_to_line <- function(pt_sf, line_sf) {
    as.numeric(sf::st_distance(pt_sf, line_sf))
  }
  
  # Helper: distance from a point to the destination
  dist_to_dest <- function(pt_sf, dest_sf, lcp_sf) {
    lcp_len <- as.numeric(sf::st_length(lcp_sf))
    d       <- as.numeric(sf::st_distance(pt_sf, dest_sf))
    abs(d - lcp_len) / lcp_len
  }
  
  # ── Build per-interval comparison table ─────────────────────────────────
  # result_difference <- as.data.frame(matrix(NA, nrow = l, ncol = 18))
  # colnames(result_difference) <- c(
  #   "tobler_pandolf",        "tobler_lloberasluckin",   "tobler_herzog",
  #   "irmischer_pandolf",     "irmischer_herzog",         "irmischer_lloberasluckin",
  #   "campbell_pandolf",      "campbell_lloberasluckin",  "campbell_herzog",
  #   "tp_pctdist",            "tls_pctdist",              "th_pctdist",
  #   "ip_pctdist",            "ils_pctdist",              "ih_pctdist",
  #   "cp_pctdist",            "cls_pctdist",              "ch_pctdist"
  # )
  # 
  # for (j in seq_len(l)) {
  #   result_difference[j,  1] <- dist_pt_to_line(t_p[j,],  p_lcp_sf)
  #   result_difference[j,  2] <- dist_pt_to_line(t_p[j,],  ls_lcp_sf)
  #   result_difference[j,  3] <- dist_pt_to_line(t_p[j,],  h_lcp_sf)
  #   result_difference[j,  4] <- dist_pt_to_line(ic_p[j,], p_lcp_sf)
  #   result_difference[j,  5] <- dist_pt_to_line(ic_p[j,], ls_lcp_sf)
  #   result_difference[j,  6] <- dist_pt_to_line(ic_p[j,], h_lcp_sf)
  #   result_difference[j,  7] <- dist_pt_to_line(c_p[j,],  p_lcp_sf)
  #   result_difference[j,  8] <- dist_pt_to_line(c_p[j,],  ls_lcp_sf)
  #   result_difference[j,  9] <- dist_pt_to_line(c_p[j,],  h_lcp_sf)
  #   result_difference[j, 10:12] <- dist_to_dest(t_p[j,],  dest1_sf, t_lcp_sf)
  #   result_difference[j, 13:15] <- dist_to_dest(ic_p[j,], dest1_sf, ic_lcp_sf)
  #   result_difference[j, 16:18] <- dist_to_dest(c_p[j,],  dest1_sf, c_lcp_sf)
  # }
  # Abbreviated version excluding pandolf
  result_difference <- as.data.frame(matrix(NA, nrow = l, ncol = 18))
  colnames(result_difference) <- c(
    "tobler_pandolf",        "tobler_lloberasluckin",   "tobler_herzog",
    "irmischer_pandolf",     "irmischer_herzog",         "irmischer_lloberasluckin",
    "campbell_pandolf",      "campbell_lloberasluckin",  "campbell_herzog",
    "tp_pctdist",            "tls_pctdist",              "th_pctdist",
    "ip_pctdist",            "ils_pctdist",              "ih_pctdist",
    "cp_pctdist",            "cls_pctdist",              "ch_pctdist"
  )
  
  for (j in seq_len(l)) {
    result_difference[j,  1] <- NA # dist_pt_to_line(t_p[j,],  p_lcp_sf)
    result_difference[j,  2] <- dist_pt_to_line(t_p[j,],  ls_lcp_sf)
    result_difference[j,  3] <- dist_pt_to_line(t_p[j,],  h_lcp_sf)
    result_difference[j,  4] <- NA # dist_pt_to_line(ic_p[j,], p_lcp_sf)
    result_difference[j,  5] <- dist_pt_to_line(ic_p[j,], ls_lcp_sf)
    result_difference[j,  6] <- dist_pt_to_line(ic_p[j,], h_lcp_sf)
    result_difference[j,  7] <-NA # dist_pt_to_line(c_p[j,],  p_lcp_sf)
    result_difference[j,  8] <- dist_pt_to_line(c_p[j,],  ls_lcp_sf)
    result_difference[j,  9] <- dist_pt_to_line(c_p[j,],  h_lcp_sf)
    result_difference[j, 10:12] <- dist_to_dest(t_p[j,],  dest1_sf, t_lcp_sf)
    result_difference[j, 13:15] <- dist_to_dest(ic_p[j,], dest1_sf, ic_lcp_sf)
    result_difference[j, 16:18] <- dist_to_dest(c_p[j,],  dest1_sf, c_lcp_sf)
  }
  write.csv(result_difference,
            paste0("./OUTPUT/result_difference_1_", i, ".csv"),
            row.names = FALSE)
  
  # ── Bind & buffer time/energy LCPs ──────────────────────────────────────
  # sf::st_union() + sf::st_buffer() replaces gBuffer(rbind(...))
  t_all_sf <- sf::st_union(rbind(t_lcp_sf, ic_lcp_sf, c_lcp_sf))
  # e_all_sf <- sf::st_union(rbind(h_lcp_sf, ls_lcp_sf, p_lcp_sf))
  e_all_sf <- sf::st_union(rbind(h_lcp_sf, ls_lcp_sf))
  
  t_buff_sf <- sf::st_buffer(t_all_sf, dist = 125)
  e_buff_sf <- sf::st_buffer(e_all_sf, dist = 125)
  
  # Export buffered corridors as shapefiles
  # sf::st_write() replaces rgdal::writeOGR()
  sf::st_write(sf::st_as_sf(t_buff_sf),
               paste0("./OUTPUT/all_time_buff_origin",   i, ".shp"),
               delete_dsn = TRUE)
  sf::st_write(sf::st_as_sf(e_buff_sf),
               paste0("./OUTPUT/all_energy_buff_origin", i, ".shp"),
               delete_dsn = TRUE)
  
  # ── Overlap between time & energy corridors ──────────────────────────────
  # sf::st_intersection() replaces rgeos::gIntersection()
  overlap_sf <- sf::st_intersection(t_buff_sf, e_buff_sf)
  #extract polygons from intersection
  if(sf::st_geometry_type(overlap_sf) == "GEOMETRYCOLLECTION") {
    overlap_sf <- sf::st_collection_extract(x=overlap_sf,type = "POLYGON")
  }
  
  sf::st_write(sf::st_as_sf(overlap_sf),
               paste0("./OUTPUT/overlap_buff_origin", i, ".shp"),
               delete_dsn = TRUE)
  
  # ── TRI extraction within buffered corridors ─────────────────────────────
  # terra::extract() replaces raster::extract()
  t_buff_terra  <- terra::vect(t_buff_sf)
  e_buff_terra  <- terra::vect(e_buff_sf)
  
  t_buff_tri_extract <- terra::extract(tri_dem, t_buff_terra)$TRI
  e_buff_tri_extract <- terra::extract(tri_dem, e_buff_terra)$TRI
  
  # Remove NAs before summary stats
  t_buff_tri_extract <- t_buff_tri_extract[!is.na(t_buff_tri_extract)]
  e_buff_tri_extract <- e_buff_tri_extract[!is.na(e_buff_tri_extract)]
  
  # Bind and export TRI table
  all_tri <- rbind(
    data.frame(value = t_buff_tri_extract, cost = "time"),
    data.frame(value = e_buff_tri_extract, cost = "energy")
  )
  write.csv(all_tri,
            paste0("./OUTPUT/result_tri_1_", i, ".csv"),
            row.names = FALSE)
  
  # ── Summary statistics for this loop iteration ───────────────────────────
  # terra::expanse() replaces raster::area() for polygons
  result_destination1[i, 1] <- terra::expanse(t_buff_terra)
  result_destination1[i, 2] <- terra::expanse(e_buff_terra)
  result_destination1[i, 3] <- terra::expanse(terra::vect(overlap_sf))
  result_destination1[i, 4] <- (terra::expanse(t_buff_terra) +
                                  terra::expanse(e_buff_terra)) / 2
  result_destination1[i, 5] <- terra::expanse(terra::vect(overlap_sf)) /
    result_destination1[i, 4]
  result_destination1[i, 6] <- mean(t_buff_tri_extract)
  result_destination1[i, 7] <- median(t_buff_tri_extract)
  result_destination1[i, 8] <- mean(e_buff_tri_extract)
  result_destination1[i, 9] <- median(e_buff_tri_extract)
}
write.csv(result_destination1,
          "./OUTPUT/results_destination_1.csv",
          row.names = FALSE)


# 1. Convert the terra DEM to the older raster format
# This is required because movecost's internal 'transition' function 
# does not yet support 'SpatRaster' objects directly in all versions.
dem_raster <- raster::raster(dem)

# 2. Run movecost with the compatible object
# Use 't' instead of 'N' (movecost uses 't' for terrain weight)
# Put "h" in quotes
pandolf_data <- movecost(dtm = dem_raster, origin = as(origins[1,], "Spatial"),
                         funct = "pandolf", W = W, L = L, V = V, t = N, time = "h")

# The transition layer is stored here:
pandolf_cs <- pandolf_data$transitionLayer  # correct slot name

# Extract the conductance matrix
pandolf_cs <- pandolf_data$conductance

###########################################################
########### STEP 3: LOOPED LCP ANALYSIS ###################
###########################################################

n <- nrow(origins)
result_destination1 <- data.frame(matrix(NA, nrow = n, ncol = 9))
colnames(result_destination1) <- c("time_area", "energy_area", "overlap_area", "mean_area", 
                                   "overlap_pct", "mean_tri_time", "med_tri_time", 
                                   "mean_tri_energy", "med_tri_energy")

for(i in 1:n) {
  origin <- origins[i,]
  
  # Calculate LCPs
  t_lcp  <- create_lcp(tobler_cs, origin, destination1)
  ic_lcp <- create_lcp(ic_off_m_cs, origin, destination1)
  c_lcp  <- create_lcp(campbell_cs, origin, destination1)
  h_lcp  <- create_lcp(herzog_cs, origin, destination1)
  ls_lcp <- create_lcp(ls_cs, origin, destination1)
  p_lcp  <- create_lcp(pandolf_cs, origin, destination1)
  
  # Buffer paths using sf (125m buffer)
  # Combine time-based paths
  t_all <- st_union(st_as_sf(rbind(t_lcp, ic_lcp, c_lcp)))
  t_buff <- st_buffer(t_all, 125)
  
  # Combine energy-based paths
  e_all <- st_union(st_as_sf(rbind(h_lcp, ls_lcp, p_lcp)))
  e_buff <- st_buffer(e_all, 125)
  
  # Intersection/Overlap
  overlap <- st_intersection(t_buff, e_buff)
  
  # Terrain Extraction (terra)
  t_tri_vals <- terra::extract(tri_dem, vect(t_buff))[[2]]
  e_tri_vals <- terra::extract(tri_dem, vect(e_buff))[[2]]
  
  # Calculations
  t_area <- as.numeric(st_area(t_buff))
  e_area <- as.numeric(st_area(e_buff))
  o_area <- if(length(overlap) > 0) as.numeric(st_area(overlap)) else 0
  
  result_destination1[i, 1:5] <- c(t_area, e_area, o_area, (t_area+e_area)/2, o_area/((t_area+e_area)/2))
  result_destination1[i, 6:9] <- c(mean(t_tri_vals, na.rm=T), median(t_tri_vals, na.rm=T),
                                   mean(e_tri_vals, na.rm=T), median(e_tri_vals, na.rm=T))
  
  # Export Spatial Data
  st_write(t_buff, paste0("./OUTPUT/time_buff_", i, ".shp"), delete_dsn = TRUE)
}

write.csv(result_destination1, "./OUTPUT/results_summary.csv")

###########################################################
########### STEP 4: ROADS AS CONDUITS (TERRA) #############
###########################################################

# Use terra::resample and logic for conduits
conduit_rast <- rast('./DATA/DEM/DEM_corridor3_conduit.tif')
dem3 <- rast('./DATA/DEM/DEM_corridor3.tif')

conduit_resampled <- resample(conduit_rast, dem3, method = "bilinear")
# Classify conduit (1 for road, NA elsewhere)
road_mask <- conduit_resampled > 1.5 & conduit_resampled < 2.0
road_mask[road_mask == 0] <- NA

# Add conductance boost
base_cs <- create_slope_cs(dem3, cost_function = "tobler", neighbours = 16)
# Create a boost layer where the road is
road_boost <- create_barrier_cs(dem3, barrier = vect(road_mask), field = 0.1, background = 0)
road_conductance <- base_cs + road_boost

###########################################################
########### STEP 5: FATIGUE THRESHOLDS ####################
###########################################################

# Note: accCost is now handled by leastcostpath::create_accumulated_cost
acc_tobler <- create_accumulated_cost(tobler_cs, origins[1,])

# Create 30-minute threshold (1800s)
fatigue_mask <- acc_tobler < 1800
# Convert to vector for intersection
fatigue_poly <- as.polygons(fatigue_mask) %>% st_as_sf() %>% filter(layer == 1)
fatigue_line <- st_cast(fatigue_poly, "LINESTRING")

# Find intersection point (replaces gIntersection)
intersection_pt <- st_intersection(st_as_sf(t_lcp), fatigue_line)

# Proceed with leg 2 calculations using this intersection_pt as the new origin.