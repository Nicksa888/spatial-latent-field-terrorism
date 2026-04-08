# ============================================================
# 1) HELPER FUNCTION TO CLEAN TEXT
# ============================================================
clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = " ")
  x[is.na(x)] <- ""
  x <- tolower(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "[[:punct:]]", " ")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x
}

# ============================================================
# 5) FUNCTION TO LOAD GADM SHAPEFILES
# ============================================================
load_gadm_country <- function(country_name, gadm_code, shape_base_dir) {
  
  folder <- file.path(shape_base_dir, paste0("gadm41_", gadm_code, "_shp"))
  shp0   <- file.path(folder, paste0("gadm41_", gadm_code, "_0.shp"))
  shp1   <- file.path(folder, paste0("gadm41_", gadm_code, "_1.shp"))
  
  if (!file.exists(shp0)) {
    stop("Missing admin0 shapefile: ", shp0)
  }
  
  if (!file.exists(shp1)) {
    stop("Missing admin1 shapefile: ", shp1)
  }
  
  adm0 <- sf::st_read(shp0, quiet = TRUE) %>%
    mutate(
      Country      = country_name,
      Country_clean = clean_text(country_name),
      gadm_code    = gadm_code
    )
  
  adm1 <- sf::st_read(shp1, quiet = TRUE) %>%
    mutate(
      Country       = country_name,
      Country_clean = clean_text(country_name),
      gadm_code     = gadm_code,
      GADM_NAME_1   = NAME_1,
      GADM_NAME_1_clean = clean_text(NAME_1)
    )
  
  list(adm0 = adm0, adm1 = adm1)
}


#######################
# add_border_distance #
#######################

# Purpose:
# Computes distance from each event to a national border or internal
# administrative border.

# Explanation:
# This function measures how far attacks are from borders. That can mean
# the outer national border or internal boundaries such as provinces,
# depending on the settings.

# What the function does:
# 1) Reads the border shapefile.
# 2) Converts both points and borders into a metric projected CRS.
# 3) Builds either:
#    - the outer national outline, or
#    - internal administrative borders.
# 4) Computes point-to-border distances in chunks.
# 5) Writes back:
#    - distance in meters,
#    - distance in kilometres,
#    - distance bins.

# Output:
# Returns the original dataset with border-distance variables added.

# Intended use:
# Used to operationalise border proximity as an operational constraint.

add_border_distance <- function(
    data,
    border_path,
    lon_col = "Longitude",
    lat_col = "Latitude",
    crs_points = 4326,
    dist_prefix = "B",
    breaks_km = c(0, 10, 50, 100, Inf),
    labels = c("0–10 km", "10–50 km", "50–100 km", "100+ km"),
    save_rds_path = NULL,
    quiet = TRUE,
    utm_epsg = 32718,        # UTM 18S (meters)
    chunk_size = 2000,
    use_s2 = FALSE,
    internal_only = FALSE    # TRUE => internal admin borders (e.g province borders)
) {
  stopifnot(file.exists(border_path))
  stopifnot(all(c(lon_col, lat_col) %in% names(data)))
  
  sf::sf_use_s2(use_s2)
  
  make_valid_any <- function(x) {
    v <- sf::st_is_valid(x)
    if (all(v, na.rm = TRUE)) return(x)
    
    if ("st_make_valid" %in% getNamespaceExports("sf")) {
      return(sf::st_make_valid(x))
    }
    
    if (requireNamespace("lwgeom", quietly = TRUE) &&
        "st_make_valid" %in% getNamespaceExports("lwgeom")) {
      return(lwgeom::st_make_valid(x))
    }
    
    suppressWarnings(sf::st_buffer(x, 0))
  }
  
  # Points -> sf -> projected #
  
  pts <- dplyr::as_tibble(data) %>%
    dplyr::mutate(
      .lon = as.numeric(.data[[lon_col]]),
      .lat = as.numeric(.data[[lat_col]])
    )
  
  ok <- is.finite(pts$.lon) & is.finite(pts$.lat)
  pts_ok <- pts[ok, , drop = FALSE]
  
  pts_sf <- sf::st_as_sf(pts_ok, coords = c(".lon", ".lat"), crs = crs_points, remove = FALSE)
  pts_sf <- sf::st_transform(pts_sf, utm_epsg)
  
  # Borders -> sf -> projected
  
  border <- sf::st_read(border_path, quiet = quiet)
  if (is.na(sf::st_crs(border))) {
    stop("Border shapefile has NA CRS. Define it (st_set_crs) before transforming.")
  }
  border <- sf::st_transform(border, utm_epsg)
  border <- make_valid_any(border)
  
  geom <- sf::st_geometry(border)
  
  # Build target line geometry for distance calculation
  # All polygon boundaries (includes internal + outer)
  all_lines <- sf::st_boundary(geom)
  all_lines_union <- sf::st_union(all_lines)
  
  if (internal_only) {
    # Outer national outline
    nat_outline <- sf::st_boundary(sf::st_union(geom))
    
    # Remove the outer outline, keep internal boundaries
    border_line <- tryCatch(
      sf::st_difference(all_lines_union, nat_outline),
      error = function(e) {
        # Fallback (still gives something sensible)
        all_lines_union
      }
    )
  } else {
    # National border (outer outline only)
    border_line <- sf::st_boundary(sf::st_union(geom))
  }
  
  # Ensure line type
  border_line <- sf::st_cast(border_line, "MULTILINESTRING", warn = FALSE)
  
  # Distances in chunks
  
  n_ok <- nrow(pts_sf)
  d_m_ok <- numeric(n_ok)
  
  idx <- split(seq_len(n_ok), ceiling(seq_len(n_ok) / chunk_size))
  
  for (ii in seq_along(idx)) {
    rows <- idx[[ii]]
    d <- sf::st_distance(pts_sf[rows, ], border_line)  # units matrix
    
    # If result is a matrix with multiple columns, take min per row
    if (is.matrix(d) && ncol(d) > 1) {
      d_num <- apply(d, 1, min)
    } else {
      d_num <- as.numeric(d)
    }
    d_m_ok[rows] <- as.numeric(d_num)
  }
  
  
  # Write back to full data
  
  out <- data
  
  dist_m_col  <- paste0(dist_prefix, "_Dist")
  dist_km_col <- paste0(dist_prefix, "_Dist_km")
  dist_bin_col <- paste0(dist_prefix, "_Dist_bin")
  
  out[[dist_m_col]]  <- NA_real_
  out[[dist_km_col]] <- NA_real_
  
  out[[dist_m_col]][ok]  <- d_m_ok
  out[[dist_km_col]][ok] <- d_m_ok / 1000
  
  out[[dist_bin_col]] <- cut(
    out[[dist_km_col]],
    breaks = breaks_km,
    labels = labels,
    include.lowest = TRUE,
    right = FALSE
  )
  
  if (!is.null(save_rds_path)) saveRDS(out, save_rds_path)
  
  out
}

# Example Usage #

# data_border <- add_border_distance(
#   data = data_kernel,
#   border_path = "Shapefiles/gadm41_country.shp",
#   lon_col = "Longitude",
#   lat_col = "Latitude",
#   dist_prefix = "B",
#   utm_epsg = 32633
# )
