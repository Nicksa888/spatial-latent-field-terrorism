rm(list = ls())

# ============================================================
# RESEARCH QUESTION
# Does each group operate within its own hidden spatial
# opportunity structure?
#
# STRONGER DESIGN
# - shared latent field
# - ELN-specific latent field
# - FARC-specific latent field
#
# ADDED FROM PRIOR SCRIPT
# - SPDE prior grid search
# - comparison across prior choices
# - M0 / M1 / M2 comparison using best prior
# ============================================================

# ============================================================
# Libraries
# ============================================================

library(INLA)
library(sf)
library(dplyr)
library(ggplot2)
library(Matrix)
library(purrr)
library(stringr)
library(geosphere)
library(rnaturalearth)
library(elevatr)
library(terra)
library(geodata)
library(raster)
library(spatstat.geom)
library(tibble)

# ============================================================
# Source functions
# ============================================================

source("C:/R Portfolio/Terrorism_Different_Civil_War/Latent_Field_Analyses/Latent_Analysis_Functions_07_04_2026.R")
source("C:/R Portfolio/Terrorism_Different_Civil_War/Microcycle/Microcycle Functions_24_03_2026.R")


# ============================================================
# 1) INPUT DATA
# ============================================================

GTD <- read.csv(
  "C:/R Portfolio/Terrorism_Different_Civil_War/GTD/GTD_Final_29_09_25.csv",
  stringsAsFactors = FALSE
)

GTD$Date <- as.Date(GTD$Date, format = "%d/%m/%Y")

GTD_Colombia <- GTD %>%
  filter(Country == "Colombia")
table(GTD_Colombia$Group)

# ============================================================
# 2) GEOGRAPHIC VARIABLES
# ============================================================

GTD_Colombia_xy <- GTD_Colombia[, c("Longitude", "Latitude")]

dem <- raster("C:/R Portfolio/Terrorism_Different_Civil_War/Victoria_Volodina_Supervisor_Meetings/color_etopo1_ice_full.tif")
crs(dem) <- "+proj=longlat +datum=WGS84"

slope_rast   <- terrain(dem, opt = "slope",   unit = "degrees")
aspect_rast  <- terrain(dem, opt = "aspect",  unit = "degrees")
flowdir_rast <- terrain(dem, opt = "flowdir")

# Extract topographic variables for Colombia
GTD_Colombia$Elevation <- raster::extract(dem, GTD_Colombia_xy)
GTD_Colombia$Slope     <- raster::extract(slope_rast, GTD_Colombia_xy)
GTD_Colombia$Aspect    <- raster::extract(aspect_rast, GTD_Colombia_xy)
GTD_Colombia$FlowDir   <- raster::extract(flowdir_rast, GTD_Colombia_xy)

# ============================================================
# 3) KEEP ONLY FARC AND ELN
# ============================================================

FARC <- GTD_Colombia %>%
  filter(Group == "FARC")

ELN <- GTD_Colombia %>%
  filter(Group == "ELN")

events <- bind_rows(FARC, ELN) %>%
  mutate(
    # Binary outcome:
    # 1 = ELN
    # 0 = FARC
    ELN_Y = if_else(Group == "Islamic State of Iraq and the Levant (ELNIL)", 1L, 0L),
    
    # Indicators for group-specific fields
    is_FARC = if_else(Group == "FARC", 1L, 0L),
    is_ELN   = if_else(Group == "Islamic State of Iraq and the Levant (ELNIL)", 1L, 0L)
  )

# ============================================================
# 4) PREPARE EVENTS DATA
# ============================================================

events <- events %>%
  mutate(
    row_id_internal = dplyr::row_number(),
    Country_clean   = clean_text(Country),
    Province_clean  = clean_text(Province),
    City_clean      = clean_text(City)
  )

events_sf <- events %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
glimpse(events)
# ============================================================
# 5) PROJECT TO A METRIC CRS
# ============================================================
# IMPORTANT:
# Because the data contain both Colombia and Colombia, a single-country UTM
# projection based on mean longitude/latitude is not appropriate.
# Use a Colombia projected CRS instead.

Colombia_epsg <- 32618

events_proj <- st_transform(events_sf, crs = Colombia_epsg)
coords_m <- st_coordinates(events_proj)

events <- events %>%
  mutate(
    X = coords_m[, 1],
    Y = coords_m[, 2]
  )

# ============================================================
# 6) LOAD NATURAL EARTH POPULATED PLACES
# ============================================================

places <- rnaturalearth::ne_download(
  scale = 10,
  type = "populated_places",
  category = "cultural",
  returnclass = "sf"
)

places_all <- places %>%
  mutate(
    place_name       = NAME,
    feature_class    = FEATURECLA,
    country_name     = ADM0NAME,
    country_clean    = clean_text(ADM0NAME),
    scalerank        = SCALERANK,
    labelrank        = LABELRANK,
    pop_max          = dplyr::coalesce(POP_MAX, 0),
    place_name_clean = clean_text(NAME)
  )

cap_points <- places_all %>%
  filter(
    feature_class %in% c(
      "Admin-0 capital",
      "Admin-0 capital alt",
      "Admin-0 region capital",
      "Admin-1 capital",
      "Admin-1 region capital"
    )
  )

# ============================================================
# 7) LOAD GADM SHAPEFILES
# ============================================================

shape_base_dir <- "C:/R Portfolio/Terrorism_Different_Civil_War/GADM_Border_Distances_Nov_2024/Shape_Files"

country_gadm_lookup <- tibble::tribble(
  ~Country,    ~gadm_code,
  "Colombia",  "COL"
)

gadm_meta <- country_gadm_lookup %>%
  filter(Country %in% unique(events$Country))

gadm_shapes <- purrr::pmap(
  list(gadm_meta$Country, gadm_meta$gadm_code),
  ~ load_gadm_country(..1, ..2, shape_base_dir)
)

names(gadm_shapes) <- gadm_meta$Country

adm0_all <- bind_rows(lapply(gadm_shapes, `[[`, "adm0")) %>%
  st_make_valid()

adm1_all <- bind_rows(lapply(gadm_shapes, `[[`, "adm1")) %>%
  st_make_valid()

# ============================================================
# 8) IDENTIFY NATIONAL CAPITALS
# ============================================================

nat_cap_points <- cap_points %>%
  filter(
    feature_class %in% c(
      "Admin-0 capital",
      "Admin-0 capital alt",
      "Admin-0 region capital"
    )
  )

nat_caps_joined <- st_join(
  nat_cap_points,
  adm0_all,
  join = st_within,
  left = FALSE
)

national_capitals <- nat_caps_joined %>%
  st_drop_geometry() %>%
  mutate(
    nat_cap_lon = LONGITUDE,
    nat_cap_lat = LATITUDE
  ) %>%
  arrange(Country, labelrank, scalerank) %>%
  group_by(Country) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    Country,
    Country_clean = clean_text(Country),
    National_Capital_Name = place_name,
    nat_cap_lon,
    nat_cap_lat
  )

# ============================================================
# 9) IDENTIFY PROVINCIAL / ADMIN1 CAPITALS
# ============================================================

places_admin1_joined <- st_join(
  places_all,
  adm1_all,
  join = st_within,
  left = FALSE
)

if (!("NAME_1" %in% names(places_admin1_joined))) {
  stop("Expected NAME_1 column not found after joining populated places to admin1 polygons.")
}

places_admin1_ranked <- places_admin1_joined %>%
  mutate(
    capital_priority = case_when(
      feature_class %in% c("Admin-1 capital", "Admin-1 region capital") ~ 1L,
      feature_class %in% c("Admin-0 capital", "Admin-0 capital alt", "Admin-0 region capital") ~ 2L,
      TRUE ~ 3L
    ),
    prov_cap_lon = LONGITUDE,
    prov_cap_lat = LATITUDE
  )

provincial_capitals <- places_admin1_ranked %>%
  st_drop_geometry() %>%
  arrange(
    Country,
    NAME_1,
    capital_priority,
    labelrank,
    scalerank,
    desc(pop_max)
  ) %>%
  group_by(Country, NAME_1) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    Country,
    Country_clean = clean_text(Country),
    GADM_NAME_1 = NAME_1,
    GADM_NAME_1_clean = clean_text(NAME_1),
    Province_capital_name = place_name,
    prov_cap_lon,
    prov_cap_lat,
    capital_priority
  )

# ============================================================
# 10) ASSIGN EVENTS TO ADMIN1 POLYGONS
# ============================================================

events_with_admin1 <- st_join(
  events_sf,
  dplyr::select(adm1_all, Country, GADM_NAME_1, GADM_NAME_1_clean),
  join = st_intersects,
  left = TRUE
) %>%
  mutate(
    event_country_clean = clean_text(Country.x),
    join_country_clean  = clean_text(Country.y)
  ) %>%
  st_drop_geometry() %>%
  group_by(row_id_internal) %>%
  slice(1) %>%
  ungroup() %>%
  rename(
    Event_Country = Country.x,
    GADM_Country  = Country.y
  )

# ============================================================
# 11) JOIN CAPITAL LOOKUPS ONTO EVENTS
# ============================================================

events_capitals <- events_with_admin1 %>%
  left_join(
    national_capitals %>%
      dplyr::select(Country, National_Capital_Name, nat_cap_lon, nat_cap_lat),
    by = c("Event_Country" = "Country")
  ) %>%
  left_join(
    provincial_capitals %>%
      dplyr::select(
        Country,
        GADM_NAME_1,
        Province_capital_name,
        prov_cap_lon,
        prov_cap_lat,
        capital_priority
      ),
    by = c("Event_Country" = "Country", "GADM_NAME_1" = "GADM_NAME_1")
  )

# ============================================================
# 12) COMPUTE DELNTANCES TO CAPITALS
# ============================================================

events_capitals <- events_capitals %>%
  mutate(
    event_lon = Longitude,
    event_lat = Latitude
  ) %>%
  rowwise() %>%
  mutate(
    Dist_to_National_Capital_km = if (!is.na(nat_cap_lon) & !is.na(nat_cap_lat)) {
      geosphere::distHaversine(
        c(event_lon, event_lat),
        c(nat_cap_lon, nat_cap_lat)
      ) / 1000
    } else {
      NA_real_
    },
    Dist_to_Provincial_Capital_km = if (!is.na(prov_cap_lon) & !is.na(prov_cap_lat)) {
      geosphere::distHaversine(
        c(event_lon, event_lat),
        c(prov_cap_lon, prov_cap_lat)
      ) / 1000
    } else {
      NA_real_
    }
  ) %>%
  ungroup()

# ============================================================
# 13) CREATE CAPITAL FLAGS
# ============================================================

capital_threshold_km <- 25

events_capitals <- events_capitals %>%
  mutate(
    National_Capital = if_else(
      !is.na(Dist_to_National_Capital_km) &
        Dist_to_National_Capital_km <= capital_threshold_km,
      1L, 0L
    ),
    Provincial_Capital = if_else(
      !is.na(Dist_to_Provincial_Capital_km) &
        Dist_to_Provincial_Capital_km <= capital_threshold_km,
      1L, 0L
    )
  )

# ============================================================
# 14) OPTIONAL NAME-BASED CAPITAL MATCH FLAGS
# ============================================================

events_capitals <- events_capitals %>%
  mutate(
    National_Capital_Name_Match = if_else(
      !is.na(National_Capital_Name) &
        City_clean == clean_text(National_Capital_Name),
      1L, 0L
    ),
    Provincial_Capital_Name_Match = if_else(
      !is.na(Province_capital_name) &
        City_clean == clean_text(Province_capital_name),
      1L, 0L
    )
  )

# ============================================================
# 15) JOIN CAPITAL VARIABLES BACK
# ============================================================

capital_vars <- events_capitals %>%
  dplyr::select(
    row_id_internal,
    GADM_NAME_1,
    National_Capital_Name,
    Province_capital_name,
    Dist_to_National_Capital_km,
    Dist_to_Provincial_Capital_km,
    National_Capital,
    Provincial_Capital,
    National_Capital_Name_Match,
    Provincial_Capital_Name_Match
  )

events <- events %>%
  left_join(capital_vars, by = "row_id_internal")

#######################
# 4) BORDER DISTANCES #
#######################

mean_lon <- mean(events$Longitude, na.rm = TRUE)
mean_lat <- mean(events$Latitude, na.rm = TRUE)

utm_zone <- floor((mean_lon + 180) / 6) + 1
epsg_utm <- if (mean_lat >= 0) 32600 + utm_zone else 32700 + utm_zone
epsg_utm

# Distance to outer national border
events <- add_border_distance(
  data = events,
  border_path = "C:/R Portfolio/Terrorism_Different_Civil_War/GADM_Border_Distances_Nov_2024/Shape_Files/gadm41_COL_shp/gadm41_COL_0.shp",
  dist_prefix = "B",
  internal_only = FALSE,
  utm_epsg = 32618,
  use_s2 = FALSE,
  save_rds_path = "C:/R Portfolio/Terrorism_Different_Civil_War/Latent_Field_Analyses/Colombia_B_distances.rds"
)

# Distance to internal state / provincial borders
events <- add_border_distance(
  data = events,
  border_path = "C:/R Portfolio/Terrorism_Different_Civil_War/GADM_Border_Distances_Nov_2024/Shape_Files/gadm41_COL_shp/gadm41_COL_1.shp",
  dist_prefix = "PB",
  internal_only = TRUE,
  utm_epsg = 32618,
  use_s2 = FALSE,
  save_rds_path = "C:/R Portfolio/Terrorism_Different_Civil_War/Latent_Field_Analyses/Colombia_PB_distances.rds"
)

# ============================================================
# 16) PREPARE MODELLING DATA
# ============================================================

events <- events %>%
  mutate(
    Attack = as.factor(Attack),
    Target = as.factor(Target)
  )

# Include all variables used later in the model
model_vars <- c(
  "ELN_Y",
  "is_FARC",
  "is_ELN",
  "Attack",
  "Target",
  "B_Dist_km",
  "PB_Dist_km",
  "Dist_to_National_Capital_km",
  "Dist_to_Provincial_Capital_km",
  "Elevation",
  "Slope",
  "X",
  "Y"
)

# keep_idx <- complete.cases(events[, model_vars])

events_model <- events %>%
  dplyr::mutate(
    Intercept = 1
  )

events_model <- events_model %>%
  mutate(
    ELN_Y = if_else(Group == "ELN", 1L, 0L),
    is_ELN = if_else(Group == "ELN", 1L, 0L),
    is_FARC = if_else(Group == "FARC", 1L, 0L)
  )

table(events_model$Group, events_model$ELN_Y)
table(events_model$is_ELN)
table(events_model$is_FARC)

# events <- events_model %>%
#   mutate(
#     Dist_to_National_Capital_km   = Dist_to_National_Capital_km.x,
#     Dist_to_Provincial_Capital_km = Dist_to_Provincial_Capital_km.x,
#     National_Capital              = National_Capital.x,
#     Provincial_Capital            = Provincial_Capital.x,
#     National_Capital_Name_Match   = National_Capital_Name_Match.x,
#     Provincial_Capital_Name_Match = Provincial_Capital_Name_Match.x,
#     GADM_NAME_1                   = GADM_NAME_1.x,
#     National_Capital_Name         = National_Capital_Name.x,
#     Province_capital_name         = Province_capital_name.x
#   ) %>%
#   dplyr::select(
#     -ends_with(".x"),
#     -ends_with(".y")
#   )

coords_m_model <- as.matrix(events_model[, c("X", "Y")])

glimpse(events_model)

table(events_model$Group, events_model$ELN_Y)

# ============================================================
# 18) DEFINE ONE SPDE MODEL
#    Use one sensible prior first, not a full grid on M2
# ============================================================

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(50000, 0.5),   # P(range < 50 km) = 0.5
  prior.sigma = c(1, 0.01)       # P(sigma > 1) = 0.01
)

# ============================================================
# 19) CREATE BASE A MATRICES
# ============================================================

A_base <- inla.spde.make.A(
  mesh = mesh,
  loc = coords_m_model
)

# Shared field applies to all rows
A_shared <- A_base

# ELN-specific field only applies where is_ELN == 1
A_eln <- Matrix::Diagonal(x = events_model$is_ELN) %*% A_base

# FARC-specific field only applies where is_FARC == 1
A_farc <- Matrix::Diagonal(x = events_model$is_FARC) %*% A_base

# Quick checks
mesh$n
dim(A_base)
table(events_model$is_ELN)
table(events_model$is_FARC)
table(events_model$ELN_Y)
table(events_model$Group, events_model$ELN_Y)

# ============================================================
# 20) CREATE INDICES
# ============================================================

shared_index <- inla.spde.make.index(
  name = "shared_field",
  n.spde = spde$n.spde
)

eln_index <- inla.spde.make.index(
  name = "eln_field",
  n.spde = spde$n.spde
)

farc_index <- inla.spde.make.index(
  name = "farc_field",
  n.spde = spde$n.spde
)

# ============================================================
# 21) M0 = FIXED EFFECTS ONLY
# ============================================================

stack_m0 <- inla.stack(
  data = list(y = events_model$ELN_Y),
  A = list(1),
  effects = list(
    data.frame(
      Intercept = events_model$Intercept,
      Attack = events_model$Attack,
      Target = events_model$Target,
      B_Dist_km = events_model$B_Dist_km,
      PB_Dist_km = events_model$PB_Dist_km,
      Dist_to_National_Capital_km = events_model$Dist_to_National_Capital_km,
      Dist_to_Provincial_Capital_km = events_model$Dist_to_Provincial_Capital_km,
      Elevation = events_model$Elevation,
      Slope = events_model$Slope
    )
  ),
  tag = "est"
)

formula_m0 <- y ~ 0 + Intercept +
  Attack +
  Target +
  B_Dist_km +
  PB_Dist_km +
  Dist_to_National_Capital_km +
  Dist_to_Provincial_Capital_km +
  Elevation +
  Slope

fit_m0 <- inla(
  formula_m0,
  data = inla.stack.data(stack_m0),
  family = "binomial",
  control.predictor = list(
    A = inla.stack.A(stack_m0),
    compute = TRUE,
    link = 1
  ),
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE
  ),
  control.inla = list(
    strategy = "adaptive"
  )
)

# ============================================================
# 22) M1 = SHARED FIELD ONLY
# ============================================================

stack_m1 <- inla.stack(
  data = list(y = events_model$ELN_Y),
  A = list(1, A_shared),
  effects = list(
    data.frame(
      Intercept = events_model$Intercept,
      Attack = events_model$Attack,
      Target = events_model$Target,
      B_Dist_km = events_model$B_Dist_km,
      PB_Dist_km = events_model$PB_Dist_km,
      Dist_to_National_Capital_km = events_model$Dist_to_National_Capital_km,
      Dist_to_Provincial_Capital_km = events_model$Dist_to_Provincial_Capital_km,
      Elevation = events_model$Elevation,
      Slope = events_model$Slope
    ),
    shared_field = shared_index
  ),
  tag = "est"
)

formula_m1 <- y ~ 0 + Intercept +
  Attack +
  Target +
  B_Dist_km +
  PB_Dist_km +
  Dist_to_National_Capital_km +
  Dist_to_Provincial_Capital_km +
  Elevation +
  Slope +
  f(shared_field, model = spde)

fit_m1 <- inla(
  formula_m1,
  data = inla.stack.data(stack_m1),
  family = "binomial",
  control.predictor = list(
    A = inla.stack.A(stack_m1),
    compute = TRUE,
    link = 1
  ),
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE
  ),
  control.inla = list(
    strategy = "adaptive"
  )
)

# ============================================================
# 23) M2 = SHARED + ELN + FARC FIELDS
# ============================================================

stack_m2 <- inla.stack(
  data = list(y = events_model$ELN_Y),
  A = list(1, A_shared, A_eln, A_farc),
  effects = list(
    data.frame(
      Intercept = events_model$Intercept,
      Attack = events_model$Attack,
      Target = events_model$Target,
      B_Dist_km = events_model$B_Dist_km,
      PB_Dist_km = events_model$PB_Dist_km,
      Dist_to_National_Capital_km = events_model$Dist_to_National_Capital_km,
      Dist_to_Provincial_Capital_km = events_model$Dist_to_Provincial_Capital_km,
      Elevation = events_model$Elevation,
      Slope = events_model$Slope
    ),
    shared_field = shared_index,
    eln_field = eln_index,
    farc_field = farc_index
  ),
  tag = "est"
)

formula_m2 <- y ~ 0 + Intercept +
  Attack +
  Target +
  B_Dist_km +
  PB_Dist_km +
  Dist_to_National_Capital_km +
  Dist_to_Provincial_Capital_km +
  Elevation +
  Slope +
  f(shared_field, model = spde) +
  f(eln_field, model = spde) +
  f(farc_field, model = spde)

fit_m2 <- inla(
  formula_m2,
  data = inla.stack.data(stack_m2),
  family = "binomial",
  control.predictor = list(
    A = inla.stack.A(stack_m2),
    compute = TRUE,
    link = 1
  ),
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE
  ),
  control.inla = list(
    strategy = "adaptive"
  )
)

# ============================================================
# 24) COMPARE MODELS
# ============================================================

model_metrics <- data.frame(
  Model = c("M0_fixed_only", "M1_shared_only", "M2_shared_plus_ELN_plus_FARC"),
  DIC   = c(fit_m0$dic$dic, fit_m1$dic$dic, fit_m2$dic$dic),
  WAIC  = c(fit_m0$waic$waic, fit_m1$waic$waic, fit_m2$waic$waic)
) %>%
  dplyr::arrange(WAIC)

print(model_metrics)

# ============================================================
# 25) SUMMARIES
# ============================================================

summary(fit_m0)
summary(fit_m1)
summary(fit_m2)

fit_m0$summary.fixed
fit_m1$summary.fixed
fit_m2$summary.fixed

fit_m1$summary.hyperpar
fit_m2$summary.hyperpar

# ============================================================
# 26) EXTRACT LATENT FIELD SUMMARIES FOR M2
# ============================================================

shared_summary <- fit_m2$summary.random$shared_field
eln_summary    <- fit_m2$summary.random$eln_field
farc_summary   <- fit_m2$summary.random$farc_field

# ============================================================
# 27) PUT FIELDS INTO DATA FRAMES
# ============================================================

shared_df <- data.frame(
  x = mesh$loc[, 1],
  y = mesh$loc[, 2],
  mean  = shared_summary$mean,
  sd    = shared_summary$sd,
  lower = shared_summary$`0.025quant`,
  upper = shared_summary$`0.975quant`,
  field = "Shared field"
)

eln_df <- data.frame(
  x = mesh$loc[, 1],
  y = mesh$loc[, 2],
  mean  = eln_summary$mean,
  sd    = eln_summary$sd,
  lower = eln_summary$`0.025quant`,
  upper = eln_summary$`0.975quant`,
  field = "ELN-specific field"
)

farc_df <- data.frame(
  x = mesh$loc[, 1],
  y = mesh$loc[, 2],
  mean  = farc_summary$mean,
  sd    = farc_summary$sd,
  lower = farc_summary$`0.025quant`,
  upper = farc_summary$`0.975quant`,
  field = "FARC-specific field"
)

all_fields_df <- bind_rows(shared_df, eln_df, farc_df)

# ============================================================
# 28) PLOT FIELDS AT MESH NODES
# ============================================================

ggplot(all_fields_df, aes(x = x, y = y, color = mean)) +
  geom_point(size = 1.8) +
  coord_equal() +
  facet_wrap(~ field) +
  theme_classic() +
  labs(
    title = "Shared and Group-Specific Latent Fields",
    x = "Projected X",
    y = "Projected Y",
    color = "Posterior mean"
  )

# ============================================================
# 29) PROJECT TO REGULAR GRID FOR SMOOTHER MAPS
# ============================================================

proj <- inla.mesh.projector(
  mesh,
  xlim = range(mesh$loc[, 1]),
  ylim = range(mesh$loc[, 2]),
  dims = c(200, 200)
)

shared_grid <- inla.mesh.project(
  projector = proj,
  field = shared_summary$mean
)

eln_grid <- inla.mesh.project(
  projector = proj,
  field = eln_summary$mean
)

farc_grid <- inla.mesh.project(
  projector = proj,
  field = farc_summary$mean
)

grid_template <- expand.grid(
  x = proj$x,
  y = proj$y
)

shared_grid_df <- grid_template %>%
  mutate(
    mean = as.vector(shared_grid),
    field = "Shared field"
  )

eln_grid_df <- grid_template %>%
  mutate(
    mean = as.vector(eln_grid),
    field = "ELN-specific field"
  )

farc_grid_df <- grid_template %>%
  mutate(
    mean = as.vector(farc_grid),
    field = "FARC-specific field"
  )

grid_fields_df <- bind_rows(shared_grid_df, eln_grid_df, farc_grid_df)

ggplot(grid_fields_df, aes(x = x, y = y, fill = mean)) +
  geom_raster() +
  coord_equal() +
  facet_wrap(~ field) +
  theme_classic() +
  labs(
    title = "Projected Shared and Group-Specific Latent Fields",
    x = "Projected X",
    y = "Projected Y",
    fill = "Posterior mean"
  )

# ============================================================
# 30) EXTRACT EVENT-LEVEL FITTED VALUES FOR M2
# ============================================================

idx_est <- inla.stack.index(stack_m2, tag = "est")$data

events_model$Pred_Prob_ELN <- fit_m2$summary.fitted.values[idx_est, "mean"]
events_model$Pred_Group    <- ifelse(events_model$Pred_Prob_ELN >= 0.5, "ELN", "FARC")

table(
  Observed  = ifelse(events_model$ELN_Y == 1, "ELN", "FARC"),
  Predicted = events_model$Pred_Group
)
