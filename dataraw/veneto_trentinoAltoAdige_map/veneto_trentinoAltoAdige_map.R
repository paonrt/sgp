### shapefile downloaded at the following link:
# https://www.istat.it/notizia/confini-delle-unita-amministrative-a-fini-statistici-al-1-gennaio-2018-2/ # nolint: line_length_linter.

### website in italian, Versione generalizzata (meno dettagliata) 2024, direct
###  link:
# https://www.istat.it/storage/cartografia/confini_amministrativi/generalizzati/2024/Limiti01012024_g.zip # nolint: line_length_linter.
### downloaded 04/09/2024

### coordinate system: WGS84 UTM32N

### go to folder Reg01012024_g, copy and paste the files in the main directory

#' @import sf
#' @import ggplot2
#' @import oce
#' @import elevatr


### get polygons from shape file
veneto_trentinoAltoAdige_polygons = sf::st_read("Reg01012024_g_WGS84.shp")

### restrict to regions of interest
veneto_trentinoAltoAdige_polygons = veneto_trentinoAltoAdige_polygons[
  veneto_trentinoAltoAdige_polygons$DEN_REG %in% c(
    "Trentino-Alto Adige", "Veneto"
  ),
]

#veneto_trentinoAltoAdige_polygons$geometry
xmin = 606006
ymin = 4965550
xmax = 819684.5
ymax = 5220292

### resolutions of the map
resolutionX = 90
resolutionY = 120

deltaX = (xmax - xmin) / resolutionX
deltaY = (ymax - ymin) / resolutionY

### get a grid of the map
points = array(
  0, dim = c(resolutionY, resolutionX, 3),
  dimnames = list(NULL, NULL, c("longitude", "latitude", "is_inPolygons"))
)

### remove points of the grid outside the map
for (lat in 1:resolutionY) {
  #print(paste(lat, " at time ", Sys.time(), sep = ""))
  for (lon in 1:resolutionX) {
    points[lat, lon, 1:2] = c(xmin + deltaX * lon, ymin + deltaY * lat)
    point = sf::st_sfc(
      sf::st_point(points[lat, lon, 1:2]),
      crs = sf::st_crs(veneto_trentinoAltoAdige_polygons)
    )
    contains = sum(unlist(sf::st_contains(
      veneto_trentinoAltoAdige_polygons, point
    )))
    if (identical(contains, 1L)) {
      points[lat, lon, 3] = 1
    } else {
      points[lat, lon, 3] = 0
    }
  }
}
true_points = as.logical(points[, , 3])
truePoints_lon = c(points[, , 1])[true_points]
truePoints_lat = c(points[, , 2])[true_points]




########## get altitude

### get utm coordinates
coordinate_utm = oce::utm2lonlat(truePoints_lon, truePoints_lat, 32)
coordinate_utm = data.frame(
  x = coordinate_utm$longitude, y = coordinate_utm$latitude
)

### get elevation
truePoints_alt = elevatr::get_elev_point(
  locations = coordinate_utm, prj = 4326, src = "aws"
)
truePoints_alt = truePoints_alt$elevation



veneto_trentinoAltoAdige_map = list()
veneto_trentinoAltoAdige_map$polygons = veneto_trentinoAltoAdige_polygons
veneto_trentinoAltoAdige_map$xyz = data.frame(
  longitude = truePoints_lon,
  latitude = truePoints_lat,
  altitude = truePoints_alt
)
usethis::use_data(veneto_trentinoAltoAdige_map)