## trkC_gain PREREG 2.6: ESI v6r NA 1 km at the public (perturbed) DataMart PLOT coordinates of GA and IN; output CN + esi only (restricted, 600).
suppressMessages(library(terra)); library(data.table)
r <- rast("/home/aaron/jobs/na30/release_v1/A_esi/data/ESI_v6r_NA_1km.tif")
for (st in c("GA", "IN")) {
  P <- fread(sprintf("/home/aaron/restricted/ncasi-modeleval/fia_data_rc/%s/%s_PLOT.csv", st, st), select = c("CN", "LAT", "LON"), colClasses = list(character = "CN"))
  P <- P[is.finite(LAT) & is.finite(LON)]; p <- project(vect(P, geom = c("LON", "LAT"), crs = "EPSG:4269"), crs(r)); P[, esi := extract(r[[1]], p)[, 2]]
  f <- sprintf("/home/aaron/perseus_me/tracks/trkC_gain/restricted/esi_plot_%s.csv", st); fwrite(P[, .(CN, esi)], f); Sys.chmod(f, "600")
  cat(st, "plots", nrow(P), "esi missing share", round(mean(is.na(P$esi)), 4), "q01 q50 q99", round(quantile(P$esi, c(.01, .5, .99), na.rm = TRUE), 2), "\n")
}
cat("raster md5", tools::md5sum("/home/aaron/jobs/na30/release_v1/A_esi/data/ESI_v6r_NA_1km.tif"), "\n")
