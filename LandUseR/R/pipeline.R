# Outline
# 1. get ranking of top n french cities by population in 2015. This defines the study population. The cities are selected such that they are always present (i.e. not Villeurbane near Lyon and not Boulogne-Billancourt near paris) and not overseas territory
# 2. final output is a data.table with five columns: CODGEO, LIBGEO, year, area, population.
# 3. year = c(1876, 1950, seq(1975,2015, by = 5))
# 4. area is a visual measure in 1876 and 1950, automatized from raster data after
# 5. manual measures are in online google sheet

#' LandUse pipeline
#'
#' @export
pipe <- function(overwrite = FALSE){
    pop_1950_2(overwrite)
    cc = cropcities(overwrite)
    cm = measure_cities(overwrite)
    mm = combine_measures(overwrite)
}
