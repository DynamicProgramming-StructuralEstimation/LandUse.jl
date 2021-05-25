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

#' Read main output table
#'
#' @export
read_output <- function(){
    readRDS(file.path(outdatadir(),"france_final.Rds"))
}

#' Show Relative population and area
#'
#' output table for top5 cities with relative population and area to Paris.
#' @export
relative_pop_area <- function(cities = top5now()){
    x = read_output()
    x[(LIBGEO %in% cities) & (year == 2015), list(LIBGEO,relative_pop = pop / max(pop), relative_area = area / max(area))]
}
