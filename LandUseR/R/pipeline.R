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
#' By default selects the version with `cutoff = 30`
#' @export
read_output <- function(cutoff = 30){
    readRDS(file.path(outdatadir(),paste0("france_final_cutoff",30,".Rds")))
}

#' Show Relative population and area
#'
#' output table for top5 cities with relative population and area to Paris.
#' @export
relative_pop_area <- function(cities = top5now(), overwrite = FALSE){
    if (overwrite){
        x = read_output()
        y = x[(LIBGEO %in% cities) & (year == 2015), list(LIBGEO,relative_pop = pop / max(pop), relative_area = area / max(area))][order(relative_pop)]
        fwrite(y, file.path(outdatadir(),"top5poparea.csv"))
        y
    } else {
        fread(file.path(outdatadir(),"top5poparea.csv"))
    }
}

#' Main Population Count Output File
#'
#' reads output of manual and satellite measures via \code{\link{read_output}}
#' and amends Paris Population to also include the Seine department as well before WW2.
#'
#' The csv output of this file is used in the model.
pop_allyears <- function(){
    x = read_output()
    p0 = data.table(readpop())
    p1 = p0[DEP == 75, list(CODGEO = "75060", LIBGEO = "Paris", population = sum(population)), by = year(date)]
    p = rbind(p0[,list(CODGEO,LIBGEO,population,year)],p1[,list(CODGEO,LIBGEO,population,year)])
    x[year == 1950 , year := 1954] # fix census year
    p = p[CODGEO %in% x[,unique(CODGEO)]]
    y = merge(x[,list(CODGEO,rank,pop, area, type,year)], p, all.y = TRUE, by = c("CODGEO","year"))

    # paris census data is incorrect - paris is larger than "paris" even before WW2
    pa = fread(file.path(LandUseR:::datadir(),"paris-seine.csv"))
    setnames(pa, "population", "pop2")
    pa[,CODGEO := "75060"]
    y = merge(y, pa, all.x = TRUE)
    y[ year < 1954 & CODGEO == "75060", population := pop2]
    y[,pop2 := NULL]

    # before 1954, take census measures, after take satellite
    y[ year < 1954, pop := population]
    y[ , rank := .SD[year == 1876, rank], by = CODGEO]
    y <- y[!is.na(pop)]
    y = y[, list(CODGEO = as.character(CODGEO),LIBGEO, rank,density_data = pop / area, reldensity_data = (pop / area) / max(pop / area), area_data = area,relarea_data = area / max(area),pop_data = pop,relpop_data = pop / max(pop)),by=year]
    fwrite(y, file.path(outdatadir(),"relpop.csv"))

    p0 = ggplot(y[rank < 21], aes(x=year, y = pop, color = LIBGEO)) + geom_line()
    p01 = ggplot(y[rank < 21], aes(x=pop_data, color = year)) + geom_density()
    p1 = ggplot(y[rank < 21][rank != 1], aes(x=year, y = relpop_data, color = LIBGEO)) + geom_line() + geom_point() + ggtitle("Population relative to Paris")
    ggsave(p1, file = file.path(dataplots(),"relpop-paris.pdf"))
    p2 = ggplot(y[rank != 1], aes(x=year, y = relpop_data, color = LIBGEO)) + geom_line() + theme(legend.position = "none")
    list(y, p0,p01,p1,p2)

}
