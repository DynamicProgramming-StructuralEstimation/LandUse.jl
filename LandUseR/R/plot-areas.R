






#' Read main output table
#'
#' @export
read_output <- function(){
    readRDS(file.path(datadir(),"areas_pops.Rds"))
}


#' read area data
#'
#' @export
read_areas_gsheet <- function(overwrite = FALSE){
    if (!file.exists(file.path(datadir(),"areas.Rds")) | overwrite){
        url = "https://docs.google.com/spreadsheets/d/e/2PACX-1vQ-HXqJgJ-RDYZGDvlDh-2kfCy8yaTk-BtmOMqzgVGlwGCSev7U_mV20B8WZoePjGAK01hLAtVjhinJ/pub?gid=302223516&single=true&output=csv"
        gs = readr::read_csv(url) %>%
            mutate(city = Commune) %>%
            dplyr::select(city, area_2016, area_1950, area_etat_major)
        readr::write_csv(gs,file.path(datadir(),"areas_gsheet.csv"))
        gd = data.table::as.data.table(gs)
        data.table::setnames(gd, c("area_2016","area_1950","area_etat_major"), c("2015","1950","1866"))
        gd = data.table::melt(gd, id.vars = "city", variable.name = "year", value.name = "built_manual" )
        # gd[, population := NA_integer_]
        gd[, year := as.integer(as.character(year))]
        gd <- gd[complete.cases(gd)]
        saveRDS(gd,file.path(datadir(),"areas.Rds"))

    } else {
        gd = readRDS(file.path(datadir(),"areas.Rds"))
    }
    gd
}


#' takes city ranking data.table and fills in populations for 1950 and 1866
#'
#' we take census 1876 for 1866 area
fill_pops <- function(ci,p) {
    pd = as.data.table(p)

    setkey(ci,year,CODGEO)
    setkey(pd,year,CODGEO)

    # 1866 / 1876
    # easy: just merge in official census data and the sum of all paris arrondissements
    li = list()
    li$EM <- merge(ci[J(1866),!"pop"],pd[J(1876),list(CODGEO,pop=population)],by = "CODGEO",all.x = TRUE)
    li$EM[J("75000"), c("pop") := list(pd[year == 1876 & DEP == 75,sum(population)]) ]

    # 1950 / 1954
    # less easy: need to manually compute pop for our definitions of "the city of x"
    p50 <- pop_1950(p)
    setkey(p50,year, CODGEO)
    li$fifties <- merge(ci[J(1950),!"pop"],p50[J(1954),list(CODGEO,pop=population)],by = "CODGEO",all.x = TRUE)

    return(rbindlist(li))
}



plot_paris <- function(d){
    pl = d %>%
        dplyr::filter(DEP == 75) %>%
        dplyr::mutate(Arrond = as.factor(as.numeric(CODGEO)-75100)) %>%
        ggplot2::ggplot(aes(x = year, y = population, group = CODGEO, color = Arrond)) + ggplot2::geom_line(size=1.1) + ggplot2::scale_y_continuous(labels = scales::comma) + ggplot2::labs(title = "Central Paris Population by Arrondissement",subtitle = "Equivalent to `Density` (Area is constant)") + ggplot2::theme_bw()
    ggplot2::ggsave(file.path(outdir(),"plots","paris.pdf"),plot = pl, width = 8, height = 6)
    pl

}



