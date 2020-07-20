

# Outline
# 1. get ranking of top n french cities by population in 2015. This defines the study population. The cities are selected such that they are always present (i.e. not Villeurbane near Lyon and not Boulogne-Billancourt near paris) and not overseas territory
# 2. final output is a data.table with five columns: CODGEO, LIBGEO, year, area, population.
# 3. year = c(1876, 1950, seq(1975,2015, by = 5))
# 4. area is a visual measure in 1876 and 1950, automatized from raster data after
# 5. manual measures are in online google sheet

#' LandUse pipeline
#'
#' @export
pipe <- function(){
    p = readpop()  # census pop counts time series
    g = read_areas_gsheet(overwrite = FALSE)

    f = cropFrance()
    f[,city_sane := gsub(pattern = " \\[FRA\\].*$",replacement = "",x=city)]
    f[city == "Nice-Cannes [FRA]", city_sane := "Nice"]
    f[city == "Lille [FRA]; Mouscron [BEL]", city_sane := "Lille"]
    f[city == "Strasbourg [FRA]; Kehl [DEU]", city_sane := "Strasbourg"]
    f[,density := pop / built]
    m1 = merge(f,g,by.x = c("city_sane","year"),by.y = c("city","year"),all.x = TRUE)
    m50 = merge(f[year == 2015,list(city,city_sane,year = 1950,pop = NA,built = NA)],g[,list(city,year,built_manual)],by.x = c("city_sane","year"),by.y = c("city","year"),all.x = TRUE)
    m50[, density := NA]
    m1860 = merge(f[year == 2015,list(city,city_sane,year = 1866,pop = NA,built = NA)],g[,list(city,year,built_manual)],by.x = c("city_sane","year"),by.y = c("city","year"),all.x = TRUE)
    m1860[,density := NA]

    paris <- p %>% dplyr::filter(grepl("Paris ",LIBGEO), year %in% c(1876,1954,1975,1990,2000,2015)) %>%
        dplyr::select(-c(REG,DEP,date)) %>%
        as.data.table
    paris[,city_sane := gsub(" ","-",gsub(" Arrondissement","",LIBGEO))]
    paris[,c("built","density","built_manual","source") := NA]
    paris[year == 1876, year := 1866]
    setnames(paris, c("LIBGEO","population"), c("city","pop"))

    # get areas of paris arrondissements
    a <- readRDS(file.path(LandUseR:::datadir(),"paris-areas.Rds")) %>% as.data.table
    setnames(a,"area","built")
    a[,built := as.numeric(built)]
    setkey(paris, CODGEO, year)
    setkey(a, CODGEO)
    paris <- a[paris[,list(CODGEO,city_sane,year,city,density, built_manual, source, pop)]]

    m = rbind(m1,m50,m1860)
    m[,source := "GHS-SMOD"]
    m[year < 1975 , source := "manual"]
    m[,city_ss := city_sane]
    m[city_sane == "Orleans", city_ss := "Orléans"]
    m[city_sane == "Saint-Etienne", city_ss := "Saint-Étienne"]

    # get CODGEO
    cg = p %>%
        dplyr::filter(LIBGEO %in% m[,unique(city_ss)], year==2015) %>%
        dplyr::select(CODGEO, LIBGEO) %>%
        dplyr::bind_rows(data.frame(CODGEO = "75000",LIBGEO = "Paris",stringsAsFactors = FALSE)) %>%
        dplyr::mutate(city_ss = LIBGEO) %>%
        dplyr::select(CODGEO,city_ss) %>%
        as.data.table
    setkey(cg,city_ss)
    setkey(m,city_ss)
    m = cg[m]
    m[,city_ss := NULL]

    # fill in manual built measures
    m[source == "manual", built := built_manual]

    out = LandUseR:::fill_pops(m,p)  # adds info before 1975

    out = rbind(m[year > 1950],out)
    out = rbind(out,paris)
    out[,density := pop / built]
    saveRDS(out, file.path(datadir(),"areas_pops.Rds"))
    fwrite(out, file.path(datadir(),"areas_pops.csv"))

    # output plot
    pl = list()
    pl$all = plot_france(out,fname = "densities-1866-FRA")
    pl$only50 = plot_france(out[year>1949],fname = "densities-1950-FRA")
    pl$only75 = plot_france(out[year>1970],fname = "densities-1975-FRA",dotted = TRUE)
    pl$pop_noparis   = plot_france(out[city_sane != "Paris"],outcome = "pop", fname = "population-1866-FRA-noparis")
    pl$pop_onlyparis = plot_france(out[city_sane == "Paris"],outcome = "pop", fname = "population-1866-FRA-paris")
    pl$built_noparis = plot_france(out[city_sane != "Paris"],outcome = "built", fname = "built-1866-FRA-noparis")
    pl$built_onlyparis = plot_france(out[city_sane == "Paris"],outcome = "built", fname = "built-1866-FRA-paris")

    # plot densities 1866 vs 1950
    # wide = dcast(out, CODGEO + LIBGEO ~ year, value.var = "density" )
    # tolabel = c("42218","76351","59350","69123","06088","75000")
    # pl = ggplot(wide, aes(x = `1866`,y = `1950`,label = LIBGEO)) + theme_bw() + geom_abline(slope = 1, intercept = 0) + coord_fixed(xlim = c(5000,120000),ylim = c(5000,120000)) + geom_point() + labs(title = "Density in 1866 vs 1950", subtitle = "people per km2")
    # pl = pl + geom_text(data = wide[wide$CODGEO %in% tolabel],vjust = 0, nudge_y = 3000)
    # ggsave(plot = pl, filename = file.path(datadir(),"..","output","plots","densities.pdf"),width = 7, height = 6)
    #
    return(list(out,pl))
}


setup_gsheet <- function(){
    googledrive::drive_upload(file.path(datadir(),"france_cities.csv"), type = "spreadsheet")
}

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

plot_paris_sat <- function(){
    cc = LandUseR:::cropcities()

    png(file.path(LandUseR:::outdir(),"data","plots","paris-sat-area.png"), width = 1000, height=700)
    par(mar = c(1.2, 0.1, 1.5, 0.1), mfrow=c(2,2))
    for (yr in paste(LandUseR:::GHS_years())){
        raster::plot(cc[[yr]][[1]]$built,main = paste("Paris % built-up",yr),xaxt = "n",yaxt = "n", legend = FALSE)
    }
    dev.off()
    png(file.path(LandUseR:::outdir(),"data","plots","Lyon-sat-area.png"), width = 1000, height=700)
    par(mar = c(1.2, 0.1, 1.5, 0.1), mfrow=c(2,2))
    # par( mfrow=c(2,2))
    for (yr in paste(LandUseR:::GHS_years())){
        raster::plot(cc[[yr]][[2]]$built,main = paste("Lyon % built-up",yr),xaxt = "n",yaxt = "n",legend = FALSE)
    }
    dev.off()
    png(file.path(LandUseR:::outdir(),"data","plots","Marseille-sat-area.png"), width = 1000, height=700)
    par(mar = c(1.2, 0.1, 1.5, 0.1), mfrow=c(2,2))
    for (yr in paste(LandUseR:::GHS_years())){
        raster::plot(cc[[yr]][[3]]$built,main = paste("Marseille % built-up",yr),xaxt = "n",yaxt = "n",legend = FALSE)
    }
    dev.off()
    png(file.path(LandUseR:::outdir(),"data","plots","Reims-sat-area.png"), width = 1000, height=700)
    par(mar = c(1.2, 0.1, 1.5, 0.1), mfrow=c(2,2))
    for (yr in paste(LandUseR:::GHS_years())){
        raster::plot(cc[[yr]][["51454"]]$built,main = paste("Reims % built-up",yr),xaxt = "n",yaxt = "n",legend = FALSE)
    }
    dev.off()

    png(file.path(LandUseR:::outdir(),"data","plots","Toulouse-sat-area.png"), width = 1000, height=700)
    par(mar = c(1.2, 0.1, 1.5, 0.1), mfrow=c(2,2))
    # par(mfrow=c(2,2))
    for (yr in paste(LandUseR:::GHS_years())){
        raster::plot(cc[[yr]][["31555"]]$built,main = paste("Toulouse % built-up",yr),xaxt = "n",yaxt = "n", legend = FALSE)
    }
    dev.off()
}





#' get a ranking of cities by size
#'
#' this is somewhat arbitrary as it includes communes which are clearly part of larger agglomerations. the population column is not the final product, in other words.
#' excludes boulogne and versaille, villeurbane and La Reunion
#' @export
city_ranking <- function(p,yr = 2015, topn=30, overwrite = FALSE){
    if (!file.exists(file.path(datadir(),paste0("top",topn,"_",yr,".Rds"))) | overwrite ){
        paris = p %>% dplyr::filter(year == yr, grepl("Paris ",LIBGEO)) %>% dplyr::summarise(CODGEO = "75000", REG = 11, DEP = 75, LIBGEO = "Paris", year = yr, population = sum(population), date = date[1])
        t30 = p %>% dplyr::filter(year == yr, !grepl("Paris ",LIBGEO), !(CODGEO %in%  c("92012","78646","69266", "97411"))) %>% dplyr::top_n(n = topn, wt = population) %>% dplyr::arrange(desc(population))

        t30 %<>%
            dplyr::bind_rows(paris) %>%
            dplyr::arrange(desc(population))  %>%
            dplyr::mutate(rank = row_number()) %>%
            dplyr::select(CODGEO,LIBGEO) %>%
            # dplyr::rename(rank_pop = population) %>%
            as.data.table()
        saveRDS(t30,file.path(datadir(),paste0("top",topn,"_",yr,".Rds")))

    } else {
        t30 = readRDS(file.path(datadir(),paste0("top",topn,"_",yr,".Rds")))
    }
    return(t30)
}
