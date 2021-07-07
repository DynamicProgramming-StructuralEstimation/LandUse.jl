
# Plotting Functions


#' Plot distance from center
#'
#' using GHS grid data, plot average density at a certain distance
#' from the center of a city.
#'
#' This function uses output of \code{\link{dist_from_center}} and
#' estimates an exponential decay model via NLS.
#'
#' @param city_name vector of city names
#'
#' @return list with plots
#'
#' @references https://douglas-watson.github.io/post/2018-09_exponential_curve_fitting/
plot_density_center <- function(city_name = c("Paris","Lyon")){
    d0 = dist_from_center()
    d0[, distance := distance / 1000]
    l = list()
    lnls = list()
    for (ic in city_name){
        d = d0[LIBGEO == ic & !is.na(distance)]

        # estimate nls model of exponential decay
        # https://douglas-watson.github.io/post/2018-09_exponential_curve_fitting/
        fit <- nls(density ~ SSasymp(distance, yf, y0, log_alpha), data = d)
        fitted <- d %>%
            tidyr::nest(-year) %>%
            mutate(
                fit = purrr::map(data, ~nls(density ~ SSasymp(distance, yf, y0, log_alpha), data = .)),
                tidied = purrr::map(fit, tidy),
                augmented = purrr::map(fit, augment)
            )

        tab = fitted %>%
                tidyr::unnest(tidied) %>%
                dplyr::select(year, term, estimate) %>%
                tidyr::spread(term, estimate) %>%
                mutate(lambda = exp(log_alpha))

        augmented <- fitted %>%
            tidyr::unnest(augmented)

        l[[ic]] = ggplot(data = augmented, aes(x=distance,y = density, color = factor(year))) +
            geom_point() +
            geom_line(aes(y = .fitted)) +
            labs(caption = paste("Data: Average density at",d[,max(quantile,na.rm=TRUE)],"points"),
                 subtitle = paste("Avg exponential parameter:",round(mean(tab$lambda),5))) +
            ggtitle(paste("Density over time in",ic)) +
            scale_x_continuous(name = "Distance to Center in km")

        ggsave(l[[ic]], filename = file.path(dataplots(),paste0("density-center-",ic,".pdf")))
    }
    l
}

#' plot manual vs satellite measurements 2015/2016
#'
#' makes a plot to compare our manual area measures in 2016 to the satellite measure in 2015
plot_sat_vs_manual <- function(){
    f = readRDS(file.path(LandUseR:::outdir(),"data","france_final.Rds"))
    cf = dcast.data.table(f[year > 2014,list(type, area, LIBGEO)], LIBGEO~ type, value.var = "area")
    cf[ , error := manual - satellite]
    p = ggplot(dcast.data.table(f[year > 2014,list(type, area, LIBGEO)], LIBGEO~ type, value.var = "area"), aes(x= manual, y = satellite)) + geom_point() + geom_abline(slope=1) + ggtitle("manual (2016) vs satellite (2015) area measures", subtitle = "solid line is slope 1") + scale_x_log10("km2 manual 2016") + scale_y_log10("km2 satellite 2015")
    ggsave(p, filename = file.path(LandUseR:::outdir(),"data","plots","areas-manual-satellite.pdf"))
    p
}

#' plot top 100 cities densities over time
#' https://github.com/floswald/LandUse.jl/issues/42
#'
#' @param w plot width in inches
#' @param h plot height in inches
#'
#' @export
plot_top100_densities <- function(save = FALSE,w=9,h=6){
    f = readRDS(file.path(LandUseR:::outdir(),"data","france_final.Rds"))
    f[,density := pop / area]
    ff = f[,list(LIBGEO,year,density,rank,pop)]
    ff = ff[complete.cases(ff)]
    p1 = ggplot(ff, aes(factor(year), y = density)) + geom_violin(draw_quantiles = 0.5) + scale_y_log10(labels = scales::comma, name = "population / km2") + ggtitle("Urban Density over time in France") + scale_x_discrete(name = "year") + theme_bw()

    # label min/max city
    ll = ff[, .SD[,list(mi = min(density),
                        ma = max(density),
                        milab = LIBGEO[which.min(density)],
                        malab = LIBGEO[which.max(density)])],by = factor(year)]
    p2 = p1 + geom_text(data = ll, mapping = aes(x = factor, y = mi, label = milab), nudge_y = -.1)
    p3 = p2 + geom_text(data = ll, mapping = aes(x = factor, y = ma, label = malab),nudge_y = .1)

    p = list()
    p$violin <- p2 + geom_text(data = ll, mapping = aes(x = factor, y = ma, label = malab),nudge_y = .1)


    # aggregate time series
    d4 = ff[,.(density = median(density),wdensity = weighted.mean(density, w = pop)),by=year]
    fall = d4[,round(max(density)/min(density),0)]
    wfall = d4[,round(max(wdensity)/min(wdensity),0)]
    p$ts = ggplot(d4,aes(x=year, y = density)) +
        geom_point(size = 3) +
        scale_x_continuous(breaks = d4[,year]) +
        scale_y_continuous(name = "median density (pop / km2)", breaks = c(3000,5000,10000,25000)) +
        geom_path(size = 1.1) + theme_bw() +
        ggtitle(paste0("Median Urban Density in France fell by Factor ",fall),subtitle = "Top 100 French cities") +
        theme(panel.grid.minor = element_blank())

    p$ts_w = ggplot(d4,aes(x=year, y = wdensity)) +
        geom_point(size = 3) +
        scale_x_continuous(breaks = d4[,year]) +
        scale_y_continuous(name = "weighted mean density (pop / km2)", breaks = c(3000,5000,10000,25000)) +
        geom_path(size = 1.1) + theme_bw() +
        ggtitle(paste0("Mean Urban Density in France fell by Factor ",wfall),subtitle = "Top 100 French cities") +
        theme(panel.grid.minor = element_blank()) + labs(caption = "mean weighted by population")

    p$ts_log = ggplot(d4,aes(x=year, y = density)) +
        geom_point(size = 3) +
        scale_x_continuous(breaks = d4[,year]) +
        scale_y_log10(name = "[log scale] median density (pop / km2)", breaks = c(3000,5000,10000,25000)) +
        geom_path(size = 1.1) + theme_bw() +
        ggtitle(paste0("Median Urban Density in France fell by Factor ",fall),subtitle = "Top 100 French cities") +
        theme(panel.grid.minor = element_blank())

    p$ts_log_w = ggplot(d4,aes(x=year, y = wdensity)) +
        geom_point(size = 3) +
        scale_x_continuous(breaks = d4[,year]) +
        scale_y_log10(name = "[log scale] weighted mean density (pop / km2)", breaks = c(3000,5000,10000,25000)) +
        geom_path(size = 1.1) + theme_bw() +
        ggtitle(paste0("Mean Urban Density in France fell by Factor ",wfall),subtitle = "Top 100 French cities") +
        theme(panel.grid.minor = element_blank()) + labs(caption = "mean weighted by population")



    # time series for tops
    dtnow = ff[LIBGEO %in% LandUseR:::top5now()]
    p$tstop = ggplot(dtnow,aes(x=year, y = density, color = LIBGEO)) +
        geom_line(size = 1.1) +
        scale_y_log10(name = "[log scale] median density (pop / km2)", breaks = c(3000,5000,10000,25000,66000)) +
        labs(color = "City") +
        scale_x_continuous(breaks = dtnow[,year]) +
        ggtitle(paste0("Urban Density by City"),subtitle = "Top 5 cities in 2015") +
        theme_bw() + theme(panel.grid.minor = element_blank())

    dtthen = ff[LIBGEO %in% LandUseR:::top5then()]
    p$tstopthen = ggplot(dtthen,aes(x=year, y = density, color = LIBGEO)) +
        geom_line(size = 1.1) +
        scale_y_log10(name = "[log scale] median density (pop / km2)", breaks = c(3000,5000,10000,25000,66000)) +
        labs(color = "City") +
        scale_x_continuous(breaks = dtthen[,year]) +
        ggtitle(paste0("Urban Density by City"),subtitle = "Top 5 cities in 1866") +
        theme_bw() + theme(panel.grid.minor = element_blank())


    # cross section
    p$xsect <- list()
    p$xsect[[1]] <- ggplot(ff[LIBGEO %in% c("Paris", "Lyon")], aes(x = log(pop), y = density, color = LIBGEO)) + geom_point() + geom_line()
    p$xsect[[2]] <- ggplot(ff[LIBGEO %in% c("Toulouse", "Lille")], aes(x = pop, y = density, color = LIBGEO)) + geom_point() + geom_line()
    p$xsect[[3]] <- ggplot(ff[LIBGEO %in% c("Nantes", "Lyon")], aes(x = pop, y = density, color = LIBGEO)) + geom_point() + geom_line()
    p$xsect[[4]] <- ggplot(ff[LIBGEO %in% c("Toulouse", "Saint-Malo")], aes(x = pop, y = density, color = LIBGEO)) + geom_point() + geom_line()

    # small vs large cities
    ff[, small := rep(.SD[year == 1876, pop < median(pop)], 6)]
    d = ff[, list(pop = median(pop), density = median(density)) , by = list(year,small)]
    p$xsect[[5]] = ggplot(d, aes(x = pop, y = density, color = small)) + geom_point() + geom_line()
    p$xsect[[6]] = ggplot(d, aes(x = pop, y = density, color = factor(year))) + geom_point() + geom_line()
    p$xsect[[7]] =

    if (save) {
        ggsave(p$violin, width = w, height = h,filename = file.path(LandUseR:::outdir(),"data","plots","densities-violins.pdf"))
        ggsave(p$ts, width = w, height = h,filename = file.path(LandUseR:::outdir(),"data","plots","densities-time.pdf"))
        ggsave(p$ts_log, width = w, height = h,filename = file.path(LandUseR:::outdir(),"data","plots","densities-time-log.pdf"))
        ggsave(p$ts_w, width = w, height = h,filename = file.path(LandUseR:::outdir(),"data","plots","densities-time-wtd.pdf"))
        ggsave(p$ts_log_w, width = w, height = h,filename = file.path(LandUseR:::outdir(),"data","plots","densities-time-log-wtd.pdf"))
        ggsave(p$tstop, width = w, height = h,filename = file.path(LandUseR:::outdir(),"data","plots","densities-time-top5now-city.pdf"))
        ggsave(p$tstopthen, width = w, height = h,filename = file.path(LandUseR:::outdir(),"data","plots","densities-time-top5then-city.pdf"))

    }

    p
}

#' Write Paris Arrondissement Areas to Disk
#'
write_paris_areas <- function(){
    sh = sf::st_read(file.path("~/git/intro-to-r/data","CONTOURS-IRIS","1_DONNEES_LIVRAISON_2018-06-00105","CONTOURS-IRIS_2-1_SHP_LAMB93_FXX-2017","CONTOURS-IRIS.shp"),stringsAsFactors=FALSE)
    a=sh %>% filter(INSEE_COM > 75100 & INSEE_COM < 75121) %>% group_by(INSEE_COM) %>% summarise(n=n()) %>% mutate(area = units::set_units(sf::st_area(.), km^2), CODGEO = INSEE_COM) %>% sf::st_set_geometry(NULL) %>% dplyr::select(CODGEO,area)
    saveRDS(a, file.path(LandUseR:::datadir(),"paris-areas.Rds"))
}


#' Plot Density for Central vs Fringe Paris
#'
#' Reads area data from \code{\link{write_paris_areas}}, merges with paris census data
#' and computes densitiy for each arrondissement over time.
plot_paris_densities <- function(){
    a = readRDS(file.path(LandUseR:::datadir(),"paris-areas.Rds"))
    p = readpop() %>%
        filter(DEP == "75") %>%
        left_join(a, by = "CODGEO") %>%
        mutate(density = population / as.numeric(area)) %>%
        mutate(center = if_else(as.integer(CODGEO) < 75107, TRUE, FALSE)) %>%
        group_by(year, center) %>%
        summarise(density = mean(density), population = sum(population))
    out = ggplot(p, aes(year,density,color = center)) +
        geom_point(aes(size = population)) +
        geom_line() + theme_bw() +
        ggtitle("Paris Central vs Fringe Density", subtitle = "Center: Arrondissements 1-6. Constant geography of 2021 Arrondissement boundaries.") +
        labs(caption = "Each point corresponds to a Census counts by arrondissement") +
        scale_y_continuous("population / km2")
    ggsave(plot = out, filename = file.path(LandUseR:::outdir(),"data","plots","paris-densities.pdf"))
    out


}


plot_sat_densities <- function(){
    z = readRDS(file.path(LandUseR:::outdir(),"data","france_final.Rds"))
    x = z[type == "satellite" ,list(year,p50,p90,p10,CODGEO,LIBGEO,rank,area,logarea=log(area),pop,logpop = log(pop))]

    parea = list()
    ppop = list()

    ppop$q10 = ggplot(x[rank < 11],aes(x = logpop, y = p10, color = LIBGEO, label = year, group = LIBGEO)) + geom_point() + geom_path() + geom_text(hjust = 0,nudge_x = 0.05) + ggtitle("10th percentile pop density",subtitle= "how did the lowest density parts of each city develop?")
    ppop$q50 = ggplot(x[rank < 11],aes(x = logpop, y = p50, color = LIBGEO, label = year, group = LIBGEO)) + geom_point() + geom_path() + geom_text(hjust = 0,nudge_x = 0.05) + ggtitle("median pop density",subtitle= "how did the median density parts of each city develop?")
    ppop$q90 = ggplot(x[rank < 11],aes(x = logpop, y = p90, color = LIBGEO, label = year, group = LIBGEO)) + geom_point() + geom_path() + geom_text(hjust = 0,nudge_x = 0.05) + ggtitle("90th percentile pop density",subtitle= "how did the densest parts of each city develop?")
    lapply(ppop, function(x) theme_set(theme_bw()))


    parea$q10 = ggplot(x[rank < 11],aes(x = logarea, y = p10, color = LIBGEO, label = year, group = LIBGEO)) + geom_point() + geom_path() + geom_text(hjust = 0,nudge_x = 0.05) + ggtitle("10th percentile pop density",subtitle= "how did the lowest density parts of each city develop?")
    parea$q50 = ggplot(x[rank < 11],aes(x = logarea, y = p50, color = LIBGEO, label = year, group = LIBGEO)) + geom_point() + geom_path() + geom_text(hjust = 0,nudge_x = 0.05) + ggtitle("median pop density",subtitle= "how did the median density parts of each city develop?")
    parea$q90 = ggplot(x[rank < 11],aes(x = logarea, y = p90, color = LIBGEO, label = year, group = LIBGEO)) + geom_point() + geom_path() + geom_text(hjust = 0,nudge_x = 0.05) + ggtitle("90th percentile pop density",subtitle= "how did the densest parts of each city develop?")
    lapply(parea, function(x) theme_set(theme_bw()))


    setkey(x,CODGEO,year)
    x[ , agrowth_pct := .SD[,100*(area[.N] - area[1])/area[1]],by = CODGEO]
    x[ , agrowth := .SD[,area[.N] - area[1]],by = CODGEO]
    x[ , pgrowth := .SD[,pop[.N] - pop[1]],by = CODGEO]
    x[ , pgrowth_pct := .SD[,100*(pop[.N] - pop[1])/pop[1]],by = CODGEO]

    # summarize: percent growth for all cols
    s = x[ , lapply(.SD,function(y) 100*(y[.N] - y[1])/y[1]), .SDcols = c("p10","p50","p90", "area","pop"),by = .(CODGEO,LIBGEO)]

    # growth plots
    pg = list()

    pg$p10 = ggplot(s,aes(x=area,y=p10)) + geom_point() + scale_x_continuous("% change in area") + scale_y_continuous("% change in density") + geom_smooth(method = "lm") + ggtitle("% change in density at 10-th quantile", subtitle = "Change 1975-2015")
    pg$p50 = ggplot(s,aes(x=area,y=p50)) + geom_point() + scale_x_continuous("% change in area") + scale_y_continuous("% change in density")+ geom_smooth(method = "lm") + ggtitle("% change in density at 50-th quantile", subtitle = "Change 1975-2015")
    pg$p90 = ggplot(s,aes(x=area,y=p90)) + geom_point() + scale_x_continuous("% change in area") + scale_y_continuous("% change in density")+ geom_smooth(method = "lm") + ggtitle("% change in density at 90-th quantile", subtitle = "Change 1975-2015")

    lapply(pg, function(x) theme_set(theme_bw()))
    return(list(growth = pg, area = parea, pop = ppop))

}


plot_paris <- function(d){
    pl = d %>%
        dplyr::filter(DEP == 75) %>%
        dplyr::mutate(Arrond = as.factor(as.numeric(CODGEO)-75100)) %>%
        ggplot2::ggplot(aes(x = year, y = population, group = CODGEO, color = Arrond)) + ggplot2::geom_line(size=1.1) + ggplot2::scale_y_continuous(labels = scales::comma) + ggplot2::labs(title = "Central Paris Population by Arrondissement",subtitle = "Equivalent to `Density` (Area is constant)") + ggplot2::theme_bw()
    ggplot2::ggsave(file.path(outdir(),"plots","paris.pdf"),plot = pl, width = 8, height = 6)
    pl

}

