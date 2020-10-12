

plot_sat_built <- function(){
    cc = LandUseR:::cropcities()  # first index level is years

    ci = bboxes_top100()[LIBGEO %in% tops()]  # 5 cities we want to see
    out = list()
    for (ir in 1:nrow(ci)){
        cci = ci[ir,]
        # indices 1,2,3,4 are years
        ss = stack(cc[[1]][[cci$CODGEO]]$built,cc[[2]][[cci$CODGEO]]$built,cc[[3]][[cci$CODGEO]]$built,cc[[4]][[cci$CODGEO]]$built)
        pngname = file.path(LandUseR:::outdir(),"data","plots",paste0(cci$LIBGEO,"-sat-area.png"))
        pdfname = file.path(LandUseR:::outdir(),"data","plots",paste0(cci$LIBGEO,"-sat-area.pdf"))

        # plot
        out[[ir]] = levelplot(ss,xlab = NULL, ylab = NULL, scales=list(draw=FALSE), names.attr=paste0(cci$LIBGEO,' % built ', LandUseR:::GHS_years()))

        # write
        png(pngname,width = 1000, height=700, res = 175)
        print(out[[ir]])
        dev.off()
        pdf(pdfname,width = 11,height = 6)
        print(out[[ir]])
        dev.off()


    }
    out
}

#' plot manual vs satellite measurements 2015/2016
#'
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
#'@export
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

