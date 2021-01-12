# resource: https://geocompr.robinlovelace.net
# resource: https://datacarpentry.org/r-raster-vector-geospatial/

# This uses produced by the European Commission for their Global Human Settlement
# database GHS.
# https://ghslsys.jrc.ec.europa.eu/ghs_smod2019.php


# urban centre database
# https://ghslsys.jrc.ec.europa.eu/ghs_stat_ucdb2015mt_r2019a.php
# UCDB <- function(overwrite = FALSE){
#   if (overwrite){
#     x = sf::st_read(file.path(datadir(),"GHS/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp"))
#     saveRDS(x,file.path(datadir(),"UCDB.Rds"))
#   } else {
#     x = sf::st_as_sf(readRDS(file.path(datadir(),"UCDB.Rds")))
#   }
#   x
# }

GHS_years <- function() c(1975,1990,2000,2015)


cropcities <- function(overwrite = FALSE){
  if (overwrite){
    x = LandUseR:::pop_1950_2()
    lgs = x[,unique(CODGEO)]
    L = list()
    for (iy in 1:length(GHS_years())){

      L[[iy]] <- list()

      yr = GHS_years()[iy]
      flog.info("cropping year %d",yr)
      r = LandUseR:::loadRasters(yr)
      pb <- progress::progress_bar$new(total = 100)
      for (lg in lgs){

        L[[iy]][[lg]] = crop_rasters(r,x[CODGEO==lg,extent][[1]])
        pb$tick()
      }

    }
    names(L) <- paste(GHS_years())
    saveRDS(L,file.path(outdir(),"data","france_cropped.Rds"))
  } else {
    L = readRDS(file.path(outdir(),"data","france_cropped.Rds"))
  }
  L
}

#' take each cities cropped raster and
#' measure total pop and built up
measure_cities <- function(overwrite = FALSE, cutoff = 30){
    if (overwrite){
        # get list of cities
        L = readRDS(file.path(LandUseR:::outdir(),"data","france_cropped.Rds"))

        # output list
        OL = list()

        # for each year and each city
        for (iy in 1:length(L)){
            yr = GHS_years()[iy]
            flog.info("measuring cities in year %d",yr)

            OL[[iy]] <- list()

            for (ic in names((L[[iy]]))){
                # cut out all cells above cutoff
                cutmask = copy(L[[iy]][[ic]]$built)
                raster::values(cutmask) <- (raster::values(L[[iy]][[ic]]$built) > cutoff)
                cityb = L[[iy]][[ic]]$built[ cutmask , drop = FALSE]   # built cut out of city
                cityp = L[[iy]][[ic]]$pop[ cutmask , drop = FALSE]   # pop cut ouf of city (same cut out)

                # get biggest connected set of built
                cc = raster::clump(cityb)
                fc = data.frame(raster::freq(cc))
                whichc = fc %>% na.omit() %>% dplyr::filter(count == max(count)) %>% dplyr::pull(value)  # get biggest cluster
                citymask = copy(cc)
                raster::values(citymask) <- (raster::values(cc) == whichc)
                raster::values(citymask)[!raster::values(citymask)] = NA  # set FALSE to NA

                # inverse mask
                # imask = citymask
                # raster::values(imask) <- is.na(raster::values(citymask))
                # raster::values(imask)[!raster::values(imask)] <- NA

                totarea = raster::cellStats(citymask, "sum",na.rm = T)   # sum of gridcells
                totarea = totarea * (250 * 250) / 1e+6  # total area in square kilometers

                totpop = sum(cityp[citymask],na.rm = T)
                meanbuilt = mean(cityb[citymask],na.rm = T)
                qbuilt = quantile(cityb[citymask],na.rm = T)
                qpop = quantile(cityp[citymask],na.rm = T)

                # output
                OL[[iy]][[ic]] = list(area = totarea, pop = totpop, meanbuilt = meanbuilt,
                                      b10 = qbuilt[1],
                                      b25 = qbuilt[2],
                                      b50 = qbuilt[3],
                                      b75 = qbuilt[4],
                                      b90 = qbuilt[5],
                                      p10 = qpop[1],
                                      p25 = qpop[2],
                                      p50 = qpop[3],
                                      p75 = qpop[4],
                                      p90 = qpop[5])
            }
        }
        names(OL) <- names(L)
        saveRDS(OL,file.path(LandUseR:::outdir(),"data","france_measured.Rds"))

    } else {
        OL = readRDS(file.path(LandUseR:::outdir(),"data","france_measured.Rds"))
    }
    return(OL)
}


combine_measures <- function(overwrite = FALSE,cutoff = 30){
    # melt manual measures
    x = LandUseR:::pop_1950_2(overwrite = overwrite)
    setnames(x,"area_EM","area_1876")
    m = melt.data.table(x,id = c("CODGEO","LIBGEO","rank","REG","DEP","extent"),
                      measure = patterns("^pop_","^area_"),
                      value.name = c("pop","area"),
                      variable.name = "iyear",variable.factor = FALSE)
    m[iyear == 1 , year := 1876]
    m[iyear == 2 , year := 1950]
    m[iyear == 3 , year := 2016]
    m[ , type := "manual"]
    m[,iyear := NULL]

    # now add satelite measures
    m2 = m[rep(1:100,4)]
    m2[,year := rep(LandUseR:::GHS_years(),each = 100)]
    m2[,c("pop","area") := NULL]
    setkey(m2,year,CODGEO)

    sat = LandUseR:::measure_cities(overwrite = overwrite, cutoff = cutoff)
    # fill in values
    ll = list()
    for (yr in GHS_years()){
        ss = sat[[paste(yr)]]
        q = rbindlist(ss)
        q[,c("year","CODGEO") := list(yr,names(ss))]
        ll[[yr]] <- merge(m2[.(yr)], q, by = c("year","CODGEO"))
    }
    m2 = rbindlist(ll)
    m2[,type := "satellite"]
    z = rbind(m,m2,fill = TRUE)
    setkey(z, year, rank)
    saveRDS(z, file.path(LandUseR:::outdir(),"data","france_final.Rds") )
    fwrite(z[,!"extent"],  file.path(LandUseR:::outdir(),"data","france_final.csv"))
    z
}


tops <- function(){
    c("Paris", "Lyon", "Marseille", "Toulouse", "Reims")
}
top5then <- function(){
  c("Paris", "Lyon", "Marseille", "Bordeaux", "Lille")
}
top5now <- function(){
  c("Paris", "Lyon", "Lille", "Marseille", "Nice")
}


fun = function(){
  r = LandUseR:::loadRasters(1975)
  x = LandUseR:::pop_1950_2()
  z = LandUseR:::crop_rasters(r,x[LIBGEO=="Brest",extent][[1]])
}

#' crop a list of rasters to new extent
#'
#' takes list r of rasters, records CRS of `pop` entry
#' and crops to new extent `e`
crop_rasters <- function(r,e){
    stopifnot(is(e,"Extent"))
    dr = raster::raster()
    raster::extent(dr) <- e
    ee = raster::extent(raster::projectExtent(dr,crs = raster::projection(r$pop)))
    y = lapply(r, function(x) raster::crop(x,ee))
    y
}


#' Cut Out rasters for French cities
#'
# cropFrance <- function(cutoff=20,overwrite = FALSE,overwrite_cropping = FALSE){
#     if (overwrite) {
#       bb <- bboxes()   # get bounding boxes for top 30 cities
#       res = bb[,list(city,year = 1975, built = NA_real_, pop = NA_real_)]
#       res = res[rep(1:nrow(bb),4)]
#       res[,year := rep(GHS_years(),each = nrow(bb))]
#       setkey(res,city,year)
#
#       for (yr in GHS_years()){
#         # load raster for that year
#         r = loadRaster(yr = yr)
#         for (ci in french_cities()){
#         # for (ci in c("Paris [FRA]","Lyon [FRA]")){
#           x = crop_smod(bb[city == ci],r,yr,cutoff = cutoff,overwrite=overwrite_cropping)
#           res[.(ci,yr), c("built", "pop") := list(x$area,x$pop)]
#           # print(res[.(ci)])
#         }
#       }
#       saveRDS(res,file.path(outdir(),"france_cropped.Rds"))
#
#     } else {
#       res = readRDS(file.path(outdir(),"france_cropped.Rds"))
#     }
#
#     res
# }



#' Crop GHS Settlement Model for One City
#'
#' @param yr for which year to do that
# crop_smod <- function(bbox,r,yr,cutoff = 20,overwrite = FALSE){
#     city = bbox[,city]
#     flog.info("cropping %s to level %s in year %d",city,cutoff,yr)
#
#     # get bouding box for city
#     # pa = LandUseR:::paris_box()
#     # pa = LandUseR:::city_box(city)
#
#
#     # setup output dir structure
#     odir = file.path(outdir(),city)
#     dir.create(odir,showWarnings = FALSE, recursive = T)
#     levs = GHS_L2()
#     ghc_cols = GHS_LS_scale()
#
#     if (!file.exists(file.path(odir,paste0("cropped",yr,".Rds"))) | overwrite){
#         flog.info("re-cropping raster to %s",city)
#
#         r_c <- crop_rasters(r,bbox[,extent][[1]])  # extents in WGS84!
#         # r_c$smod <- raster::as.factor(r_c$smod)
#         # print(unique(values(r_c$smod)))
#         r_c$builtpop <- raster::stack(r_c$built,r_c$pop)
#         names(r_c$builtpop) <- c("built","pop")
#         r_c$built <- NULL
#         r_c$pop <- NULL
#
#         saveRDS(r_c, file = file.path(odir,paste0("cropped",yr,".Rds")))
#     } else {
#         r_c <- readRDS(file.path(odir,paste0("cropped",yr,".Rds")))
#     }
#
#     r_c$smod <- raster::ratify(r_c$smod)   # need to preserve categorical data!!!! no reprojection possible!!!!
#     # rat <- levels(r_c$smod)[[1]]
#     # rat$ID <- LandUseR:::GHS_L2()$ID
#     # rat$label <- as.character(LandUseR:::GHS_L2()$label)
#     #
#     # # rat$cols <- LandUseR:::GHS_LS_scale()
#     # levels(r_c$smod) <- rat
#
#     levs = subset(levs, ID %in% unique(raster::values(r_c$smod)))
#     ghc_cols = subset(ghc_cols, ID %in% unique(raster::values(r_c$smod)))
#
#     # png(file.path(odir,paste0("SMOD_",city,yr,".png")),width = 9,height=6,units = "in",res = 200)
#     # plot(r_c$smod,legend = FALSE, col = cols$color,main=paste(city,"SMOD as is",yr))
#     # legend(x = 'bottomright', legend = levs$label, fill = cols$color)
#     # dev.off()
#
#     # plot fine grids
#     png(file.path(odir,paste0("pop-built-",city,yr,".png")),width = 10,height=5,units = "in",res = 200)
#     par(mfcol = c(1,2))
#     raster::plot(r_c$builtpop$pop,xaxt="n",yaxt="n", legend = FALSE, paste(city,"population"))
#     raster::plot(r_c$builtpop$pop, legend.only=TRUE, horizontal = TRUE,
#          legend.args = list(text='Population', side = 1, line = 2))
#     raster::plot(r_c$builtpop$built,xaxt="n",yaxt="n", legend = FALSE, paste(city,"% Built"))
#     raster::plot(r_c$builtpop$built, legend.only=TRUE, horizontal = TRUE,
#                  legend.args = list(text='% Built', side = 1, line = 2))
#     par(mfcol = c(1,1))
#     dev.off()
#
#
#
#     # aggregate pop and built to 1km grids
#     built1k <- raster::aggregate(r_c$builtpop$built, fact = 4, mean)
#     pop1k <- raster::aggregate(r_c$builtpop$pop, fact = 4, sum)
#     bp1k <- raster::stack(built1k, pop1k)
#     names(bp1k) <- c("built","pop")
#
#     # force extents to be equal
#     raster::extent(bp1k) <- raster::extent(r_c$smod)
#
#     # get definition of SMOD "city"
#     s_mask = r_c$smod
#     raster::values(s_mask) = raster::values(r_c$smod > cutoff)
#     stopifnot(raster::extent(s_mask) == raster::extent(bp1k))
#
#     # clip to mask
#     # bp1k_c = bp1k[s_mask, drop = FALSE]
#
#     cc <- r_c$smod[s_mask,drop = FALSE]
#     cc = raster::clump(cc)
#     fc = data.frame(raster::freq(cc))
#     whichc = fc %>% na.omit() %>% dplyr::filter(count == max(count)) %>% dplyr::pull(value)
#     mask_city = cc
#     raster::values(mask_city) = raster::values(cc == whichc)
#     raster::values(mask_city)[raster::values(mask_city) == 1] = TRUE
#     raster::values(mask_city)[raster::values(mask_city) != 1] = NA
#     # values(mask_city) = values(cc > 0)
#     # values(mask_city)[values(mask_city) == 1] = TRUE
#     # values(mask_city)[values(mask_city) != 1] = NA
#
#
#     # inverse mask
#     imask = mask_city
#     raster::values(imask) <- is.na(raster::values(mask_city))
#     raster::values(imask)[!raster::values(imask)] <- NA
#
#     totarea = raster::cellStats(mask_city, sum)   # total area of main city from SMOD
#     bp1k_mask = bp1k[s_mask]
#     totpop = sum(bp1k_mask[,"pop"],na.rm=T)
#     meanbuilt = mean(bp1k_mask[,"built"],na.rm=T)
#
#     # plot city model only
#     png(file.path(odir,paste0("SMOD_",city,yr,".png")),width = 9,height=6,units = "in",res = 200)
#     raster::plot(r_c$smod,legend = FALSE, col = ghc_cols$color,main=paste(city,"Settlement Model",yr))
#     legend(x = 'bottomright', legend = levs$label, fill = ghc_cols$color)
#     dev.off()
#
#
#     # plot city model and mask
#     png(file.path(odir,paste0("SMOD_",city,"_classify",yr,".png")),width = 9,height=6,units = "in",res = 200)
#     par(mfcol = c(1,2))
#     raster::plot(r_c$smod,legend = FALSE, col = ghc_cols$color,main=paste(city,"SMOD as is",yr))
#     legend(x = 'bottomright', legend = levs$label, fill = ghc_cols$color)
#     raster::plot(r_c$smod,legend = FALSE, col = ghc_cols$color,main = "Main City Classification",alpha = 0.1)
#     raster::plot(mask_city,legend = FALSE, col = ghc_cols$color[length(ghc_cols$color)], add = TRUE)
#     dev.off()
#     par(mfcol = c(1,1))
#
#     # plot both original and masked data
#     png(file.path(odir,paste0("maskpop-",city,yr,".png")),width = 9,height=6,units = "in",res = 200)
#     par(mfcol = c(1,2))
#     raster::plot(bp1k$pop,legend = TRUE,main=paste(city,"pop as is",yr))
#     # plot(bp1k$pop, legend.only=TRUE, horizontal = TRUE,
#     #      legend.args = list(text='Population', side = 1, line = 2))
#     raster::plot(bp1k$pop, main = "Cut-out population",legend = FALSE)
#     raster::plot(imask,add=TRUE,col="black",legend=FALSE)
#     dev.off()
#     par(mfcol = c(1,1))
#
#     png(file.path(odir,paste0("maskbuilt-",city,yr,".png")),width = 9,height=6,units = "in",res = 200)
#     par(mfcol = c(1,2))
#     raster::plot(bp1k$built,legend = TRUE,main=paste(city,"Built as is",yr))
#     # plot(bp1k$pop, legend.only=TRUE, horizontal = TRUE,
#     #      legend.args = list(text='Population', side = 1, line = 2))
#     raster::plot(bp1k$built, main = "Cut-out Built",legend = FALSE)
#     raster::plot(imask,add=TRUE,col="black",legend=FALSE)
#     dev.off()
#     par(mfcol = c(1,1))
#
#     # get densitty in each gridcell
#     # dgrid = copy(built_paris)
#     # values(dgrid) <- values(pop_paris) / values(built_paris)
#
#     l <- list(area = totarea, pop = totpop, meanbuilt=meanbuilt)
#     return(l)
# }


#' France UCDB: Cities
#'
# franceUCDB <- function(topn = 30,overwrite = FALSE){
#   if (overwrite) {
#     x = UCDB()
#     # WTF is Cambrin?? https://en.wikipedia.org/wiki/Cambrin
#     y = x %>% dplyr::filter(CTR_MN_ISO == "FRA") %>%
#       sf::st_set_geometry(NULL) %>% data.table
#
#     y[, city := as.character(UC_NM_MN)]
#     y = y[!(city %in% c("Cambrin [FRA]","Brest [FRA]","Nantes [FRA]", "Rennes [FRA]","Saint-Denis [FRA]")), .(P15,P00, P90, P75, B15, B00, B90, B75, city)][order(P15,decreasing = T)][1:topn]
#
#     setcolorder(y,c("city","P75","P90","P00","P15","B75","B90","B00","B15"))
#     m = melt.data.table(y, id = "city", measure =patterns("^P","^B"), value.name = c("pop","built"), variable.name = "nyear", variable.factor = FALSE)
#     m[nyear == "4", year := 2015]
#     m[nyear == "3", year := 2000]
#     m[nyear == "2", year := 1990]
#     m[nyear == "1", year := 1975]
#     m[,nyear := NULL]
#     m[,density := pop / built]
#     m[city == "Clermont-Ferran [FRA]", city := "Clermont-Ferrand [FRA]"]
#     setcolorder(m,c("city","year"))
#     fwrite(m,file.path(datadir(),"france_cities_ucdb.csv"))
#   } else {
#     if (!file.exists(file.path(datadir(),"france_cities_ucdb.csv"))) {
#       stop("must run with overwrite=TRUE first")
#     } else {
#       m = fread(file.path(datadir(),"france_cities_ucdb.csv"))
#     }
#   }
#   m
# }

#' List of French cities in study
#'
#' @export
french_cities <- function(){
  f = france()
  unique(f$city)
}


# plot_france <- function(y,outcome = "density",fname = "densities-1975-FRA",dotted = FALSE){
#   others = "Paris|Lyon|Marseille|Montpellier"
#   y[city %like% others,lines :=  "Paris"]
#   y[!(city %like% others),lines :=  "Other"]
#   years = y[,unique(year)]
#   setkey(y,city,year)
#   topn = y[,length(unique(city))]
#   linet = grep(others,y[,unique(city)]) # indices of dashed lines
#   linev = rep("solid",topn)
#   if (dotted){
#     linev[linet] <- "dotted"
#   }
#
#
#   pl = ggplot2::ggplot(y,ggplot2::aes(x = year,y = .data[[outcome]],color = city,linetype = city)) +
#     ggplot2::geom_line() +
#     ggplot2::geom_point() +
#     ggplot2::theme_bw() +
#     ggplot2::scale_color_hue(name = "City") +
#     ggplot2::scale_linetype_manual(name = "City", values = linev) +
#     ggplot2::scale_x_continuous(breaks = years) +
#     ggplot2::scale_y_continuous(name = outcome, labels = scales::comma)
#   if (outcome == "density") {
#     pl = pl + ggplot2::labs(title = paste0("France ",outcome," over time"),
#                            caption = "density = number of people / total built up surface (in km2)")
#   } else {
#     pl = pl + ggplot2::labs(title = paste0("France ",outcome," over time"))
#   }
#
#   ggplot2::ggsave(plot = pl, filename = file.path(outdir(),"plots",paste0(fname,".pdf")),width = 10, height = 6)
#   ggplot2::ggsave(plot = pl, filename = file.path(outdir(),"plots",paste0(fname,".png")),width = 10, height = 6)
#   pl
# }


#' Get French Boudning Boxes
#'
#' @export
bboxes <- function(){
    bb = LandUseR:::france()
    bb <- bb[year == 1975]

    bn = data.table(city = bb[,city], extent = vector("list", length = nrow(bb)))

    bn[city == "Lille [FRA]; Mouscron [BEL]",  extent := raster::extent(c(2.944861, 3.298341 , 50.543319, 50.778156))]
    bn[city == "Lille [FRA]; Mouscron [BEL]",  extent := raster::extent(c(2.2145981, 2.414104 ,49.847181 , 49.959225))]
    bn[city == "Amiens [FRA]",                 extent := raster::extent(c(2.2145981, 2.414104 ,49.847181 , 49.959225))]
    bn[city == "Le Havre [FRA]",               extent := raster::extent(c(0.033290, 0.362656 ,49.454226 , 49.558283))]
    bn[city == "Rouen [FRA]",                  extent := raster::extent(c(0.940602, 1.183368 ,49.283603 , 49.520694))]
    bn[city == "Reims [FRA]",                  extent := raster::extent(c(3.943876, 4.132186 ,49.193367 , 49.310117))]
    bn[city == "Caen [FRA]",                   extent := raster::extent(c(-0.489701, -0.270721 ,49.142144 , 49.243900))]
    bn[city == "Metz [FRA]",                   extent := raster::extent(c(6.059380, 6.270723 ,49.055194 , 49.175658))]
    bn[city == "Nancy [FRA]",                  extent := raster::extent(c(6.038759, 6.395210 ,48.574603 , 48.781997))]
    bn[city == "Paris [FRA]",                  extent := raster::extent(c(1.61, 3.08 ,48.46 , 49.118))]
    bn[city == "Strasbourg [FRA]; Kehl [DEU]", extent := raster::extent(c(7.642639,7.853926 ,48.484619 , 48.672193))]
    # bn[city == "Brest [FRA]",                extent := raster::extent(c(-4.601281,-4.360628 ,48.352529 , 48.468810))]
    bn[city == "Dunkerque [FRA]",              extent := raster::extent(c(2.179341,2.517408, 50.979502, 51.056436))]
    bn[city == "Perpignan [FRA]",              extent := raster::extent(c(2.799684, 2.989531, 42.629930,42.761191))]
    # bn[city == "Rennes [FRA]",               extent := raster::extent(c(-1.778073,-1.574094 ,48.046103 , 48.173615))]
    bn[city == "Le Mans [FRA]",                extent := raster::extent(c(0.091035,0.281172 ,47.917295 , 48.062800))]
    bn[city == "Orleans [FRA]",                extent := raster::extent(c(1.767287,2.073181 ,47.796967 , 47.971889))]
    bn[city == "Mulhouse [FRA]",               extent := raster::extent(c(7.204776,7.403701 ,47.697986 ,47.842130))]
    bn[city == "Angers [FRA]",                 extent := raster::extent(c(-0.655981,-0.447919 ,47.404412 ,47.544447))]
    bn[city == "Tours [FRA]",                  extent := raster::extent(c(0.591164,0.798459 ,47.325082 ,47.483604))]
    bn[city == "Dijon [FRA]",                  extent := raster::extent(c(4.951104,5.166130 ,47.230547 ,47.384402))]
    # bn[city == "Nantes [FRA]",               extent := raster::extent(c(-1.691281,-1.419254 ,47.131828 ,47.309569))]
    bn[city == "Nantes [FRA]",                 extent := raster::extent(c(-1.691281,-1.419254 ,47.131828 ,47.309569))]
    bn[city == "Limoges [FRA]",                extent := raster::extent(c(1.145045,1.358150 ,45.771763 ,45.909983))]
    bn[city == "Clermont-Ferrand [FRA]",       extent := raster::extent(c(3.000352,3.223957 ,45.701680 ,45.847506))]
    bn[city == "Lyon [FRA]",                   extent := raster::extent(c(4.682204,5.052525 ,45.628088 ,45.836050))]
    bn[city == "Saint-Etienne [FRA]",          extent := raster::extent(c(4.215847,4.480892 ,45.369009 ,45.503158))]
    bn[city == "Grenoble [FRA]",               extent := raster::extent(c(5.605581,5.887174 ,45.068137 ,45.261293))]
    bn[city == "Bordeaux [FRA]",               extent := raster::extent(c(-0.830520,-0.377359 ,44.671694 ,44.967807))]
    bn[city == "Montpellier [FRA]",            extent := raster::extent(c(3.649175,4.095652 ,43.447040 ,43.701412))]
    bn[city == "Nice-Cannes [FRA]",            extent := raster::extent(c(6.861982,7.399601 ,43.522038 ,43.764086))]
    bn[city == "Toulouse [FRA]",               extent := raster::extent(c(1.216777,1.655748 ,43.397305 ,43.738317))]
    bn[city == "Marseille [FRA]",              extent := raster::extent(c(5.256408,5.606574 ,43.221860 ,43.393227))]
    bn[city == "Toulon [FRA]",                 extent := raster::extent(c(5.746459,6.189482 ,43.051562 ,43.208111))]
    bn[city == "Pau [FRA]",                    extent := raster::extent(c(-0.483260,-0.273838,43.276050 ,43.391593))]

    bn
}


# plot_france_UCDB <- function(topn = 30){
#   y = france(topn)
#   plot_france(y)
# }

# overlay_GHSPOP_UDB <- function(){
#   # paris city database
#   x = UCDB()
#   pa = x %>% dplyr::filter(grepl("Paris",UC_NM_MN))
#   ma = x %>% dplyr::filter(grepl("Marseille",UC_NM_MN))
#
#   for (yr in paste(c(1975, 1990, 2000, 2015))){
#
#     pop <- raster::raster(file.path(datadir(),paste0("GHS/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3.tif")))
#     if(yr==2015){
#       built <- raster::raster(file.path(datadir(),paste0("GHS/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3.tif")))
#     } else {
#       built <- raster::raster(file.path(datadir(),paste0("GHS/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3.tif")))
#     }
#
#
#
#     # same crs
#     pat = sf::st_transform(pa, raster::projection(pop))
#
#     # crop raster to shape of border
#     pop_paris = raster::crop(pop,raster::extent(pat))
#     built_paris = raster::crop(built,raster::extent(pat))
#
#     # plot
#     png(file.path(outdir(),paste0("plots/UCDB-paris-",yr,".png")),width = 19,height=8,units = "in",res = 200)
#     par(mfcol=c(1,2))
#     raster::plot(pop_paris,main = paste0("Paris Population ",yr," overlaid with city boundary 2015"), legend=FALSE, axes=FALSE)
#     raster::plot(pop_paris,legend.only=TRUE,horizontal = TRUE,
#                                 legend.args=list(text='Number of People',side=1,line = 2))
#     plot(sf::st_geometry(pat),add = TRUE, fill = NA, border = "red", lw = 2)
#     raster::plot(built_paris,main = paste0("Paris Built ",yr), legend=FALSE, axes=FALSE)
#     raster::plot(built_paris,legend.only=TRUE,horizontal = TRUE,
#                  legend.args=list(text='Percent Built Up',side=1,line = 2))
#     plot(sf::st_geometry(pat),add = TRUE, fill = NA, border = "red", lw = 2)
#     dev.off()
#   }
#   # https://cran.r-project.org/web/packages/magick/vignettes/intro.html
#   img <- magick::image_read(grep("UCDB-paris-",list.files(file.path(outdir(),"plots"),full.names = TRUE),value = TRUE))
#   # ani <- magick::image_animate(img,fps = 0.5)
#   magick::image_write_gif(img, file.path(outdir(),paste0("plots/UCDB-paris.gif")),delay = 2)
#   # magick::image_convert(grep("UCDB-paris-",list.files(file.path(outdir(),"plots"),full.names = TRUE),value = TRUE))
# }
#
# GHS_L2 <- function(){
#     data.frame(ID = c(10,11,12,13,21,22,23,30), label = c("Water","Uninhabited","Rural Dispersed","Village","Suburbs","Semi-Dense Town","Dense Town","City"))
# }
# GHS_LS_scale <- function(){
#     d = data.frame(ID = c(10,11,12,13,21,22,23,30), color =
#     c(rgb(122, 182, 245, maxColorValue = 255),
#       rgb(205, 245, 122, maxColorValue = 255),
#       rgb(171, 205, 102, maxColorValue = 255),
#       rgb(55, 86, 35, maxColorValue = 255),
#       rgb(255, 255, 0, maxColorValue = 255),
#       rgb(168, 112, 0, maxColorValue = 255),
#       rgb(115, 38, 0, maxColorValue = 255),
#       rgb(255, 0, 0, maxColorValue = 255)
#       )
#     )
#     d$color = as.character(d$color)
#     d
# }
#
# rayt <- function(){
#     el = matrix(raster::extract(pop_paris, raster::extent(pop_paris)),nrow = ncol(pop_paris), ncol = nrow(pop_paris))
#     el %>% sphere_shade(texture = "imhof4") %>%   add_shadow(ray_shade(el, zscale = 15, maxsearch = 300), 0.5) %>%  add_shadow(ambmat, 0.5) %>% plot_3d(el, zscale = 10, fov = 0, theta = 135, zoom = 0.75, phi = 45, windowsize = c(1000, 800))
#
#
#     phivechalf = 30 + 60 * 1/(1 + exp(seq(-7, 20, length.out = 180)/2))
#     phivecfull = c(phivechalf, rev(phivechalf))
#     thetavec = -90 + 60 * sin(seq(0,359,length.out = 360) * pi/180)
#     zoomvec = 0.45 + 0.2 * 1/(1 + exp(seq(-5, 20, length.out = 180)))
#     zoomvecfull = c(zoomvec, rev(zoomvec))
#
#     render_movie(filename = "paris-movie", type = "custom",
#                              frames = 360,  phi = phivecfull, zoom = zoomvecfull, theta = thetavec)
# }



loadRasters <- function(yr,WGS84 = FALSE){
  stopifnot(yr %in% c(1975,1990,2000,2015))
  li = list()
  if (WGS84){
      ending = "_WGS84.tif"
  } else {
      ending = ".tif"
  }
  # li$smod <- raster::raster(file.path(LandUseR:::datadir(),paste0("GHS/GHS_SMOD_POP",yr,"_GLOBE_R2019A_54009_1K_V1_0_18_3/GHS_SMOD_POP",yr,"_GLOBE_R2019A_54009_1K_V1_0_18_3",ending)))
  # li$smod <- raster::ratify(li$smod)   # need to preserve categorical data!!!! no reprojection possible!!!!
  # rat <- raster::levels(li$smod)[[1]]
  # rat$cat <- LandUseR:::GHS_L2()
  # # rat$cols <- LandUseR:::GHS_LS_scale()
  # levels(li$smod) <- rat

  pop <- raster::raster(file.path(LandUseR:::datadir(),paste0("GHS/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3",ending)))
  pop2 <- raster::raster(file.path(LandUseR:::datadir(),paste0("GHS/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_17_3/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_17_3",ending)))
  li$pop <- raster::merge(pop,raster::crop(pop2,raster::extent(-500000, -41000, 5e+06, 6e+06)))  # the crop() cuts off 3/4 of that tile, leaving only most eastern strip
  if(yr==2015){
    built <- raster::raster(file.path(LandUseR:::datadir(),paste0("GHS/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3",ending)))
    built2 <- raster::raster(file.path(LandUseR:::datadir(),paste0("GHS/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_3/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_3",ending)))
    li$built <- raster::merge(built,raster::crop(built2,raster::extent(-500000, -41000, 5e+06, 6e+06)))
  } else {
    built <- raster::raster(file.path(datadir(),paste0("GHS/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3",ending)))
    built2 <- raster::raster(file.path(datadir(),paste0("GHS/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_17_3/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_17_3",ending)))
    li$built <- raster::merge(built,raster::crop(built2,raster::extent(-500000, -41000, 5e+06, 6e+06)))

  }
  li
}

reproject <- function(to = "+proj=longlat +datum=WGS84 +no_defs"){
  for (yr in c(1975,1990,2000,2015)){
    futile.logger::flog.info("reprojecting rasters for year %d",yr)

    x = raster::raster(file.path(LandUseR:::datadir(),paste0("GHS/GHS_SMOD_POP",yr,"_GLOBE_R2019A_54009_1K_V1_0_18_3/GHS_SMOD_POP",yr,"_GLOBE_R2019A_54009_1K_V1_0_18_3.tif")))
    y = raster::projectRaster(x, crs = to)
    writeRaster(y, filename = file.path(LandUseR:::datadir(),paste0("GHS/GHS_SMOD_POP",yr,"_GLOBE_R2019A_54009_1K_V1_0_18_3/GHS_SMOD_POP",yr,"_GLOBE_R2019A_54009_1K_V1_0_18_3_WGS84")),format = "GTiff",overwrite = TRUE)

    x = raster::raster(file.path(datadir(),paste0("GHS/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3.tif")))
    y = raster::projectRaster(x, crs = to)
    writeRaster(y, filename = file.path(datadir(),paste0("GHS/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3_WGS84")),format = "GTiff",overwrite = TRUE)

    if(yr==2015){
      x <- raster::raster(file.path(datadir(),paste0("GHS/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3.tif")))
      y = raster::projectRaster(x, crs = to)
      writeRaster(y, filename = file.path(datadir(),paste0("GHS/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3_WGS84")),format = "GTiff",overwrite = TRUE)

    } else {
      x <- raster::raster(file.path(datadir(),paste0("GHS/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3.tif")))
      y = raster::projectRaster(x, crs = to)
      writeRaster(y, filename = file.path(datadir(),paste0("GHS/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3_WGS84")),format = "GTiff",overwrite = TRUE)
    }
  }
}



#' Crop GHS Settlement Model at Threshold for IDF
#'
#'
IDF_smod <- function(cutoff = 20){
    l = list()
    pa = LandUseR:::paris_box()
    levs = LandUseR:::GHS_L2()
    cols = LandUseR:::GHS_LS_scale()

    for (yr in paste(c(1975, 1990, 2000, 2015))){
        smod <- raster::raster(file.path(LandUseR:::datadir(),paste0("GHS/GHS_SMOD_POP",yr,"_GLOBE_R2019A_54009_1K_V1_0_18_3/GHS_SMOD_POP",yr,"_GLOBE_R2019A_54009_1K_V1_0_18_3.tif")))

        pop <- raster::raster(file.path(datadir(),paste0("GHS/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3/GHS_POP_E",yr,"_GLOBE_R2019A_54009_250_V1_0_18_3.tif")))
        if(yr==2015){
          built <- raster::raster(file.path(datadir(),paste0("GHS/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3.tif")))
        } else {
          built <- raster::raster(file.path(LandUseR:::datadir(),paste0("GHS/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS",yr,"_GLOBE_R2018A_54009_250_V2_0_18_3.tif")))
        }


        pat = sf::st_transform(pa, raster::projection(smod))
        projected_smod = raster::projectRaster(smod, crs = "+proj=longlat +datum=WGS84 +no_defs")
        projected_built = raster::projectRaster(built, crs = "+proj=longlat +datum=WGS84 +no_defs")
        # pat = sf::st_transform(pa, raster::projection(pop))

        smod_paris = raster::crop(smod,raster::extent(pat))
        projected_smod_paris = raster::crop(projected_smod,raster::extent(pa))
        pop_paris = raster::crop(pop,raster::extent(smod_paris))
        built_paris = raster::crop(built,raster::extent(smod_paris))
        smod_paris <- raster::as.factor(smod_paris)   # declare as categorical data
        levs = subset(levs, ID %in% unique(values(smod_paris)))
        cols = subset(cols, ID %in% unique(values(smod_paris)))

        png(file.path(outdir(),paste0("plots/SMOD_paris",yr,".png")),width = 9,height=6,units = "in",res = 200)
        plot(smod_paris,legend = FALSE, col = cols$color,main=paste("Paris SMOD as is",yr))
        legend(x = 'bottomright', legend = levs$label, fill = cols$color)
        plot(projected_smod_paris,legend = FALSE, col = cols$color,main=paste("Paris SMOD as is",yr))
        legend(x = 'bottomright', legend = levs$label, fill = cols$color)
        dev.off()


        # aggregate  to 1km grids
        pop_paris = raster::aggregate(pop_paris, fact = 4, fun = sum)
        built_paris = raster::aggregate(built_paris, fact = 4, fun = mean)

        # get largest cluster to get main city
        s_mask = smod_paris
        raster::values(s_mask) = raster::values(smod_paris > cutoff)
        city <- smod_paris[s_mask,drop = FALSE]
        raster::plot(city)
        cc = clump(city)
        mask_city = cc
        values(mask_city) = values(cc == 1)
        values(mask_city)[values(mask_city) == 1] = TRUE
        values(mask_city)[values(mask_city) != 1] = NA
        totarea = cellStats(mask_city, sum)   # total area of main city

        # plot city model and mask
        png(file.path(outdir(),paste0("plots/SMOD_paris_classify",yr,".png")),width = 9,height=6,units = "in",res = 200)
        par(mfcol = c(1,2))
        plot(smod_paris,legend = FALSE, col = cols$color,main=paste("Paris SMOD as is",yr),xaxt = "n",yaxt = "n")
        legend(x = 'bottomright', legend = levs$label, fill = cols$color)
        plot(mask_city, main = "Classified as City",legend = FALSE, col = cols$color[length(cols$color)])
        dev.off()
        par(mfcol = c(1,1))


        # clip other data to that outline
        pop_paris = raster::crop(pop_paris,raster::extent(mask_city))
        built_paris = raster::crop(built_paris,raster::extent(mask_city))
        pop2 <- mask(pop_paris,mask_city)
        built2 <- mask(built_paris,mask_city)

        # plot both original and masked data
        png(file.path(outdir(),paste0("plots/maskpop-Paris",yr,".png")),width = 9,height=6,units = "in",res = 200)
        par(mfcol = c(1,2))
        plot(pop_paris,legend = FALSE,xaxt = "n",yaxt = "n")
        plot(pop_paris, legend.only=TRUE, horizontal = TRUE,
             legend.args = list(text='Population', side = 1, line = 2))
        plot(pop2, main = "Cut-out population",legend = FALSE,xaxt = "n",yaxt = "n")
        plot(pop2, legend.only=TRUE, horizontal = TRUE,
             legend.args = list(text='Population', side = 1, line = 2))
        dev.off()
        par(mfcol = c(1,1))

        png(file.path(outdir(),paste0("plots/maskbuilt-Paris",yr,".png")),width = 9,height=6,units = "in",res = 200)
        par(mfcol = c(1,2))
        plot(built_paris,legend = FALSE, main=paste("Paris built as is",yr),xaxt = "n",yaxt = "n")
        plot(built_paris, legend.only=TRUE, horizontal = TRUE,
             legend.args = list(text='Built Up', side = 1, line = 2))
        plot(built2, main = "Cut-out population",legend = FALSE,xaxt = "n",yaxt = "n")
        plot(built_paris, legend.only=TRUE, horizontal = TRUE,
             legend.args = list(text='Built Up', side = 1, line = 2))
        dev.off()
        par(mfcol = c(1,1))

        # get total built and total pop for masked areas
        mask_pop = pop_paris
        values(mask_pop)[!values(mask_city)] = 0
        plot(mask_pop)
        totpop = round(cellStats(mask_pop,sum))
        density = round(totpop / totarea)
        l[[yr]] <- list(smod = smod_paris, pop = pop_paris,main_city = mask_city, main_pop = mask_pop)

        png(file.path(outdir(),paste0("plots/SMOD_paris",yr,".png")),width = 9,height=6,units = "in",res = 200)
        plot(smod_paris, legend = FALSE, col = cols$color, main = paste0("Paris Settlement Model ",yr,". density: ",density), sub = paste0("population: ",totpop,". Area: ",totarea))
        legend(x = 'bottomright', legend = levs$label, fill = cols$color)
        dev.off()
    }

    # get
    return(l)
}



paris_box <- function(){
    if (!file.exists(file.path(datadir(),"paris-box.Rds"))){
      comms = sf::st_read(file.path(LandUseR:::datadir(),"communes-20190101/communes-20190101.json"))
      paris = sf::st_crop(comms, xmin = 1.61, ymin = 48.46, ymax = 49.118,xmax = 3.08)
      saveRDS(paris,file = file.path(datadir(),"paris-box.Rds"))
    } else {
      paris = readRDS(file.path(datadir(),"paris-box.Rds"))
    }
      # find bounding box for each city
    paris
}


plotIDF <- function(){
    # coastlines_sp <- rnaturalearth::ne_coastline(returnclass = "sp")
    #
    # worldmap <- rnaturalearth::ne_countries(scale = 'medium', type = 'map_units',returnclass = 'sf')
    # ggplot() + geom_sf(data = worldmap) + theme_bw()
    #
    # france <- worldmap[worldmap$name == 'France',]
    # ggplot() + geom_sf(data = france) + theme_bw()

    # read a raster file
    ra = list()
    ra_df = list()
    ra$built1975 = raster::raster(file.path(datadir(),"GHS/GHS_BUILT_LDS1975_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS1975_GLOBE_R2018A_54009_250_V2_0_18_3.tif"))
    ra$built2014 = raster::raster(file.path(datadir(),"GHS/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_18_3.tif"))

    ra$pop1975 = raster::raster(file.path(datadir(),"GHS/GHS_POP_E1975_GLOBE_R2019A_54009_250_V1_0_18_3/GHS_POP_E1975_GLOBE_R2019A_54009_250_V1_0_18_3.tif"))
    ra$pop2015 = raster::raster(file.path(datadir(),"GHS/GHS_POP_E2015_GLOBE_R2019A_54009_250_V1_0_18_3/GHS_POP_E2015_GLOBE_R2019A_54009_250_V1_0_18_3.tif"))


    # does not work: 16million rows
    # ggplot() + geom_raster(data = ra_df, aes(x,y,alpha = people))

    # cropping france
    # france_cropped <- st_crop(worldmap, xmin = 0, xmax = 2,
    #                           ymin = 47, ymax = 48)
    # ggplot() + geom_sf(data = france_cropped) + theme_bw()


    # plot pop and built up count for paris
    comms = st_read(file.path(datadir(),"communes-20190101/communes-20190101.json"))
    paris = st_crop(comms, xmin = 1.7, ymin = 48.46, ymax = 49.118,xmax = 3.08)  # find bounding box for each city

    # transform to raster's crs
    paris = st_transform(paris, raster::projection(ra$built1975))
    paris_pop = st_transform(paris, raster::projection(ra$pop1975))

    # compute total built up area by usign the clump function
    # set a threshold for what we consider "city"
    thresh = 50

    # crop raster to current city
    rparis = raster::crop(ra[[1]],raster::extent(paris))
    rparis_pop = raster::crop(ra[[2]],raster::extent(paris))
    png(file.path(outdir(),"plots/paris1975.png"),width = 9,height=6,units = "in",res = 200)
    plot(rparis,main = "Paris fraction built up")
    dev.off()
    png(file.path(outdir(),"plots/paris1975_pop.png"),width = 9,height=6,units = "in",res = 200)
    plot(rparis_pop,main = "Paris cell population")
    dev.off()

    # impose threshold
    I = rparis < thresh
    rparis[I] = 0
    png(file.path(outdir(),"plots/paris1975_thresh.png"),width = 9,height=6,units = "in",res = 200)
    plot(rparis,main = paste("Paris fraction built up (imposed threshold of min",thresh,"%)"))
    dev.off()

    # compute


    # crop all rasters to IDF
    for (g in 1:length(ra)){
        ra[[g]] = raster::crop(ra[[g]],raster::extent(comms))
        #
        ra_df[[g]] = raster::as.data.frame(ra[[g]],xy =TRUE)
        ra_df[[g]] = raster::as.data.frame(ra[[g]],xy =TRUE)
    }
    names(ra_df) <-names(ra)
    names(ra_df$built1975)[3] <- "built_up_density"
    names(ra_df$built2014)[3] <- "built_up_density"
    names(ra_df$pop1975)[3] <- "pop_density"
    names(ra_df$pop2015)[3] <- "pop_density"

    # plots
    p = list()
    p$built1975 = ggplot() + geom_sf(data = comms) + theme_bw() + geom_raster(data = ra_df$built1975, aes(x,y,alpha = built_up_density))

    plist = list(built = ra_df[1:2], pop = ra_df[3:4])
    names(plist$built) <- c("1975","2014")
    names(plist$pop) <- c("1975","2015")

    ggplot() + geom_sf(data = comms) + theme_bw() + geom_raster(data = ra_df$`z`, aes(x=x,y=y,alpha = pop_density))
    pop_plots = lapply(names(plist$pop), function(z) {
        ggplot() + geom_sf(data = comms, fill = NA) + theme_bw() + geom_raster(data = plist$pop[[z]], aes(x=x,y=y,alpha = pop_density)) + ggtitle(paste("IDF Population Density in",z))
    })
    names(pop_plots) <- names(plist$pop)

    built_plots = lapply(names(plist$built), function(z) {
        ggplot() + geom_sf(data = comms, fill = NA) + theme_bw() + geom_raster(data = plist$built[[z]], aes(x=x,y=y,alpha = built_up_density)) + ggtitle(paste("IDF Built Up Density in",z))
    })
    names(built_plots) <- names(plist$built)
    ggsave(file.path(outdir(),"plots/pop1975.png"),plot = pop_plots$`1975`)
    ggsave(file.path(outdir(),"plots/pop2015.png"),plot = pop_plots$`2015`)
    ggsave(file.path(outdir(),"plots/built1975.png"),plot = built_plots$`1975`)
    ggsave(file.path(outdir(),"plots/built2014.png"),plot = built_plots$`2014`)

}






# archive
#
# europe <- worldmap[worldmap$continent == 'Europe',]
# ggplot() + geom_sf(data = europe) + theme_bw()
#
# # cropping
# europe_cropped <- st_crop(worldmap, xmin = -20, xmax = 45,
#                           ymin = 30, ymax = 73)
# ggplot() + geom_sf(data = europe_cropped) + theme_bw()
#
# # simpliyfing polygons
# coastlines_sim2 <- rgeos::gSimplify(coastlines_sp,
#                              tol = .1,
#                              topologyPreserve = TRUE)
# coastlines_sim2_df <- st_as_sf(coastlines_sim2)
# ggplot(coastlines_sim2_df) + geom_sf()

