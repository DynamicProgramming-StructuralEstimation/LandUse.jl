# resource: https://geocompr.robinlovelace.net
# resource: https://datacarpentry.org/r-raster-vector-geospatial/

# This uses data produced by the European Commission for their Global Human Settlement
# database GHS.
# https://ghslsys.jrc.ec.europa.eu/ghs_smod2019.php

GHS_years <- function() c(1975,1990,2000,2015)


#' Crop boxes out of GHS rasters
#'
cropboxes <- function(overwrite = FALSE){
  if (overwrite){
    x = LandUseR:::pop_1950_2()
    lgs = x[,unique(CODGEO)]
    L = list()
    for (iy in 1:length(GHS_years())){

      L[[iy]] <- list()

      yr = GHS_years()[iy]
      flog.info("cropping year %d",yr)
      r = loadRasters(yr = yr)

      pb <- progress::progress_bar$new(total = 100)
      for (lg in lgs){

        L[[iy]][[lg]] = crop_rasters(r,x[CODGEO==lg,extent][[1]])
        pb$tick()
      }

    }
    names(L) <- paste(GHS_years())
    saveRDS(L,file.path(outdir(),"data","france-cropped-boxes.Rds"))
  } else {
    L = readRDS(file.path(outdir(),"data","france-cropped-boxes.Rds"))
  }
  L
}

#' Measure City Extents
#'
#' Takes cropped cities and counts the number of grid cells which are classified as "urban area".
#'
#' @param overwrite TRUE/FALSE
#' @param cutoff numeric indicating the cutoff value (percent) of built-up area above which we classify as "urban"
#'
#' We read the list of cropped city bounding boxes from \code{\link{cropboxes}} in order to classify them as either urban area or not.
#' The classification of urban or not is governed by two constraints:
#' \enumerate{
#' \item the `cutoff` parameter. If a grid cell is measured with built-up surface of more than `cutoff`, the cell is retained as `urban`.
#' \item Contiguity of grid cells. Given that classified grid cells could be separated in space, we check for the largest connect set. The check uses the "rook move" as criterion.
#' }
#'  As a second
measure_cities <- function(overwrite = FALSE, cutoff = 30){
    if (overwrite){
        # get list of cities
        L = readRDS(file.path(LandUseR:::outdir(),"data","france-cropped-boxes.Rds"))

        # get centers
        ct = merge_centers(overwrite = FALSE)

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
                imask <- citymask
                raster::values(imask) <- is.na(raster::values(citymask))
                raster::values(imask)[!raster::values(imask)] <- NA

                # store final cut out areas for populaiton and built up back on list
                L[[iy]][[ic]]$built_cut = cityb[citymask, drop = FALSE]
                L[[iy]][[ic]]$pop_cut = cityp[citymask, drop = FALSE]

                # store inverse city mask - everything that is *not* city
                L[[iy]][[ic]]$inverse_mask <- imask

                # store city center on list
                L[[iy]][[ic]]$center = ct[CODGEO == ic, c(center_x, center_y)]
                L[[iy]][[ic]]$cityname = ct[CODGEO == ic, LIBGEO]

                totarea = raster::cellStats(citymask, "sum",na.rm = T)   # sum of gridcells
                totarea = totarea * (250 * 250) / 1e+6  # total area in square kilometers

                totarea_rural = raster::cellStats(imask, "sum",na.rm = T)   # sum of gridcells
                totarea_rural = totarea_rural * (250 * 250) / 1e+6  # total area in square kilometers

                totpop = sum(cityp[citymask],na.rm = T)
                meanbuilt = mean(cityb[citymask],na.rm = T)
                qbuilt = quantile(cityb[citymask],na.rm = T)
                qpop = quantile(cityp[citymask],na.rm = T)

                # output
                OL[[iy]][[ic]] = list(area = totarea, rural_area = totarea_rural, pop = totpop, meanbuilt = meanbuilt,
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
        saveRDS( L,file.path(LandUseR:::outdir(),"data","france-cropped-cities.Rds"))

    } else {
        OL = readRDS(file.path(LandUseR:::outdir(),"data","france_measured.Rds"))
        L = readRDS(file.path(LandUseR:::outdir(),"data","france-cropped-cities.Rds"))
    }
    return(list(measured = OL,cropped = L))
}


#' Compute Distance of Fringe from City Center
#'
#' Read cut-out city rasters from \code{\link{measure_cities}} and compute
#' distance of each point in raster to the point called \emph{city center}.
#' City Center is classified by IGN as \emph{Chef Lieu} of a place.
#'
#' @param overwrite TRUE/FALSE whether to load from disk or recreate
#' @param ndists numeric value indicating how many distance bins
#'
#' @return returns a data.table with distance, density by city and year.
dist_from_center <- function(overwrite = FALSE, ndists = 51){

    if (overwrite) {
        # read cut out city list
        L = measure_cities()$cropped
        OL = list()

        d0 = data.table()
        # for each year and each city
        for (iy in 1:length(L)){
            yr = GHS_years()[iy]

            OL[[iy]] <- list()
            flog.info("measuring distances in year %d",yr)

            for (ic in names((L[[iy]]))){

                dists = raster::distanceFromPoints(L[[iy]][[ic]]$pop_cut, L[[iy]][[ic]]$center)
                dists = dists[!is.na(L[[iy]][[ic]]$pop_cut),drop = FALSE]  # keep only cells that are not NA
                md = raster::cellStats(dists, "max",na.rm = T)
                dt = seq(0,md,length.out = ndists)
                rc = cut(raster::values(dists), breaks = dt, label = FALSE)
                # now compute density at each quantile
                d = data.table(quantile = rc, d = raster::values(L[[iy]][[ic]]$pop_cut) )
                d <- d[,list(density = mean(d,na.rm=T), pop = sum(d,na.rm=T)),by = quantile]
                d[, c("distance", "CODGEO", "LIBGEO","year") := list(dt[quantile], ic, L[[iy]][[ic]]$cityname, yr)]
                d0 = rbind(d0,d)
            }
        }
        d0 = d0[complete.cases(d0)]
        saveRDS(d0,file.path(LandUseR:::outdir(),"data","density-distance.Rds"))
        fwrite(d0,file.path(LandUseR:::outdir(),"data","density-distance.csv"))
        return(d0)
    } else {
        readRDS(file.path(LandUseR:::outdir(),"data","density-distance.Rds"))
    }
}


#' Combine City Measures
#'
#' combine satellite measures from \code{\link{measure_cities}} and manual measures
#' from \code{\link{pop_1950_2}} to create the final output
#' dataset.
#'
#' @param cutoff numeric value above which classify as urban
#' @param overwrite TRUE/FALSE
combine_measures <- function(overwrite = FALSE,cutoff = 30){
    # melt manual measures
    x = LandUseR:::pop_1950_2(overwrite = overwrite)
    setnames(x,"area_EM","area_1876")

    # overwrite rank column because that is non-contiguous
    x[, rank := order(pop_1876, decreasing = TRUE)]
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
        ss = sat$measured[[paste(yr)]]
        q = rbindlist(ss)
        q[, c("year","CODGEO") := list(yr,names(ss))]
        ll[[paste0(yr)]] <- merge(m2[.(yr)], q, by = c("year","CODGEO"))
    }
    m2 = rbindlist(ll)
    m2[,type := "satellite"]
    z = rbind(m,m2,fill = TRUE)
    setkey(z, year, rank)
    saveRDS(z, file.path(LandUseR:::outdir(),"data",paste0("france_final_cutoff",cutoff,".Rds")) )
    saveRDS(z, file.path(instindir(),paste0("france_final_cutoff",cutoff,".Rds")) )  # write to package internatl store
    fwrite(z[,!"extent"],  file.path(LandUseR:::outdir(),"data",paste0("france_final_cutoff",cutoff,".csv")))
    z
}


#' Cutoff Sensitivity Analysis for `measure_cities`
#'
#' @param cuts desired cutoff levels
#'
#' compute output measures for different values of `cutoff`
cutoff_sensitivity <- function(cuts = c(25,30,35,50),overwrite = FALSE){

    if (overwrite) {
        L = list()
        for (cu in cuts) {
            futile.logger::flog.info("Computing City measures for cutoff %d",cu)
            L[[paste(cu)]] = combine_measures(overwrite = TRUE, cutoff = cu)
        }
    } else {
        L = list()
        for (cu in cuts) {
            futile.logger::flog.info("Reading cutoff %d from disk",cu)
            L[[paste(cu)]] = readRDS(file.path(LandUseR:::outdir(),"data",paste0("france_final_cutoff",cu,".Rds")) )
        }
    }
    L
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


#' get CRS of GHS raster data
#'
get_raster_CRS <- function(){
  fi = file.path(LandUseR:::datadir(),paste0("GHS/GHS_POP_E",1975,"_GLOBE_R2019A_54009_250_V1_0_17_3/GHS_POP_E",1975,"_GLOBE_R2019A_54009_250_V1_0_17_3",".tif"))
  if (!file.exists(fi)) {
    stop(sprintf("file %s does not exist. dropbox synced?",fi))
  } else {
    pop <- raster::raster(fi)
    raster::crs(pop)
  }
}


#' Load GHS Raster Data
#'
#' @param yr valid GHS year to load: 1975,1990,2000,2015
loadRasters <- function(yr){
  stopifnot(yr %in% c(1975,1990,2000,2015))
  li = list()
  ending = ".tif"

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
