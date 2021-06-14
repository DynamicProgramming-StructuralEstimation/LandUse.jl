

# corine land cover https://land.copernicus.eu/pan-european/corine-land-cover

#' Crop CLC 2018 raster to 2015 French Cities
#'
#' load CLC raster and cut out our bounding boxes of French cities
CLC_bboxes <- function(overwrite = FALSE){
    if (overwrite) {
        ms = measure_cities()$cropped$`2015`
        if (!file.exists(file.path(LandUseR:::datadir(),"france_CLC.tif"))) {
            futile.logger::flog.info("reading CLC raster")
            CLC = raster::raster(file.path(LandUseR:::datadir(),"u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"))
            # france
            CLCF = raster::crop(CLC,raster::extent(3100000, 4300000, 2100000, 3200000))
            # CLCF <- raster::projectRaster(CLCF, crs = raster::crs(ms[[1]]$inverse_mask), method = "ngb")  # does not work! it interpolates data values and fucks up everything
            raster::writeRaster(CLCF, file.path(LandUseR:::datadir(),"france_CLC.tif"))
        } else {
            CLCF = raster::raster( file.path(LandUseR:::datadir(),"france_CLC.tif"))
        }

        OL = list()
        futile.logger::flog.info("reading CLC raster")

        for (ix in 1:length(ms)){
            # bounding box in CLC crs
            bb = sf::st_bbox(ms[[ix]]$inverse_mask)  # sf to the rescue
            bb_ll = sf::st_bbox(
                sf::st_transform(
                    sf::st_as_sfc(bb), crs = raster::projection(CLCF))
            )
            e = raster::extent(c(bb_ll["xmin"],bb_ll["xmax"],bb_ll["ymin"],bb_ll["ymax"]))
            OL[[ix]] <- list()
            OL[[ix]]$full <- raster::crop(CLCF, e)

            # inverse_mask in CLC crs
            ivm = raster::projectRaster(ms[[ix]]$inverse_mask, to = OL[[ix]]$full, method = "ngb")
            OL[[ix]]$cut  <- OL[[ix]]$full[ivm, drop = FALSE]
            OL[[ix]]$cityname  <- ms[[ix]]$cityname
        }
        names(OL) <- names(ms)
        saveRDS(OL,file.path(LandUseR:::outdir(),"data","CLC-France-cropped.Rds"))
        OL
    } else {
        readRDS(file.path(LandUseR:::outdir(),"data","CLC-France-cropped.Rds"))
    }
}

#' Read CLC classification
#'
CLC_read_legend <- function(){
    x = fread(file.path(LandUseR:::datadir(),"u2018_clc2018_v2020_20u1_raster100m/Legend/CLC2018_CLC2018_V2018_20_QGIS.txt"))
    setnames(x, c("code","red","green","blue","alpha","label"))
    x[,color := rgb(red,green,blue, maxColorValue = alpha)]
    x[,mcode := factor(c(1:44,48))]  # get out of u2018_clc2018_v2020_20u1_raster100m/Legend/clc_legend_qgis_raster.qml
    x[mcode == 21, label := "Agr. Land w/ sign. areas of natural vegetation"]
    x[,legend := factor(code,labels = label)]
    x
}



#' Measure Landuse Around Cities with CLC
#'
#' takes cropped CLC rasters for top 100 cities and counts the proportion of land
#' outside city falling into each category
#'
CLC_measure <- function(){
    le = CLC_read_legend()
    cuts = CLC_bboxes(FALSE)

    r = lapply(cuts, function(z){ r <- raster::values(z$cut);
                                  r <- r[!is.na(r)];
                                  fr = factor(r, levels = levels(le[,mcode]))
                                  prop.table(table(fr))
        })
    df = data.table(mcode = names(r[[1]]), avg_prop = rowMeans(matrix(unlist(r),ncol = length(r), byrow = F)))
    df2 = merge(df,le[,list(mcode,legend,color)])
    df2 = df2[order(avg_prop,decreasing = TRUE)]

    colnames = df2[, as.character(color)]
    names(colnames) <- df2[, as.character(legend)]
    gg = ggplot(df2[1:15], aes(x = reorder(legend,avg_prop) ,y = avg_prop)) + geom_bar(stat = "identity", aes(fill = legend)) + coord_flip() + ggtitle("Average Land Use Outside top 100 Cities")  + scale_fill_manual(values = colnames) + theme_bw()  + theme(legend.position = "none") + scale_x_discrete(name = "") + scale_y_continuous("Proportion")

    ggsave(gg, file = file.path(dataplots(),"CLC-landuse-top100.pdf"), width = 8 , height = 4.5)
    gg
}


#' Plot CLC Landuse patterns
#'
CLC_plots <- function(){
    le = CLC_read_legend()
    cuts = CLC_bboxes(FALSE)

    OL = list()
    for (insee in 1:4) {
        r <- cuts[[insee]]$cut
        v = raster::values(r)
        raster::values(r) <- raster::as.factor(v)
        rat <- raster::levels(r)[[1]]
        rat = merge.data.table(rat, le[,list(mcode,legend,color)],by.x = "VALUE", by.y = "mcode")
        setDT(rat)
        setcolorder(rat, "ID")
        setkey(rat,ID)
        cols = rat[,color]
        rat[, c("VALUE","color") := NULL]

        levels(r) <- rat

        OL[[insee]] <- rasterVis::levelplot(r,col.regions = cols,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
        lattice::trellis.device(pdf, file=file.path(dataplots(),paste0("CLC-",cuts[[insee]]$cityname,".pdf")),height=4.5, width=8)
        print(OL[[insee]])
        dev.off()
    }
    OL
}


