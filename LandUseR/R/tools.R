

setup_gsheet <- function(){
    googledrive::drive_upload(file.path(datadir(),"france_cities.csv"), type = "spreadsheet")
}
dboxdir <- function(){file.path(Sys.getenv("HOME"),"Dropbox","research","LandUse")}
datadir <- function(){file.path(dboxdir(),"data")}
dataplots <- function(){file.path(outdir(),"data","plots")}

entddir <- function(){file.path(datadir(),"ENTD")}

outdir <- function(){file.path(dboxdir(),"output")}
plotdir <- function(){file.path(outdir(),"plots")}

gsheet <- function(){"https://docs.google.com/spreadsheets/d/1IotLzprM5Os-06MpfU895VOgHieYkQfY3WOLBBiTKZA/edit#gid=0"}

export_sheet <- function(topn = 110){
    p = readpop()
    p = p %>% dplyr::filter(date == "1876-01-01") %>% dplyr::arrange(desc(population))
    # p = p %>% dplyr::filter(year == 2015) %>% dplyr::arrange(desc(population))
    paris = p %>%
        dplyr::filter(CODGEO > "75100" & CODGEO < "751210") %>%
        dplyr::summarise(CODGEO = "75060", REG = REG[1], DEP = DEP[1], LIBGEO = "Paris", year = year[1], population = sum(population), date = date[1])
    p0 = p %>% dplyr::filter(date == "1876-01-01",(CODGEO < "75101" | CODGEO > "75120")) %>% dplyr::arrange(desc(population)) %>% dplyr::top_n(topn,population)
    # p0 = p %>% dplyr::filter(year == 2015,(CODGEO < "75101" | CODGEO > "75120")) %>% dplyr::arrange(desc(population)) %>% dplyr::top_n(topn,population)
    r = dplyr::bind_rows(p0,paris) %>% dplyr::arrange(desc(population)) %>% dplyr::mutate(rank = 1:nrow(.), area_EM = 0.0, area_1950 = 0.0,area_2016 = 0.0, pop_1876 = population, pop_1950 = 0)
    r = r %>%
        dplyr::select(rank,CODGEO,REG,DEP,LIBGEO,pop_1876, pop_1950, area_EM, area_1950,area_2016)
        # dplyr::filter(!(CODGEO %in% c("59512","78646","59599","93066","92044","92012","92051")))

    # ss = googlesheets4::read_sheet(LandUseR:::gsheet())
    googlesheets4::sheets_write(r,ss = LandUseR:::gsheet(),sheet = "manual-collection-upd")
    # googlesheets4::sheets_read(LandUseR:::gsheet(),sheet = "manual-collection")
}

import_sheet <- function(sheet = NULL){
    googlesheets4::read_sheet(LandUseR:::gsheet(),sheet = sheet)
}

get_manuals <- function(reload = FALSE){
    if (reload){
        s = LandUseR:::import_sheet(sheet = "master-list") %>% dplyr::filter(area_2016 != 0.0 & is.na(which_city_2016))
        sd = data.table(s)
        sd[,which_city_2016 := NULL]
        saveRDS(sd,file = file.path(datadir(),"top100.Rds"))
    } else {
        sd = readRDS(file = file.path(datadir(),"top100.Rds"))
    }
    sd

}


writepop <- function(){
    pop = readxl::read_excel(file.path(datadir(),"base-pop-historiques-1876-2015.xls"),skip = 5)
    names(pop)[5:33] <- paste0(c(2015,2014	,2013,2012	,2011	,2010	,2009	,2008	,2007	,2006	,1999	,1990	,1982	,1975	,1968	,1962	,1954	,1936	,1931	,1926	,1921	,1911	,1906	,1901	,1896	,1891	,1886	,1881	,1876), "-01-01")
    pop %<>%
        tidyr::gather(year, population , -c(CODGEO,REG,DEP,LIBGEO)) %>%
        dplyr::mutate(date = as.Date(year), DEP = as.integer(DEP), REG = as.integer(REG), year = lubridate::year(date))
    saveRDS(pop,file.path(datadir(),"base-pop-historiques-1876-2015.Rds"))
    return(pop)
}

#' Read Census Data
#'
#' @export
readpop <- function(){
    readRDS(file.path(datadir(),"base-pop-historiques-1876-2015.Rds"))
}


plot_reims <- function(){
    p = readpop()
    y=p %>% dplyr::filter(CODGEO=="51454", year %in% c(1876,2015))
    y %<>% dplyr::mutate(area = c(50.52,3.29),density = population / area)
    d1 = data.table(variable = c("population","area","density"), increase = unlist(c(y[1,"population"] / y[2,"population"], y[1,"area"] / y[2,"area"], NA)))


    d2 = data.table(variable = factor(c("population","area","density"), levels = c("population","area","density")), increase = unlist(c(y[1,"population"] / y[2,"population"], y[1,"area"] / y[2,"area"], -y[2,"density"] / y[1,"density"])))

    cols <- c("area" = "red", "population" = "green", "density" = "grey")

    p1 = ggplot(d1,aes(x=variable,increase)) + geom_col(aes(fill=variable),alpha=0.6,color = "black") + theme_bw() + scale_fill_manual(values = cols) + scale_y_continuous(name = "Increased by Factor",limits = c(-10,16),breaks = c(-7,0,3,10,15),minor_breaks = NULL) + geom_hline(yintercept = 0,size = 1) + ggtitle("Reims from 1866 to 2015")  + scale_x_discrete(name = "") + theme(legend.position = "none")

    p2 = ggplot(d2,aes(x=variable,increase)) + geom_col(aes(fill=variable),alpha=0.6,color = "black") + theme_bw() + scale_fill_manual(values = cols) + scale_y_continuous(name = "Increased by Factor",limits = c(-10,16),breaks = c(-7,0,3,10,15),minor_breaks = NULL) + geom_hline(yintercept = 0,size = 1) + ggtitle("Reims from 1866 to 2015") + scale_x_discrete(name = "")+ theme(legend.position = "none")

    ggsave(plot = p1, filename=file.path(dataplots(),"reims1.pdf"),width=7,height=5)
    ggsave(plot = p2, filename=file.path(dataplots(),"reims2.pdf"),width=7,height=5)

   return(list(p1,p2))
}
