
# file path setup


setup_gsheet <- function(){
    googledrive::drive_upload(file.path(datadir(),"france_cities.csv"), type = "spreadsheet")
}
dboxdir <- function(){file.path(Sys.getenv("HOME"),"Dropbox","research","LandUse")}
datadir <- function(){file.path(dboxdir(),"data")}
dataplots <- function(){file.path(outdir(),"data","plots")}
datatables <- function(){file.path(outdir(),"data","tables")}

entddir <- function(){file.path(datadir(),"ENTD")}

outdir <- function(){file.path(dboxdir(),"output")}
outdatadir <- function(){file.path(dboxdir(),"output","data")}
instoutdir <- function(){system.file("data",package = "LandUseR")}
instindir <- function(){file.path(here::here(),"inst","data")}
