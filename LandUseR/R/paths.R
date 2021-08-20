
# file path setup

#' Get Dropbox Root
#'
#' This returns the root of our shared dropbox folder on your machine.
#' On my computer this would be `~/Dropbox/research/LandUse`
#' @export
dboxdir <- function() {
    db = Sys.getenv("R_LANDUSE")
    if (db == ""){
        stop("You need to set ENV var R_LANDUSE to the entry point\n of your `LandUse` dropbox. for me that is\n
        R_LANDUSE=~/Dropbox/research/LandUse .\n
        you set this by editing file ~/.Renviron .\n mine has this line for this purpose: \n R_DROPBOX=~/Dropbox/research/LandUse\n
        the easiest option is to just call this function `usethis::edit_r_environ()` from the `usethis` package")
    } else {
        return(db)
    }
}

datadir <- function(){file.path(dboxdir(),"data")}
dataplots <- function(){file.path(outdir(),"data","plots")}
datatables <- function(){file.path(outdir(),"data","tables")}
entddir <- function(){file.path(datadir(),"ENTD")}
outdir <- function(){file.path(dboxdir(),"output")}
outdatadir <- function(){file.path(dboxdir(),"output","data")}
instoutdir <- function(){system.file("data",package = "LandUseR")}
instindir <- function(){file.path(here::here(),"inst","data")}

setup_gsheet <- function(){
    googledrive::drive_upload(file.path(datadir(),"france_cities.csv"), type = "spreadsheet")
}


