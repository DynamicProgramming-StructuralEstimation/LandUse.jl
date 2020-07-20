# enquete national du transport et des deplacements 2008
# https://www.statistiques.developpement-durable.gouv.fr/enquete-nationale-transports-et-deplacements-entd-2008

# questionnaire: https://www.statistiques.developpement-durable.gouv.fr/sites/default/files/2018-12/ENTD_2007-2008_Questionnaire.pdf


#' Read ENTD data from website or disk
#'
readENTD <- function(reload = FALSE){
    if (reload){
        l = list()
        l$famdemo <- fread("https://www.statistiques.developpement-durable.gouv.fr/sites/default/files/2018-12/Q_tcm_menage_0.csv")
        l$inddemo <- fread("https://www.statistiques.developpement-durable.gouv.fr/sites/default/files/2019-01/Q_tcm_individu.csv")
        l$work <- fread("https://www.statistiques.developpement-durable.gouv.fr/sites/default/files/2019-01/Q_ind_lieu_teg.csv")
        l$work[,Map := lubridate::hms(Map)]
        l$work[,Metro := lubridate::hms(Metro)]
        l$work[,RER := lubridate::hms(RER)]
        l$work[,TC := lubridate::hms(TC)]
        saveRDS(l, file.path(entddir(),"ENTD2008.Rds"))
    } else {
        l = readRDS(file.path(entddir(),"ENTD2008.Rds"))
    }
    return(l)
}

