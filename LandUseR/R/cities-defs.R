
# This file has functionality to define city populations manually.




#' Compute French 1954 City Population
#'
#' This function aggregates `communes` into greater urban areas according to our own definition.
#' We establish extent of a city by examining aereal photography and maps from 1950.
#' Larger cities may contain separate communes as defined by INSEE (and hence, with a separate entry in the census counts).
pop_1950_2 <- function(overwrite = FALSE){
    if (overwrite){
        x = get_manuals()   # get manual area measurements for top 100 cities in france
        r = merge_centers(overwrite = FALSE) # get extent and centers of each city
        x=merge(x,r,by=c("LIBGEO","CODGEO")) #Merge by name of city
        p0 = readpop() # get census population counts 1860-2016
        p = as.data.table(p0)  # get census population counts 1860-2016
        p54 = p[year == 1954]   # get the fifties count

        # add paris
        paris = p54[DEP == 75,list(CODGEO = "75060",
                                   REG = 11, DEP = 75,
                                   LIBGEO = "Paris", year = 1954,
                                   population = sum(LandUseR:::get_paris_pop_1950(p)$population),date = date[1])]
        p54 <- p54[DEP != 75]
        p54 <- rbind(p54,paris)
        # paris is done.

        # the `population` field is incomplete for some remaining cities. it only contains the main city.
        setkey(p54,CODGEO)
        setkey(x,CODGEO)
        xp = merge(x,p54[,.(CODGEO,population)])
        setkey(xp,rank)
        setkey(p54,DEP)

        # complete city (cc) table: which CODGEO or LIBGEO are part of a 1950 urban area?
        cc <- xp[,.(CODGEO,DEP,LIBGEO,population)]
        cc[, components := vector("list")]
        setkey(cc,DEP,LIBGEO)  # easier to read than CODGEO
        # cc[.(75,"Paris") , components ]  # thats done. for cities with only one component, do nothing (no need to put NULL or anything)

        # https://www.toulouse-metropole.fr/collectivite/communes
        # https://www.plan.toulouse.fr/map/?t=TOULOUSE_METROPOLE
        # https://fr.wikipedia.org/wiki/Nice#Localisation
        cc[.(69,"Lyon") , components                 := c("Lyon","Villeurbanne","Caluire-et-Cuire", "Oullins")]  # for others, insert components
        cc[.(44,"Nantes") , components               := c("Nantes", "Rezé")]
        cc[.(67,"Strasbourg") , components           := c("Strasbourg","Schiltigheim","Bischheim", "Hoenheim")]
        cc[.(33,"Bordeaux") , components             := c("Bordeaux","Talence","Bègles","Le Bouscat")]  #As with Tours, Caudéran was a commune until 1965
        cc[.(59,"Lille") , components                := c("Lille","La Madeleine")]
        cc[.(51,"Reims") , components                := c("Reims", "Cormontreuil")] #As with Tours the town of Neuvillette was included in Reims since 1970 but is not in the list
        cc[.(76,"Le Havre") , components             := c("Le Havre", "Sainte-Adresse")]
        cc[.(83,"Toulon") , components               := c("Toulon", "La Valette-du-Var")]
        cc[.(63,"Clermont-Ferrand") , components     := c("Clermont-Ferrand", "Chamalières")]
        cc[.(57,"Metz") , components                 := c("Metz", "Montigny-lès-Metz", "Longeville-lès-Metz")]
        cc[.(43,"Troyes") , components               := c("Troyes", "Sainte-Savine","Les Noës-près-Troyes","Saint-André-les-Vergers","La Rivière-de-Corps","La Chapelle-Saint-Luc","Saint-Julien-les-Villas")]
        cc[.(62,"Boulogne-sur-Mer") , components     := c("Boulogne-sur-Mer", "Saint-Martin-Boulogne","Outreau","Le Portel")]
        cc[.(2,"Saint-Quentin") , components         := c("Saint-Quentin", "Harly","Gauchy")]
        cc[.(62,"Calais") , components               := c("Calais", "Coulogne")]
        cc[.(59,"Douai") , components                := c("Douai", "Dechy")]
        cc[.(62,"Arras") , components                := c("Arras", "Achicourt")]
        cc[.(59,"Valenciennes") , components         := c("Valenciennes", "Marly","Saint-Saulve","La Sentinelle","Anzin","Trith-Saint-Léger","Beuvrages","Raismes","Bruay-sur-l'Escaut","Petite-Forêt","Aulnoy-lez-Valenciennes")]
        cc[.(42,"Saint-Chamond") , components        := c("Saint-Chamond", "L'Horme")]
        cc[.(3,"Montluçon") , components             := c("Montluçon", "Désertines")]
        cc[.(76,"Elbeuf") , components               := c("Elbeuf", "Caudebec-lès-Elbeuf","Saint-Aubin-lès-Elbeuf")]
        cc[.(59,"Cambrai") , components              := c("Cambrai", "Proville","Neuville-Saint-Rémy")]
        cc[.(3,"Moulins") , components               := c("Moulins", "Yzeure")]
        cc[.(28,"Chartres") , components             := c("Chartres", "Mainvilliers", "Luisant")]
        cc[.(51,"Châlons-en-Champagne") , components := c("Châlons-en-Champagne", "Saint-Memmie")]
        cc[.(43,"Le Puy-en-Velay") , components      := c("Le Puy-en-Velay", "Vals-près-le-Puy")]
        cc[.(14,"Lisieux") , components              := c("Lisieux", "Saint-Désir")]
        cc[.(36,"Châteauroux") , components          := c("Châteauroux", "Déols")]
        cc[.(29,"Morlaix") , components              := c("Morlaix", "Saint-Martin-des-Champs")]


        # incomplete cities
        ic = cc[!unlist(lapply(components,is.null)), ]
        setkey(ic,DEP)
        # sum over all populations of `components` by city in ic
        # and plug into xp via CODGEO of main city
        for (j in 1:nrow(ic)){
            irow = ic[j]  # grab the j-th row of ic
            xp[ .(irow[,.(CODGEO,DEP)])   , pop_1950 := p54[.(irow[,DEP])][LIBGEO %in% irow[, components[[1]]], sum(population)] , on = c("CODGEO","DEP")]
        }

        # for complete cities, we can just take the `population column`
        ic = cc[unlist(lapply(components,is.null)), CODGEO]
        stopifnot(xp[CODGEO %in% ic, all(pop_1950 == 0)])
        xp[CODGEO %in% ic, pop_1950 := population]

        #Exceptions: Mâcon and Vienne have inglobated towns which are in different departements
        #Mâcon
        pop1_M=xp[LIBGEO=="Mâcon", pop_1950]
        pop2_M=p54[LIBGEO=="Saint-Laurent-sur-Saône",population]
        xp[LIBGEO=="Mâcon", pop_1950:=(pop1_M+pop2_M)]

        #Vienne
        pop1_V=xp[LIBGEO=="Vienne", pop_1950]
        pop2_V=p54[LIBGEO=="Sainte-Colombe" & DEP==69 ,population]
        pop3_V=p54[LIBGEO=="Saint-Romain-en-Gal",population]
        xp[LIBGEO=="Vienne", pop_1950:=(pop1_V+pop2_V+pop3_V)]

        xp[,population := NULL]

        saveRDS(xp, file.path(datadir(),"top100bbox.Rds"))
        saveRDS(cc, file.path(datadir(),"top100bbox-components.Rds"))

    } else {
        xp = readRDS(file.path(datadir(),"top100bbox.Rds"))
    }
    xp
}

#' Print France 1950 Population to Latex
#'
#'
print_pop_1950 <- function(){
    x = readRDS(file.path(datadir(),"top100bbox-components.Rds"))
    y = x[!(unlist(lapply(components,is.null))), !"population"]
    ta = knitr::kable(y, format = "latex",
                      booktabs = TRUE,
                      longtable = FALSE,
                      caption = "France 1950 Population Classification. Cities containing more than one INSEE administrative area by 1950.",
                      label = "france-1950-pop"
    ) %>% kableExtra::column_spec(4, width = "10cm")
    # cat(kableExtra::landscape(ta), file = file.path(LandUseR:::datatables(), "france-pop-1950.tex"))
    cat(ta, file = file.path(LandUseR:::datatables(), "france-pop-1950.tex"))
}



#' sum city population over communes
#'
#' assumes that list of communes li has the main city as first entry!
sum_city <- function(dt,li){
    tmp = dt[LIBGEO == li[1]]
    tmp[,population := sum(dt[,population])]
    return(tmp)
}



#' Get Paris Population in 1950
#'
#' Our definition of Paris is much greater than the administrative boundary of the City of Paris in 1950.
#'
#' We manually check which surrounding communities fall within our classification of the Paris area in 1950 and
#' add the name of those communities to our list of constituting jurisdictions making up the greater Paris population.
#' Paris region candidate communes from https://fr.wikipedia.org/wiki/Réorganisation_de_la_région_parisienne_en_1964.
#' These are *all* the communes. Some will not be part of our delineation of Paris.
#' need to check manually one by one against our screen shot of the Paris area in 1950.
get_paris_pop_1950 <- function(p){

    haute_seine_92 <- c("Antony", "Asnières-sur-Seine", "Bagneux", "Bois-Colombes", "Boulogne-Billancourt", "Bourg-la-Reine", "Châtenay-Malabry", "Châtillon", "Chaville", "Clamart", "Clichy", "Colombes", "Courbevoie", "Fontenay-aux-Roses", "Garches", "La Garenne-Colombes", "Gennevilliers", "Issy-les-Moulineaux", "Levallois-Perret", "Malakoff", "Marnes-la-Coquette", "Meudon", "Montrouge", "Nanterre", "Neuilly-sur-Seine", "Le Plessis-Robinson", "Puteaux", "Rueil-Malmaison", "Saint-Cloud", "Sceaux", "Sèvres", "Suresnes", "Vanves", "Vaucresson", "Ville-d'Avray", "Villeneuve-la-Garenne")

    p92 = p %>%
        dplyr::filter(DEP == 92, year == 1954)
    p92 = p92 %>% dplyr::filter(haute_seine_92 %in% LIBGEO)
    stopifnot(nrow(p92) == length(haute_seine_92))

    # list of communes to take out of paris
    non_paris_haute_seine_92 <- c("Marnes-la-Coquette","Sèvres","Vaucresson","Ville-d'Avray")

    p92_paris <- p92 %>% dplyr::filter(!(LIBGEO %in% non_paris_haute_seine_92 ))

    val_de_marne_94 <- c("Ablon-sur-Seine", "Alfortville", "Arcueil", "Boissy-Saint-Léger", "Bonneuil-sur-Marne", "Bry-sur-Marne", "Cachan", "Champigny-sur-Marne", "Charenton-le-Pont", "Chennevières-sur-Marne", "Chevilly-Larue", "Choisy-le-Roi", "Créteil", "Fontenay-sous-Bois", "Fresnes", "Gentilly", "L'Haÿ-les-Roses", "Ivry-sur-Seine", "Joinville-le-Pont", "Le Kremlin-Bicêtre", "Limeil-Brévannes", "Maisons-Alfort", "Mandres-les-Roses", "Marolles-en-Brie", "Nogent-sur-Marne", "Noiseau", "Orly", "Ormesson-sur-Marne", "Périgny", "Le Perreux-sur-Marne", "Le Plessis-Trévise", "La Queue-en-Brie", "Rungis", "Saint-Mandé", "Saint-Maur-des-Fossés", "Saint-Maurice", "Santeny", "Sucy-en-Brie", "Thiais", "Valenton", "Villecresnes", "Villeneuve-le-Roi", "Villejuif", "Villeneuve-Saint-Georges", "Villiers-sur-Marne", "Vincennes", "Vitry-sur-Seine")

    p94 = p %>%
        dplyr::filter(DEP == 94, year == 1954)
    p94 = p94 %>% dplyr::filter(val_de_marne_94 %in% LIBGEO)
    stopifnot(nrow(p94) == length(val_de_marne_94))

    # list of communes to take out of paris
    non_paris_val_de_marne_94 <- c("Boissy-Saint-Léger","Bonneuil-sur-Marne","Limeil-Brévannes","Mandres-les-Roses","Marolles-en-Brie","Noiseau","Orly", "Ormesson-sur-Marne","Périgny", "Le Plessis-Trévise", "La Queue-en-Brie", "Rungis", "Santeny", "Sucy-en-Brie", "Valenton", "Villecresnes", "Villeneuve-le-Roi", "Villeneuve-Saint-Georges")

    p94_paris <- p94 %>% dplyr::filter(!(LIBGEO %in% non_paris_val_de_marne_94 ))


    seine_saint_denis_93 <- c("Aubervilliers", "Aulnay-sous-Bois", "Bagnolet", "Le Blanc-Mesnil", "Bobigny", "Bondy", "Le Bourget", "Clichy-sous-Bois", "Coubron", "La Courneuve", "Drancy", "Dugny", "Epinay-sur-Seine", "Gagny", "Gournay-sur-Marne", "L'Ile-Saint-Denis", "Les Lilas", "Livry-Gargan", "Montfermeil", "Montreuil", "Neuilly-Plaisance", "Neuilly-sur-Marne", "Noisy-le-Grand", "Noisy-le-Sec", "Pantin", "Les Pavillons-sous-Bois", "Pierrefitte-sur-Seine", "Le Pré-Saint-Gervais", "Le Raincy", "Romainville", "Rosny-sous-Bois", "Saint-Denis", "Saint-Ouen", "Sevran", "Stains", "Tremblay-lès-Gonesse", "Vaujours", "Villemomble", "Villepinte", "Villetaneuse")

    p93 = p %>%
        dplyr::filter(DEP == 93, year == 1954)
    nas = match(seine_saint_denis_93,p93$LIBGEO)
    not_there = seine_saint_denis_93[which(is.na(nas))]   # offenders are not there
    in_paris_seine_saint_denis_93 <- seine_saint_denis_93[which(!is.na(nas))]

    p93_paris <- p93 %>%
        dplyr::filter((LIBGEO %in% in_paris_seine_saint_denis_93 ))

    p75 = p %>%
        dplyr::filter(DEP == 75, year == 1954)

    paris_1954 <- dplyr::bind_rows(p75,p93_paris,p94_paris,p92_paris)
    return(paris_1954)
}


#' Print Paris 1950 Population to Latex
#'
#'
print_paris_pop_1950 <- function(){
    p0 = readpop()
    p = as.data.table(p0)
    x = LandUseR:::get_paris_pop_1950(p)
    ta = knitr::kable(x, format = "latex",
                      booktabs = TRUE,
                      longtable = TRUE,
                      caption = "Paris 1950 Population Classification",
                      label = "paris-1950-pop"
                      )
    cat(kableExtra::kable_styling(ta,latex_options = c("repeat_header")), file = file.path(LandUseR:::datatables(), "paris-pop-1950.tex"))
}


#' top 100 city bounding boxes
#'
#' manually defines bounding boxes for areas containing the
#' top 100 french cities in WGS84 lat-lon (google maps) CRS.
#' Also adds the centre for Paris (notre dame) manually in WGS84 lat-lon.
bboxes_top100 <- function(){
    bn = LandUseR:::get_manuals()[,list(LIBGEO,CODGEO)]
    bn[ , extent := vector("list", length = 100) ]

    bn[LIBGEO == "Lille",                  extent := raster::extent(c(2.944861, 3.298341 , 50.543319, 50.778156))]
    bn[LIBGEO == "Amiens",                 extent := raster::extent(c(2.2145981, 2.414104 ,49.847181 , 49.959225))]
    bn[LIBGEO == "Le Havre",               extent := raster::extent(c(0.033290, 0.362656 ,49.454226 , 49.558283))]
    bn[LIBGEO == "Rouen",                  extent := raster::extent(c(0.940602, 1.183368 ,49.283603 , 49.520694))]
    bn[LIBGEO == "Reims",                  extent := raster::extent(c(3.943876, 4.132186 ,49.193367 , 49.310117))]
    bn[LIBGEO == "Caen",                   extent := raster::extent(c(-0.489701, -0.270721 ,49.142144 , 49.243900))]
    bn[LIBGEO == "Metz",                   extent := raster::extent(c(6.059380, 6.270723 ,49.055194 , 49.175658))]
    bn[LIBGEO == "Nancy",                  extent := raster::extent(c(6.038759, 6.395210 ,48.574603 , 48.781997))]
    bn[LIBGEO == "Paris",                  extent := raster::extent(c(1.61, 3.08 ,48.46 , 49.118))]
    bn[LIBGEO == "Strasbourg",             extent := raster::extent(c(7.642639,7.853926 ,48.484619 , 48.672193))]
    bn[LIBGEO == "Brest",                  extent := raster::extent(c(-4.601281,-4.360628 ,48.352529 , 48.468810))]
    bn[LIBGEO == "Dunkerque",              extent := raster::extent(c(2.179341,2.517408, 50.979502, 51.056436))]
    bn[LIBGEO == "Perpignan",              extent := raster::extent(c(2.799684, 2.989531, 42.629930,42.761191))]
    bn[LIBGEO == "Rennes",                 extent := raster::extent(c(-1.778073,-1.574094 ,48.046103 , 48.173615))]
    bn[LIBGEO == "Le Mans",                extent := raster::extent(c(0.091035,0.281172 ,47.917295 , 48.062800))]
    bn[LIBGEO == "Orléans",                extent := raster::extent(c(1.767287,2.073181 ,47.796967 , 47.971889))]
    bn[LIBGEO == "Mulhouse",               extent := raster::extent(c(7.204776,7.403701 ,47.697986 ,47.842130))]
    bn[LIBGEO == "Angers",                 extent := raster::extent(c(-0.655981,-0.447919 ,47.404412 ,47.544447))]
    bn[LIBGEO == "Tours",                  extent := raster::extent(c(0.591164,0.798459 ,47.325082 ,47.483604))]
    bn[LIBGEO == "Dijon",                  extent := raster::extent(c(4.951104,5.166130 ,47.230547 ,47.384402))]
    bn[LIBGEO == "Nantes",                 extent := raster::extent(c(-1.691281,-1.419254 ,47.131828 ,47.309569))]
    bn[LIBGEO == "Limoges",                extent := raster::extent(c(1.145045,1.358150 ,45.771763 ,45.909983))]
    bn[LIBGEO == "Clermont-Ferrand",       extent := raster::extent(c(3.000352,3.223957 ,45.701680 ,45.847506))]
    bn[LIBGEO == "Lyon",                   extent := raster::extent(c(4.682204,5.052525 ,45.628088 ,45.836050))]
    bn[LIBGEO == "Saint-Étienne",          extent := raster::extent(c(4.215847,4.480892 ,45.369009 ,45.503158))]
    bn[LIBGEO == "Grenoble",               extent := raster::extent(c(5.605581,5.887174 ,45.068137 ,45.261293))]
    bn[LIBGEO == "Bordeaux",               extent := raster::extent(c(-0.830520,-0.377359 ,44.671694 ,44.967807))]
    bn[LIBGEO == "Montpellier",            extent := raster::extent(c(3.649175,4.095652 ,43.447040 ,43.701412))]
    bn[LIBGEO == "Nice",                   extent := raster::extent(c(6.861982,7.399601 ,43.522038 ,43.764086))]
    bn[LIBGEO == "Toulouse",               extent := raster::extent(c(1.216777,1.655748 ,43.397305 ,43.738317))]
    bn[LIBGEO == "Marseille",              extent := raster::extent(c(5.256408,5.606574 ,43.221860 ,43.393227))]
    bn[LIBGEO == "Toulon",                 extent := raster::extent(c(5.746459,6.189482 ,43.051562 ,43.208111))]
    bn[LIBGEO == "Pau",                    extent := raster::extent(c(-0.483260,-0.273838,43.276050 ,43.391593))]
    bn[LIBGEO == "Nîmes",                  extent := raster::extent(c(4.295740,4.417038,43.773776,43.884270))]
    bn[LIBGEO == "Besançon",               extent := raster::extent(c(5.935203,6.089627,47.205978,47.298278))]
    bn[LIBGEO == "Cherbourg-en-Cotentin",  extent := raster::extent(c(-1.708173,-1.536951,49.591829,49.671588))]
    bn[LIBGEO == "Troyes",                 extent := raster::extent(c(4.007275,4.159603,48.249062,48.337258))]
    bn[LIBGEO == "Boulogne-sur-Mer",       extent := raster::extent(c(1.560535,1.657945,50.689725,50.745426))]
    bn[LIBGEO == "Saint-Quentin",          extent := raster::extent(c(3.245462,3.345438,49.818110,49.877047))]
    bn[LIBGEO == "Lorient",                extent := raster::extent(c(-3.404595,-3.309053,47.700910,47.790783))]
    bn[LIBGEO == "Béziers",                extent := raster::extent(c(3.192872,3.294548,43.315345,43.376839))]
    bn[LIBGEO == "Calais",                 extent := raster::extent(c(1.801553,1.924938,50.916189,50.977223))]
    bn[LIBGEO == "Avignon",                extent := raster::extent(c(4.743520,4.950697,43.920633,44.015632))]
    bn[LIBGEO == "Bourges",                extent := raster::extent(c(2.340618,2.465599,47.050629,47.130391))]
    bn[LIBGEO == "Poitiers",               extent := raster::extent(c(0.283150,0.411924,46.543529,46.638859))]
    bn[LIBGEO == "Angoulême",              extent := raster::extent(c(0.090945,0.247140,45.619035,45.694662))]
    bn[LIBGEO == "Aix-en-Provence",        extent := raster::extent(c(5.381604,5.484651,43.496109,43.568511))]
    bn[LIBGEO == "Sète",                   extent := raster::extent(c(3.655132,3.775345,43.384262,43.464131))]
    bn[LIBGEO == "Bayonne",                extent := raster::extent(c(-1.586652,-1.413312,43.460803,43.547358))]
    bn[LIBGEO == "Laval",                  extent := raster::extent(c(-0.839517,-0.682987,48.052485,48.098505))]
    bn[LIBGEO == "Rochefort",              extent := raster::extent(c(-0.989032,-0.867720,45.921615,45.962694))]
    bn[LIBGEO == "Douai",                  extent := raster::extent(c(3.013404,3.156428,50.332013,50.410827))]
    bn[LIBGEO == "Montauban",              extent := raster::extent(c(1.312652,1.418559,43.962669,44.069628))]
    bn[LIBGEO == "Arras",                  extent := raster::extent(c(2.708754,2.846578,50.255199,50.320020))]
    bn[LIBGEO == "Le Creusot",             extent := raster::extent(c(4.379429,4.488826,46.777775,46.817244))]
    bn[LIBGEO == "Saint-Malo",             extent := raster::extent(c(-2.030608,-1.960790,48.613382,48.689334))]
    bn[LIBGEO == "Niort",                  extent := raster::extent(c(-0.504714,-0.358977,46.295936,46.366651))]
    bn[LIBGEO == "Vienne",                 extent := raster::extent(c(4.840998,4.941629,45.496903,45.553525))]
    bn[LIBGEO == "Valenciennes",           extent := raster::extent(c(3.432456,3.582964,50.301349,50.416407))]
    bn[LIBGEO == "Carcassonne",            extent := raster::extent(c(2.288291,2.418513,43.190463,43.239289))]
    bn[LIBGEO == "Saint-Chamond",          extent := raster::extent(c(4.484294,4.560297,45.450598,45.496264))]
    bn[LIBGEO == "Castres",                extent := raster::extent(c(2.183609,2.280618,43.580460,43.622711))]
    bn[LIBGEO == "Périgueux",              extent := raster::extent(c(0.672202,0.796224,45.174278,45.196678))]
    bn[LIBGEO == "Colmar",                 extent := raster::extent(c(7.292660,7.420742,48.061770,48.111136))]
    bn[LIBGEO == "Montluçon",              extent := raster::extent(c(2.561184,2.638652,46.313330,46.385367))]
    bn[LIBGEO == "Valence",                extent := raster::extent(c(4.848979,4.969352,44.844188,44.977663))]
    bn[LIBGEO == "Charleville-Mézières",   extent := raster::extent(c(4.671147,4.778368,49.727527,49.798895))]
    bn[LIBGEO == "Roanne",                 extent := raster::extent(c(4.019713,4.109817,46.002517,46.079755))]
    bn[LIBGEO == "Nevers",                 extent := raster::extent(c(3.115817,3.220079,46.959705,47.034447))]
    bn[LIBGEO == "Chalon-sur-Saône",       extent := raster::extent(c(4.796634,4.895622,46.734189,46.839596))]
    bn[LIBGEO == "Elbeuf",                 extent := raster::extent(c(0.989142,1.055293,49.262341,49.323525))]
    bn[LIBGEO == "Cambrai",                extent := raster::extent(c(3.190184,3.285074,50.141238,50.198118))]
    bn[LIBGEO == "Saint-Omer",             extent := raster::extent(c(2.208486,2.339291,50.701660,50.769378))]
    bn[LIBGEO == "Moulins",                extent := raster::extent(c(3.306299,3.379103,46.537155,46.602965))]
    bn[LIBGEO == "Armentières",            extent := raster::extent(c(2.833394,2.925537,50.667501,50.710963))]
    bn[LIBGEO == "Dieppe",                 extent := raster::extent(c(1.056407,1.124279,49.891384,49.941075))]
    bn[LIBGEO == "Arles",                  extent := raster::extent(c(4.582892,4.651447,43.656457,43.710221))]
    bn[LIBGEO == "Tarbes",                 extent := raster::extent(c(0.015167,0.112835,43.207913,43.269965))]
    bn[LIBGEO == "Alès",                   extent := raster::extent(c(4.052094,4.139736,44.090918,44.178667))]
    bn[LIBGEO == "La Rochelle",            extent := raster::extent(c(-1.243887,-1.080047,46.127976,46.193562))]
    bn[LIBGEO == "Blois",                  extent := raster::extent(c(1.274601,1.415218,47.541809,47.624354))]
    bn[LIBGEO == "Chartres",               extent := raster::extent(c(1.431147,1.552891,48.392358,48.489297))]
    bn[LIBGEO == "Beauvais",               extent := raster::extent(c(2.043774,2.134711,49.398221,49.457329))]
    bn[LIBGEO == "Quimper",                extent := raster::extent(c(-4.142595,-4.007847,47.958748,48.020226))]
    bn[LIBGEO == "Châlons-en-Champagne",   extent := raster::extent(c(4.298364,4.422376,48.928553,48.998451))]
    bn[LIBGEO == "Le Puy-en-Velay",        extent := raster::extent(c(3.854840,3.936414,45.009143,45.067228))]
    bn[LIBGEO == "Chambéry",               extent := raster::extent(c(5.861817,5.993300,45.540713,45.624977))]
    bn[LIBGEO == "Narbonne",               extent := raster::extent(c(2.955498,3.040298,43.149743,43.201157))]
    bn[LIBGEO == "Sedan",                  extent := raster::extent(c(4.861946,4.991722,49.664227,49.731291))]
    bn[LIBGEO == "Lisieux",                extent := raster::extent(c(0.204394,0.288547,49.123454,49.172273))]
    bn[LIBGEO == "Agen",                   extent := raster::extent(c(0.572322,0.708979,44.172293,44.213214))]
    bn[LIBGEO == "Châteauroux",            extent := raster::extent(c(1.637382,1.748112,46.769781,46.843388))]
    bn[LIBGEO == "Abbeville",              extent := raster::extent(c(1.804252,1.878326,50.090177,50.125386))]
    bn[LIBGEO == "Mâcon",                  extent := raster::extent(c(4.775973,4.848893,46.273947,46.338645))]
    bn[LIBGEO == "Albi",                   extent := raster::extent(c(2.087710,2.230343,43.885606,43.971746))]
    bn[LIBGEO == "Saumur",                 extent := raster::extent(c(-0.116165,-0.049581,47.231638,47.284655))]
    bn[LIBGEO == "Châtellerault",          extent := raster::extent(c(0.490926,0.586213,46.771402,46.860404))]
    bn[LIBGEO == "Morlaix",                extent := raster::extent(c(-3.875495,-3.779227,48.556620,48.600425))]
    bn[LIBGEO == "Vannes",                 extent := raster::extent(c(-2.812036,-2.691376,47.633021,47.703814))]
    bn[LIBGEO == "Saint-Nazaire",          extent := raster::extent(c(-2.296981,-2.112192,47.233045,47.339467))]
    bn
}

#' add center of cities to bounding boxes
#'
#' all except paris (which has center of each arrondissement instead)
merge_centers <- function(overwrite = FALSE){
    if (overwrite) {
        bb = bboxes_top100()
        raster_crs = get_raster_CRS()
        cl = sf::st_read(file.path(datadir(), "IGN-chef-lieux", "ign_metropole_adminexpress_chefs_lieux_z.shp"))
        cl = cl %>%
            filter(insee_com %in% bb[LIBGEO != "Paris", CODGEO]) %>%
            sf::st_transform(raster_crs) %>%
            select(nom_chf, insee_com)
        cc = cl %>% sf::st_coordinates()
        cl = cl %>%
            mutate(center_x = cc[ ,"X"], center_y = cc[ ,"Y"]) %>%
            sf::st_set_geometry(NULL)

        bb_cl = merge(bb, cl, by.x = "CODGEO", by.y = "insee_com", all.x = FALSE)
        paris = sf::st_as_sf(data.frame(x = 2.349020, y = 48.853481),coords = c("x","y"), crs = 4326) %>%
            sf::st_transform(raster_crs) %>%
            sf::st_coordinates()
        bb_cl = rbind(bb_cl,bb[LIBGEO == "Paris",list(LIBGEO,CODGEO,extent, center_x = paris[,"X"], center_y = paris[,"Y"],nom_chf = "Notre Dame")] )
        saveRDS(bb_cl, file.path(datadir(), "IGN-chef-lieux.RDS"))
        bb_cl

    } else {
        readRDS(file.path(datadir(), "IGN-chef-lieux.RDS"))
    }
}

# archive
