
# landuse

<!-- badges: start -->
<!-- badges: end -->

The goal of landuse is to measure the expansion of cities over time. 

## Data 1975 - 2015

the main datasource is EU Commission's [Global Human Settlement](https://ghslsys.jrc.ec.europa.eu/download.php?ds=pop) project. 


## Data 1950 and before

* Need to resort to manual measurement of a list of French cities. 
* We use [geoportail](https://www.geoportail.gouv.fr/donnees/cartes-1950)
* We can use maps from the 1950s
* We also use maps from the 1860's produced for the Army by the Etat Major

### How to decide whether something is part of a city or not?

We need to stick to a strict protocol to get anything useful out of this. 

1. Use mainly the "Cartes 1950" map. sometimes the arial fotos are better, or at least they can give you some good guidance about where the city ends.
1. we want to get a contiguous shape covering the city area. Often the city does *not* suddenly end. Ideally we had a grid placed over the map, and would assess the *density of built up area* in each grid cell: if there are no houses, that number is 0, if there is no unbuilt space left, that number is 100. Some maps do actually have a grid on them, so use this to make the judgement. 
1. Cutoff to be included into the city: at least 50% of all space in a grid cell is taken up by buildings/roads/etc, i.e. are *not* free open space.
1. If you make an assumption about the city limit ("I assume city x does not extend beyond that river", for example), keep that assumption in each year. All such assumptions need to be discussed.
1. Grenoble is a good example.
1. Screenshots: 
    * perform your measurement on screen. 
    * when done, make a screenshot, making sure you incldue the scale info in the bottom left corner, as well as the result of the area measurement.
    * screenshot location: all screenshots go into `LandUse/data/manual-measurement/screenshots`
    * naming convention: name each screenshot like `grenoble-EM-1.32km2.png`. 
        * `grenoble` is small caps city name
        * `EM` stands for `Etat Major`. Else we put `1950` for maps from the 1950s.
        * `1.32km2` is the measured area in square km.
1. The same procedure applies for the earlier maps, i.e the Etat Major or Cassini's maps.

### Workflow

1. Go to [geoportail](https://www.geoportail.gouv.fr/donnees/cartes-1950)
1. Search for one of th cities from below list
1. Add layer "Etat Major" (top left "cartes" switch)
1. Zoom in so entire city is visible. buildings are in red. 
1. click on the wrench top right to get tools. click on "mesures" and "mesurer une surface"
1. Delineate the city surface. clicking on the final connecting dot closes the polygon and displays the measure.
1. take a screen shot as described above
1. save as described above.
1. Add layer "Cartes 1950"
1. repeat measurement and screenshot
1. save screenshots according to naming convention above
1. Add measured areas to excel file in dropbox


* here is a list of the top 100 french cities

Rank | Commune | Department | Region | Population, 2013 | Area 1950s | Area Cassini | Area Etat Major
-- | -- | -- | -- | -- | -- | -- | --
1 | Paris | Paris | Île-de-France | 2,420,069 |   |   |  
2 | Marseille | Bouches-du-Rhône | Provence-Alpes-Côte d'Azur | 855,393 |   |   |  
3 | Lyon | Rhône | Auvergne-Rhône-Alpes | 500,715 |   |   |  
4 | Toulouse | Haute-Garonne | Occitanie | 458,298 |   |   |  
5 | Nice | Alpes-Maritimes | Provence-Alpes-Côte d'Azur | 342,295 |   |   |  
6 | Nantes | Loire-Atlantique | Pays de la Loire | 292,718 |   |   |  
7 | Strasbourg | Bas-Rhin | Grand Est | 275,718 |   |   |  
8 | Montpellier | Hérault | Occitanie | 272,084 |   |   |  
9 | Bordeaux | Gironde | Nouvelle-Aquitaine | 243,626 |   |   |  
10 | Lille | Nord | Hauts-de-France | 231,491 |   |   |  
11 | Rennes | Ille-et-Vilaine | Brittany | 211,373 |   |   |  
12 | Reims | Marne | Grand Est | 182,592 |   |   |  
13 | Le Havre | Seine-Maritime | Normandy | 172,074 |   |   |  
14 | Saint-Étienne | Loire | Auvergne-Rhône-Alpes | 172,023 |   |   |  
15 | Toulon | Var | Provence-Alpes-Côte d'Azur | 163,760 |   |   |  
16 | Grenoble | Isère | Auvergne-Rhône-Alpes | 160,215 |   |   |  
17 | Dijon | Côte-d'Or | Bourgogne-Franche-Comté | 153,003 |   |   |  
18 | Nîmes | Gard | Occitanie | 150,564 |   |   |  
19 | Angers | Maine-et-Loire | Pays de la Loire | 150,125 |   |   |  
20 | Villeurbanne | Rhône | Auvergne-Rhône-Alpes | 147,192 |   |   |  
21 | Le Mans | Sarthe | Pays de la Loire | 144,244 |   |   |  
22 | Saint-Denis | Réunion | Réunion | 142,442 |   |   |  
23 | Aix-en-Provence | Bouches-du-Rhône | Provence-Alpes-Côte d'Azur | 141,545 |   |   |  
24 | Clermont-Ferrand | Puy-de-Dôme | Auvergne-Rhône-Alpes | 141,463 |   |   |  
25 | Brest | Finistère | Brittany | 139,386 |   |   |  
26 | Limoges | Haute-Vienne | Nouvelle-Aquitaine | 135,098 |   |   |  
27 | Tours | Indre-et-Loire | Centre-Val de Loire | 134,803 |   |   |  
28 | Amiens | Somme | Hauts-de-France | 132,699 |   |   |  
29 | Perpignan | Pyrénées-Orientales | Occitanie | 120,959 |   |   |  
30 | Metz | Moselle | Grand Est | 118,634 |   |   |  
31 | Besançon | Doubs | Bourgogne-Franche-Comté | 116,952 |   |   |  
32 | Boulogne-Billancourt | Hauts-de-Seine | Île-de-France | 116,794 |   |   |  
33 | Orléans | Loiret | Centre-Val de Loire | 114,375 |   |   |  
34 | Mulhouse | Haut-Rhin | Grand Est | 112,063 |   |   |  
35 | Rouen | Seine-Maritime | Normandy | 110,755 |   |   |  
36 | Saint-Denis | Seine-Saint-Denis | Île-de-France | 109,343 |   |   |  
37 | Caen | Calvados | Normandy | 107,229 |   |   |  
38 | Argenteuil | Val-d'Oise | Île-de-France | 106,817 |   |   |  
39 | Saint-Paul | Réunion | Réunion | 104,332 |   |   |  
40 | Montreuil | Seine-Saint-Denis | Île-de-France | 104,139 |   |   |  
41 | Nancy | Meurthe-et-Moselle | Grand Est | 104,072 |   |   |  
42 | Roubaix | Nord | Hauts-de-France | 95,866 |   |   |  
43 | Tourcoing | Nord | Hauts-de-France | 93,974 |   |   |  
44 | Nanterre | Hauts-de-Seine | Île-de-France | 92,227 |   |   |  
45 | Avignon | Vaucluse | Provence-Alpes-Côte d'Azur | 90,305 |   |   |  
46 | Vitry-sur-Seine | Val-de-Marne | Île-de-France | 90,075 |   |   |  
47 | Créteil | Val-de-Marne | Île-de-France | 89,989 |   |   |  
48 | Dunkirk | Nord | Hauts-de-France | 89,882 |   |   |  
49 | Poitiers | Vienne | Nouvelle-Aquitaine | 87,427 |   |   |  
50 | Asnières-sur-Seine | Hauts-de-Seine | Île-de-France | 86,020 |   |   |  
51 | Courbevoie | Hauts-de-Seine | Île-de-France | 85,523 |   |   |  
52 | Versailles | Yvelines | Île-de-France | 85,272 |   |   |  
53 | Colombes | Hauts-de-Seine | Île-de-France | 84,577 |   |   |  
54 | Fort-de-France | Martinique | Martinique | 84,174 |   |   |  
55 | Aulnay-sous-Bois | Seine-Saint-Denis | Île-de-France | 82,634 |   |   |  
56 | Saint-Pierre | Réunion | Réunion | 81,415 |   |   |  
57 | Rueil-Malmaison | Hauts-de-Seine | Île-de-France | 79,762 |   |   |  
58 | Pau | Pyrénées-Atlantiques | Nouvelle-Aquitaine | 77,575 |   |   |  
59 | Aubervilliers | Seine-Saint-Denis | Île-de-France | 77,452 |   |   |  
60 | Le Tampon | Réunion | Réunion | 76,090 |   |   |  
61 | Champigny-sur-Marne | Val-de-Marne | Île-de-France | 75,961 |   |   |  
62 | Antibes | Alpes-Maritimes | Provence-Alpes-Côte d'Azur | 75,456 |   |   |  
63 | Béziers | Hérault | Occitanie | 74,811 |   |   |  
64 | La Rochelle | Charente-Maritime | Nouvelle-Aquitaine | 74,344 |   |   |  
65 | Saint-Maur-des-Fossés | Val-de-Marne | Île-de-France | 74,133 |   |   |  
66 | Cannes | Alpes-Maritimes | Provence-Alpes-Côte d'Azur | 73,325 |   |   |  
67 | Calais | Pas-de-Calais | Hauts-de-France | 72,520 |   |   |  
68 | Saint-Nazaire | Loire-Atlantique | Pays de la Loire | 68,513 |   |   |  
69 | Mérignac | Gironde | Nouvelle-Aquitaine | 68,386 |   |   |  
70 | Drancy | Seine-Saint-Denis | Île-de-France | 68,241 |   |   |  
71 | Colmar | Haut-Rhin | Grand Est | 67,956 |   |   |  
72 | Ajaccio | Corse-du-Sud | Corsica | 67,507 |   |   |  
73 | Bourges | Cher | Centre-Val de Loire | 67,189 |   |   |  
74 | Issy-les-Moulineaux | Hauts-de-Seine | Île-de-France | 65,662 |   |   |  
75 | Levallois-Perret | Hauts-de-Seine | Île-de-France | 65,264 |   |   |  
76 | La Seyne-sur-Mer | Var | Provence-Alpes-Côte d'Azur | 64,523 |   |   |  
77 | Quimper | Finistère | Brittany | 63,532 |   |   |  
78 | Noisy-le-Grand | Seine-Saint-Denis | Île-de-France | 62,834 |   |   |  
79 | Villeneuve-d'Ascq | Nord | Hauts-de-France | 62,616 |   |   |  
80 | Neuilly-sur-Seine | Hauts-de-Seine | Île-de-France | 62,346 |   |   |  
81 | Valence | Drôme | Auvergne-Rhône-Alpes | 61,767 |   |   |  
82 | Antony | Hauts-de-Seine | Île-de-France | 61,727 |   |   |  
83 | Cergy | Val-d'Oise | Île-de-France | 61,708 |   |   |  
84 | Vénissieux | Rhône | Auvergne-Rhône-Alpes | 61,636 |   |   |  
85 | Pessac | Gironde | Nouvelle-Aquitaine | 60,763 |   |   |  
86 | Troyes | Aube | Grand Est | 59,671 |   |   |  
87 | Clichy | Hauts-de-Seine | Île-de-France | 59,255 |   |   |  
88 | Ivry-sur-Seine | Val-de-Marne | Île-de-France | 58,933 |   |   |  
89 | Chambéry | Savoie | Auvergne-Rhône-Alpes | 58,653 |   |   |  
90 | Lorient | Morbihan | Brittany | 57,961 |   |   |  
91 | Les Abymes | Guadeloupe | Guadeloupe | 57,960 |   |   |  
92 | Montauban | Tarn-et-Garonne | Occitanie | 57,921 |   |   |  
93 | Sarcelles | Val-d'Oise | Île-de-France | 57,533 |   |   |  
94 | Niort | Deux-Sèvres | Nouvelle-Aquitaine | 57,393 |   |   |  
95 | Villejuif | Val-de-Marne | Île-de-France | 57,184 |   |   |  
96 | Saint-André | Réunion | Réunion | 56,156 |   |   |  
97 | Hyères | Var | Provence-Alpes-Côte d'Azur | 55,713 |   |   |  
98 | Saint-Quentin | Aisne | Hauts-de-France | 55,698 |   |   |  
99 | Beauvais | Oise | Hauts-de-France | 55,252 |   |   |  
100 | Épinay-sur-Seine | Seine-Saint-Denis | Île-de-France | 54,857 |   |   |  

