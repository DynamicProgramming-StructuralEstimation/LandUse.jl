
# `LandUseR`: Urban Extent Measurement for Coeurdacier, Oswald and Teigner

<!-- badges: start -->
<!-- badges: end -->

This `R` package contains data measurement tools used in the paper *Structural Change, Land Use and Urban Expansion* by Coeurdacier, Oswald and Teignier. In particular, this contains code for

1. Measuring city area extents and urban population via GHSL.
2. Assessing current land use outside of cities via CLC.

Basically, this is the code accompanying sections B.6 and B.7 of our online appendix. Section B.8 (Indidivual Commuting Data) uses confidential data on CASD, code for which is contained in a separate R package.

## Contents

1. Overview
2. Installation
3. Replication Data Requirements and File Directory Structure

## Overview

This R package operates our city area and population measurement exercise, as described in the online appendix. All documented functions are accessible via the Reference link in the top menu bar. Please also have a look at the Articles Menu.

## Installation

This package is not on CRAN, so you need to install directly from github via

```R
library(remotes)  # install if you don't have it
remotes::install_github("floswald/LandUse.jl", subdir = "LandUseR")
```

## Replication Data Requirements



### Data 1975 - 2015

Here the datasource is EU Commission's [Global Human Settlement](https://ghslsys.jrc.ec.europa.eu/download.php?ds=pop) project. Please check the online appendix for a thorough description.


### Data 1950 and before

* Need to resort to manual measurement of a list of French cities. 
* We use [geoportail](https://www.geoportail.gouv.fr/donnees/cartes-1950)
* We can use maps from the 1950s
* We also use maps from the 1860's produced for the Army by the Etat Major


