
# `LandUseR`: Urban Extent Measurement for Coeurdacier, Oswald and Teigner

<!-- badges: start -->
<!-- badges: end -->

This `R` package contains data measurement tools used in the paper *Structural Change, Land Use and Urban Expansion* by Coeurdacier, Oswald and Teignier. In particular, this contains code for

1. Measuring city area extents and urban population via GHSL.
2. Assessing current land use outside of cities via CLC.

Basically, this is the code accompanying sections B.6 and B.7 of our online appendix. Please refer to this document for precise instructions and links to the relevant data downloads (we provide a full data package below). Section B.8 (Individual Commuting Data) uses confidential data on CASD, code for which is contained in a separate R package.

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

In order to use the code in this package to replicate our results you need to download our data folder from dropbox (2.5GB). This contains all the data needed. Luckily there are no license restrictions and you are free to use the data as is. Then there are 2 steps:

1. you need to tell this package the location of the downloaded data on your computer via an environment variable. Simple instructions can be followed (after you installed the package) via

```R
library(LandUseR)
dboxdir()   # will throw an error with instructions
```

2. you need to create an `output` directory with 2 subdirectories in the same location that contains the data from point 1. That is, your computer needs to contain something like this:

```
LandUse  (you set an environment variable in point 1. to this location)
 |
 |-- data  (you downloaded this from dropbox)
 |-- output
       |
       |-- plots
       |-- tables
```       
