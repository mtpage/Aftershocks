# Aftershocks
MatLab code to find Omori parameters within global tectonic regions

RunAftershockCodes.m recreates Figure 5 from the Page et al. global Omori parameter paper, currently submitted to BSSA.  This calculates Omori a-values and p-values within global tectonic regions ("strec" regions, see below), which are modifications of the Flinn-Engdahl regions (see Garcia et al. (2012)).

This code uses the other MatLab codes in this repository, as well as the MatLab data file NEIC_Catalog_1990-2015.mat.  This data file contains the NEIC catalog from 1990 to 2015 in 10-column format (year / month / day / hour / minute / sec / lat / long / depth / magnitude) and the strec region for each earthquake in the catalog.  These strec regions are the "FE Seismotectonic Domain" regions given by the NEIC strec code, which is available here: https://github.com/usgs/strec
