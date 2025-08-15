# Fluctuating asymmetry analysis in bee wings points towards developmental health issues based on land cover

## Abstract

In bilaterally symmetrical organisms, fluctuating asymmetry indicates developmental health issues typically caused by the environment of the organism. Within bees, fluctuating asymmetry of the wing has been used to indicate health defects based on environmental gradients such as areas with urban to forested land cover, with the interpretation that developed areas –due to the urban heat island effect, habitat fragmentation, and greater pollutants– will cause greater developmental defects to bees. Typical practice employs use of landmarking wing junctions to represent asymmetry, with the formula (|R - L|) as a metric to find asymmetry measurements. However, this linear method results in limitations on finding asymmetry from circularity, curvature, and area of the wings. This project attempts a more comprehensive method, which uses the outline of the wing cells themselves, saved as vector files, to directly compare between the left and right wing of the bee. It involves tracing wing cells from the bee species Panurginus Polytrichus, collected in three different land cover sites –Developed, Agricultural, and Forest– to analyze whether or not those sites will result in differing fluctuating asymmetry values, thereby linking developmental health defects to specific regions. The species Panurginus Polytrichus was chosen due to the high number of individuals collected within the dataset used for the experiment. Choosing a single species of bee allows for more control in variation within the experiment. This exploratory study so far covers the methods of obtaining bee wing asymmetry values, through tracing and plotting the data in R, as well as visualizing the data as a fluctuating asymmetry value through regression analysis and Normalized Root Mean Square Error. My hypothesis is that if a greater overall fluctuating asymmetry value is found in developed sites, then the environments of those developed sites cause birth defects in bees, leading to a lower quality of life overall.

### Wing Data Description Table

#### In Data/Data
|Data|Description|
|:---|:---|
|"occurrences.csv"|Initial dataset of bee wings traced for project; includes wing side, wing cell, land cover, and asymmetry variables.|
|"occurrences_complete.csv"|Complete dataset of all bee wings traced for project; includes wing side, wing cell, land cover, and asymmetry variables.|

#### In Data/Output
|Data|Description|
|:---|:---|
|"Rplot.png|Plot consisting of the bee wings grouped by land cover, measured over percent overlap value. Percent overlap indicating the amount of overlap between left and right wing cell, for those with greater than 80% overlap.|
|"NRMSE by Land Cover & Wing Cell.png"|Normalized root mean square error of variables used for wing cell asymmetry measurement, grouped by land cover and divided by wing cell. The higher the value the greater the asymmetry.|
|"Regression by Land Cover & Wing Cell.png"|Coefficient of determination values of variables used for wing cell asymmetry measurement, grouped by land cover and divided by wing cell. A higher value indicates more symmetry between wing cells.|

### Rstudio Code Description Table

#### In Data/Script/Archive
##### Version history of my attempts to code for fluctuating asymmetry with SVG files.
|Data|Description|
|:---|:---|
|"beedata_test.R"|Intitial testing of code for project. Converts svgs to bitmaps and attempts to read their xml data. Never got to elliptical fourier analysis.|
|"beedata_test_momocs.R"|Analysis of wing cells using Momocs R Package. Converts svg files into bitmaps and reconstructs them as ellipses using EFA harmonics. Attempts to find asymmetry from overlapping cell area. Not sufficient for analysis of asymmetry.|


#### In Data/Script
##### Code used in transforming SVG data into EFA objects.
|Data|Description|
|:---|:---|
|"beedata_test_svgparser.R"|Code for inputting SVG file into Rstudio, replicating it as an Elliptical Fourier Object using harmonics, and Procrustes translating it to its opposite side pair. Also computes variables for fluctuating asymmetry including area, circularity, length, width, and percent overlap.|
|"beedata_test_forloop.R"|Code of working for loop to parse and reconstruct all SVG files in a directory as EFA objects, and then write all of their asymmetry variables into the occurrences csv document.|
|"beedata_test_forloop2.R"|Code of working for loop, same as "beedata_test_forloop.R" but with added code to rotate svg files beyond Procrustes translation, as the initial code did not properly align the two wings every time.|

## General workflow for running analyses

----------------------------------------------------------------------------
1. Download full Panurginus Polytrichus svg folder zip file from Zenodo.
2. Download occurrences_complete.csv from Github
3. Open beedata_test_forloop.R and scroll down to "non-flipping code." run code using svg folder as directory
     1. check percent overlap to determine if any of the wing pairs were improperly aligned.
     2. Check for proper alignment directly by changing directory to the svg file with low percent overlap, run "testing orientation of wings" and plot efi of both wings over each other
5. Save dataset and open beedata_graphs.R
6. Filter dataset to not include Notes., Total.Complete.Agr, Total.Complete.Dev, and Total.Complete.For
7. Plot graphs of the data you wish to view

#### Wing Vector Shape Analysis

Contains all data of wing cell svg files.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16619668.svg)](https://doi.org/10.5281/zenodo.16619668)

#### Momocs Package Used for EFA

Package for EFA construction in R. Useful to look into this package to discern more streamlined ways of coding Fourier recreations in R.
https://momx.github.io/Momocs/


