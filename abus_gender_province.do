*==============================================================================
*==============================================================================
*=     File-Name:      abus_gender_province.do                               == 
*=     Date:           08/11/2022                                            ==
*=     Author:         Murat Abus (mabus@tamu.edu)                           ==
*=     Purpose:        Electoral District Maps for Candidate Gender          == 
*=     Input File 1:   abus_gender_province.dta                              ==
*=     Input File 2:   tur_polbnda_adm1.shp, tur_polbnda_adm1.dbf            ==
*=     Output File:    Various Maps                                          ==
*=     Data Output:    n/a                                                   ==
*=     Previous file:   n/a                                                   ==
*=     Software:       Stata IC 16.1                                         ==
*=     Machine:        Lenovo ThinkPad E570                                  == 
*=     System:         Ubuntu 22.04.1 LTS / GNOME 42.2                       ==
*==============================================================================
*==============================================================================

/* This .do file contains the code to reproduce the chloropleth maps of distribution of characteristics relating to female candidates in the 2018 Turkish parliamentary election. These maps are in addition to the district level maps that were used in the analysis and are meant to offer an electoral district level visualization to the reader. These maps are found in Appendix E of the following article:
Abus, Murat. 2023. "A Theory of Gender's Effect on Vote SHift with a Test Based on Turkush Elections." Mediterranean Politics. DOI:http://dx.doi.org/10.1080/13629395.2023.2194154                 */ 


****************************PREPARING THE DATA*****************************

/* First, we load the shapefile and using Stata's built-in commands, spset it,
meaning we link the data and the shape together to prepare the map for putting in our data of interest in the second stage. The following code assumes that all the files are in the same folder as checked by the user with the "pwd" command. I created a file with the substantive data as "abus_gender_provimce.dta". Thıs dataset contains the number of female candidates, ranking of female candidates, and the number of female candidates in Winnable slots for the five parliamentary parties at the province (electoral district) level.                          */         
version
clear
pwd

spshape2dta tur_polbnda_adm1, replace //First-level administrative units (provinces)

//Please note that this folder only contains the shapefiles for provinces to save space. Boundaries for all administrative units in Turkey are available at "https://data.humdata.org/dataset/turkey-administrative-boundaries-levels-0-1-2" 

unzipfile "tur_polbnda_adm1.zip", replace

use "tur_polbnda_adm1.dta", clear

// Generating key variables for merging later.
sort adm1
gen id = _n 
rename adm1 pcode


/*This way, after sorting on the pcode variable, we attach an id variable that is automatically numbered from 1 upwards. These two together will be our key variables for merging the datasets below.                                     */

describe
save "tur_polbnda_adm1.dta", replace   

spset 

/* Now the shape file with the district boundaries of Turkey is linked and we can check using the user-contributed "grmap" command. Before the first usage, you may need to install it.                                                  */

grmap

// We need to change our coordinate system since Stata assumes planar as default.

spset, modify coordsys(latlong, kilometers)

save "tur_polbnda_adm1.dta", replace

/* The second stage is to prepare our substantive data for merging with the spatial data. I have created this data file to easily merge with the spatial data using the pcode of the provinces, 81 in total.                          */

use "abus_gender_province", clear

describe

// Note that I have already included the pcode and id variables in the substantive dataset for merging.   

use "tur_polbnda_adm1.dta", clear

use "abus_gender_province.dta", clear    

// The master dataset is the second one. We are including the spatial data to our dataset.
           
merge 1:1 pcode id using "tur_polbnda_adm1.dta"

//We have a merge on all 81 provinces (_merge==3) (:-)). We can drop the _merge variable.

drop _merge

//We can save the dataset with a different name to distinguish it in further analyses.
save "spatial_abus_gender_province.dta", replace

*********************************MAP CODE***********************************
grmap metro, fcolor(gs8) ndfcolor(gs16) legenda(off)

// We will make a new folder for these appendix graphs 

mkdir provincemaps

set graphics off 

//Figure 9 (Mapping Number of Female Candidates for IYI)
grmap IYIfem, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure9.gph, replace
graph export provincemaps/figure9.pdf, replace
graph export provincemaps/figure9.svg, replace

grmap IYIfem, clmethod(unique) fcolor(Reds) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure9color.gph, replace
graph export provincemaps/figure9color.pdf, replace
graph export provincemaps/figure9color.svg, replace

// Figure 10 (Mapping Number of Female Candidates for HDP)
grmap HDPfem, clmethod(unique) fcolor(Greys2) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure10.gph, replace
graph export provincemaps/figure10.pdf, replace
graph export provincemaps/figure10.svg, replace

grmap HDPfem, clmethod(unique) fcolor(Reds2) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure10color.gph, replace
graph export provincemaps/figure10color.pdf, replace
graph export provincemaps/figure10color.svg, replace

// Combining the 2 graphs into one image:

set graphics on 

graph combine provincemaps/figure9.gph provincemaps/figure10.gph, rows(2) altshrink

gr save provincemaps/combined15.gph, replace 
gr export provincemaps/combined15.pdf, replace
gr export provincemaps/combined15.svg, replace

graph combine provincemaps/figure9color.gph provincemaps/figure10color.gph, rows(2) altshrink

gr save provincemaps/combinedcolor15.gph, replace 
gr export provincemaps/combinedcolor15.pdf, replace
gr export provincemaps/combinedcolor15.svg, replace

set graphics off

// Figure 11 (Mapping Number of Female Candidates for AKP)
grmap AKPfem, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure11.gph, replace
graph export provincemaps/figure11.pdf, replace
graph export provincemaps/figure11.svg, replace

grmap AKPfem, clmethod(unique) fcolor(Reds) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure11color.gph, replace
graph export provincemaps/figure11color.pdf, replace
graph export provincemaps/figure11color.svg, replace

// Figure 12 (Mapping Number of Female Candidates for CHP)
grmap CHPfem, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure12.gph, replace
graph export provincemaps/figure12.pdf, replace
graph export provincemaps/figure12.svg, replace

grmap CHPfem, clmethod(unique) fcolor(Reds) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure12color.gph, replace
graph export provincemaps/figure12color.pdf, replace
graph export provincemaps/figure12color.svg, replace

// Figure 13 (Mapping Number of Female Candidates for MHP)
grmap MHPfem, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure13.gph, replace
graph export provincemaps/figure13.pdf, replace
graph export provincemaps/figure13.svg, replace

grmap MHPfem, clmethod(unique) fcolor(Reds) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save provincemaps/figure13color.gph, replace
graph export provincemaps/figure13color.pdf, replace
graph export provincemaps/figure13color.svg, replace

// Combining the 3 graphs into one image:

set graphics on

graph combine provincemaps/figure11.gph provincemaps/figure12.gph provincemaps/figure13.gph, rows(3) altshrink

gr save provincemaps/combined16.gph, replace 
gr export provincemaps/combined16.pdf, replace
gr export provincemaps/combined16.svg, replace

graph combine provincemaps/figure11color.gph provincemaps/figure12color.gph provincemaps/figure13color.gph, rows(3) altshrink

gr save provincemaps/combinedcolor16.gph, replace 
gr export provincemaps/combinedcolor16.pdf, replace
gr export provincemaps/combinedcolor16.svg, replace

set graphics off

// Figure 14 (Mapping Highest Position of Female Candidates fpr IYI)
grmap IYIrank, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure14.gph, replace
graph export provincemaps/figure14.pdf, replace
graph export provincemaps/figure14.svg, replace

grmap IYIrank, clmethod(unique) fcolor(Blues) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure14color.gph, replace
graph export provincemaps/figure14color.pdf, replace
graph export provincemaps/figure14color.svg, replace

// Figure 15 (Mapping Highest Position of Female Candidates for HDP)
grmap HDPrank, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure15.gph, replace
graph export provincemaps/figure15.pdf, replace
graph export provincemaps/figure15.svg, replace

grmap HDPrank, clmethod(unique) fcolor(Blues) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure15color.gph, replace
graph export provincemaps/figure15color.pdf, replace
graph export provincemaps/figure15color.svg, replace

// Combining the 2 graphs into one image:

set graphics on

graph combine provincemaps/figure14.gph provincemaps/figure15.gph, rows(2) altshrink

gr save provincemaps/combined17.gph, replace 
gr export provincemaps/combined17.pdf, replace
gr export provincemaps/combined17.svg, replace

graph combine provincemaps/figure14color.gph provincemaps/figure15color.gph, rows(2) altshrink

gr save provincemaps/combinedcolor17.gph, replace 
gr export provincemaps/combinedcolor17.pdf, replace
gr export provincemaps/combinedcolor17.svg, replace

set graphics off

// Figure 16 (Mapping Highest Position of Female Candidates for AKP),
grmap AKPrank, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure16.gph, replace
graph export provincemaps/figure16.pdf, replace
graph export provincemaps/figure16.svg, replace

grmap AKPrank, clmethod(unique) fcolor(Blues) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure16color.gph, replace
graph export provincemaps/figure16color.pdf, replace
graph export provincemaps/figure16color.svg, replace

// Figure 17 (Mapping Highest Position of Female Candidates for CHP)
grmap CHPrank, clmethod(unique) fcolor(Greys2) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure17.gph, replace
graph export provincemaps/figure17.pdf, replace
graph export provincemaps/figure17.svg, replace

grmap CHPrank, clmethod(unique) fcolor(Blues2) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure17color.gph, replace
graph export provincemaps/figure17color.pdf, replace
graph export provincemaps/figure17color.svg, replace

// Figure 18 (Mapping Highest Position of Female Candidates for MHP)
grmap MHPrank, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure18.gph, replace
graph export provincemaps/figure18.pdf, replace
graph export provincemaps/figure18.svg, replace

grmap MHPrank, clmethod(unique) fcolor(Blues) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save provincemaps/figure18color.gph, replace
graph export provincemaps/figure18color.pdf, replace
graph export provincemaps/figure18color.svg, replace

// Combining the 3 graphs into one image:

set graphics on

graph combine provincemaps/figure16.gph provincemaps/figure17.gph provincemaps/figure18.gph, rows(3) altshrink

gr save provincemaps/combined18.gph, replace 
gr export provincemaps/combined18.pdf, replace
gr export provincemaps/combined18.svg, replace

graph combine provincemaps/figure16color.gph provincemaps/figure17color.gph provincemaps/figure18color.gph, rows(3) altshrink

gr save provincemaps/combinedcolor18.gph, replace 
gr export provincemaps/combinedcolor18.pdf, replace
gr export provincemaps/combinedcolor18.svg, replace

set graphics off

// Figure 19 (Mapping Number of Winnable Slots for Women - IYI)
grmap IYIwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure19.gph, replace
graph export provincemaps/figure19.pdf, replace
graph export provincemaps/figure19.svg, replace

grmap IYIwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure19color.gph, replace
graph export provincemaps/figure19color.pdf, replace
graph export provincemaps/figure19color.svg, replace

// Figure 20 (Mapping Number of Winnable Slots for Women - HDP)
grmap HDPwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure20.gph, replace
graph export provincemaps/figure20.pdf, replace
graph export provincemaps/figure20.svg, replace

grmap HDPwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure20color.gph, replace
graph export provincemaps/figure20color.pdf, replace
graph export provincemaps/figure20color.svg, replace

// Combining the 2 graphs into one image:

set graphics on

graph combine provincemaps/figure19.gph provincemaps/figure20.gph, rows(2) altshrink

gr save provincemaps/combined19.gph, replace 
gr export provincemaps/combined19.pdf, replace
gr export provincemaps/combined19.svg, replace

graph combine provincemaps/figure19color.gph provincemaps/figure20color.gph, rows(2) altshrink

gr save provincemaps/combinedcolor19.gph, replace 
gr export provincemaps/combinedcolor19.pdf, replace
gr export provincemaps/combinedcolor19.svg, replace

set graphics off 

// Figure 21 (Mapping Number of Winnable Slots for Women - AKP)
grmap AKPwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure21.gph, replace
graph export provincemaps/figure21.pdf, replace
graph export provincemaps/figure21.svg, replace

grmap AKPwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure21color.gph, replace
graph export provincemaps/figure21color.pdf, replace
graph export provincemaps/figure21color.svg, replace

// Figure 22 (Mapping Number of Winnable Slots for Women - CHP)
grmap CHPwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure22.gph, replace
graph export provincemaps/figure22.pdf, replace
graph export provincemaps/figure22.svg, replace

grmap CHPwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure22color.gph, replace
graph export provincemaps/figure22color.pdf, replace
graph export provincemaps/figure22color.svg, replace

// Figure 23 (Mapping Number of Winnable Slots for Women - MHP)
grmap MHPwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure23.gph, replace
graph export provincemaps/figure23.pdf, replace
graph export provincemaps/figure23.svg, replace

grmap MHPwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save provincemaps/figure23color.gph, replace
graph export provincemaps/figure23color.pdf, replace
graph export provincemaps/figure23color.svg, replace

// Combining the 3 graphs into one image:

set graphics on

graph combine provincemaps/figure21.gph provincemaps/figure22.gph provincemaps/figure23.gph, rows(3) altshrink

gr save provincemaps/combined20.gph, replace 
gr export provincemaps/combined20.pdf, replace
gr export provincemaps/combined20.svg, replace

graph combine provincemaps/figure21color.gph provincemaps/figure22color.gph provincemaps/figure23color.gph, rows(3) altshrink

gr save provincemaps/combinedcolor20.gph, replace 
gr export provincemaps/combinedcolor20.pdf, replace
gr export provincemaps/combinedcolor20.svg, replace

****************************END of MAP CODE***********************************
****************************END of DO FILE************************************ 
