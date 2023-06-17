*==============================================================================
*==============================================================================
*=     File-Name:      abus_medit_politics.do                                == 
*=     Date:           01/03/2023                                            ==
*=     Author:         Murat Abus (mabus@tamu.edu)                           ==
*=     Purpose:        Replicating Abus' Mediterranean Politics Paper        == 
*=     Input File 1:   abus_incongruity_dataset.dta                          ==
*=     Input File 2:   tur_polbnda_adm2.shp, tur_polbnda_adm2.dbf            ==
*=     Output File:    incongruity.png                                       ==
*=     Data Output:    n/a                                                   ==
*=     Previous file:   n/a                                                   ==
*=     Software:       Stata IC 16.1                                         ==
*=     Machine:        Lenovo ThinkPad E570                                  == 
*=     System:         Ubuntu 22.04.1 LTS / GNOME 42.5                       ==
*==============================================================================
*==============================================================================

/* This .do file contains code to replicate the spatial analysis as reported in: 

Abus, Murat. 2023. "A Theory of Gender's Effect on Vote Shift with a Test based on Turkish Elections." Mediterranean Politics DOI:http://dx.doi.org/10.1080/13629395.2023.2194154                           */

// Working directory should be set to the location of the files for replication.

// cd "/......../" Set your working directory

****************************PREPARING THE DATA*****************************

/* First, we load the shapefile and using Stata's built-in commands, spset it,
meaning we link the data and the shape together to prepare the map for putting in our data of interest in the second stage. The following code assumes that all the files are in the same folder as checked by the user with the "pwd" command. I created a file with the substantive data as "abus_incongruity_dataset.dta"                                               */

version 
clear
pwd

unzipfile "tur_polbnda_adm2.zip", replace

spshape2dta tur_polbnda_adm2, replace //Second-level administrative units (districts)

//Please note that this folder only contains the shapefiles for districts to save space. Boundaries for all administrative units in Turkey are available at "https://data.humdata.org/dataset/turkey-administrative-boundaries-levels-0-1-2" 
use "tur_polbnda_adm2.dta", clear

//Generating key variables for merging later.
encode pcode, generate(dcode)
sort dcode
gen id = _n 
drop pcode
rename dcode pcode

/*This way, after sorting on the dcode variable, we attach an id variable that is automatically numbered from 1 upwards. These two together will be our key variables for merging the datasets below.                                     */

describe
save "tur_polbnda_adm2.dta", replace   

spset 

/* Now the shape file with the district boundaries of Turkey is linked and we can check using the user-contributed "grmap" command. Before the first usage, you may need to install it.                                                  */

grmap

// We need to change our coordinate system since Stata assumes planar as default.

spset, modify coordsys(latlong, kilometers)

save "tur_polbnda_adm2.dta", replace

/* The second stage is to prepare our substantive data for merging with the spatial data. I have created this data file to easily merge with the spatial data using the pcode of the districts, 973 in total.                          */

use "abus_incongruity_dataset", clear

describe

// Note that I have already included the pcode and id variables in the substantive dataset for merging.   

use "tur_polbnda_adm2.dta", clear

use "abus_incongruity_dataset.dta", clear    

// The master dataset is the second one. We are including the spatial data to our dataset.
           
merge 1:1 pcode id using "tur_polbnda_adm2.dta"

//We have a merge on all 973 districts (_merge==3). We can drop the _merge variable.

drop _merge

//We can save the dataset with a different name to distinguish it in further analyses.
save "spatial_abus_incongruity_dataset.dta", replace

// Summary Statistics in Table 4 in the paper:
summarize educgengap sedi IYI akpdiff chpdiff mhpdiff hdpdiff margin lnpop ///
dideprat1 dibulge1 mosquedens   

// This gives the descriptive statistics about female candidate rankings
summarize AKPrank CHPrank MHPrank IYIrank HDPrank
summarize AKPfem CHPfem MHPfem IYIfem HDPfem
summarize AKPwinslot CHPwinslot MHPwinslot IYIwinslot HDPwinslot

/* Now that our spatial dataset is ready, we can use it to first map variables of interest as used in the paper. The code below will reproduce all the figures of the paper:                                                                */

**************************BEGINNING OF MAP CODE****************************

//Figure 1 (Mapping the Dependent Variable):
grmap IYI, clmethod(eqint) clnumber(4) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Vote Share) 
graph save figure1.gph, replace
graph export figure1.pdf, replace
graph export figure1.svg, replace

grmap IYI, clmethod(eqint) clnumber(4) fcolor(Blues) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Vote Share) 
graph save figure1color.gph, replace
graph export figure1color.pdf, replace
graph export figure1color.svg, replace

//Figure 2 (Mapping the Main Independent Variable)
grmap educgengap, clmethod(custom) clnumber(5) clbreaks(-.127 .15 .30 .45 .764 1.06) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Educational Gender Gap) legcount
graph save figure2.gph, replace
graph export figure2.pdf, replace
graph export figure2.svg, replace

grmap educgengap, clmethod(custom) clnumber(5) clbreaks(-.127 .15 .30 .45 .764 1.06) fcolor(Reds) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Educational Gender Gap) legcount
graph save figure2color.gph, replace
graph export figure2color.pdf, replace
graph export figure2color.svg, replace

//Figure 3 (Mapping Vote Change of Parliamentary Parties):

mkdir graphs

set graphics off //We will combine the four graphs into one as in the paper.

// Graph 1: Change in AKP Vote Share, 2015-2018
grmap akpdiff, clmethod(custom) clnumber(5) clbreaks(-36.5 -20 -10 0 5 12.9) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Change in Vote Share) legcount

graph save graphs/graph1.gph, replace

*****Saving Graph 1 in Color*****
grmap akpdiff, clmethod(custom) clnumber(4) clbreaks(-36.5 -15 0 5 12.9) fcolor(RdYlGn) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Change in Vote Share) legcount

graph save graph1color.gph, replace
graph export graph1color.pdf, replace
graph export graph1color.svg, replace

// Graph 2: Change in CHP Vote Share, 2015-2018
grmap chpdiff, clmethod(custom) clnumber(5) clbreaks(-18.6 -10 -5 0 5 11.3) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Change in Vote Share) legcount

graph save graphs/graph2.gph, replace

*****Saving Graph 2 in Color*****
grmap chpdiff, clmethod(custom) clnumber(4) clbreaks(-18.6 -10 0 5 11.3) fcolor(RdYlGn) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Change in Vote Share) legcount

graph save graph2color.gph, replace
graph export graph2color.pdf, replace
graph export graph2color.svg, replace

// Graph 3: Change in MHP Vote Share, 2015-2018
grmap mhpdiff, clmethod(custom) clnumber(5) clbreaks(-23.9 -10 0 5 10 32.3) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Change in Vote ShareVote Share) legcount

graph save graphs/graph3.gph, replace 

*****Saving Graph 3 in Color*****
grmap mhpdiff, clmethod(custom) clnumber(4) clbreaks(-23.9 -10 0 15 32.3) fcolor(RdYlGn) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Change in Vote Share) legcount

graph save graph3color.gph, replace
graph export graph3color.pdf, replace
graph export graph3color.svg, replace

// Graph 4: Change in HDP Vote Share, 2015-2018
grmap hdpdiff, clmethod(custom) clnumber(5) clbreaks(-22.2 -10 0 1 5 18.1) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Change in Vote Share) legcount

graph save graphs/graph4.gph, replace 

*****Saving Graph 4 in Color*****
grmap hdpdiff, clmethod(custom) clnumber(4) clbreaks(-22.2 -10 0 6 18.1) fcolor(RdYlGn) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Change in Vote Share) legcount

graph save graph4color.gph, replace
graph export graph4color.pdf, replace
graph export graph4color.svg, replace

// Combining the 4 graphs into one image:

set graphics on 

graph combine graphs/graph1.gph graphs/graph2.gph graphs/graph3.gph ///
graphs/graph4.gph, rows(4) altshrink

gr save graphs/combined.gph, replace 
gr export graphs/combined.pdf, replace
gr export graphs/combined.svg, replace

graph combine graph1color.gph graph2color.gph graph3color.gph ///
graph4color.gph, rows(4) altshrink

gr save combined.gph, replace 
gr export combined.pdf, replace
gr export combined.svg, replace

*****************************END OF MAP CODE********************************

*****************************THE MODELS AND THE TABLES***********************

// The working directory needs to be where all files generated through the dataset generation steps are located.

// Creating the weights matrix as contiguity matrix for spatial regression:
 
spmatrix create contiguity W if IYI !=., normalize(row)

// Checking for spatial dependence (Moran's I):
regress IYI 
estat moran, errorlag(W)
regress HDP
estat moran, errorlag(W)
/* We reject the null hypothesis that error is i.i.d., and there is evidence of 
spatial dependence and utilization of spatial model is justified.             */

// Main Model for IYI  

spregress IYI c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress IYI c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress IYI c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress IYI c.hdpdiff##c.educgengap HDP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// For the table 
eststo: quietly spregress IYI c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress IYI c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress IYI c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress IYI c.hdpdiff##c.educgengap HDP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

esttab using "table3.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("AKP CHP MHP HDP") title("IYI Party Vote Share in 2018")

// Main Model for HDP 

spregress HDP HDP15 c.educgengap##c.HDPrank AKPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.educgengap##c.HDPrank CHPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.educgengap##c.HDPrank MHPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.educgengap##c.HDPrank IYIrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// For the table
eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPrank AKPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPrank CHPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPrank MHPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPrank IYIrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

esttab using "table4.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("AKP CHP MHP IYI") title("HDP Vote Share in 2018")

// The following is for the related model to explain the HDP vote share through the number of female candidates in the Party List Compared to Other Parties 

spregress HDP HDP15 c.educgengap##c.HDPwinslot AKPwinslot mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.educgengap##c.HDPwinslot CHPwinslot mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.educgengap##c.HDPwinslot MHPwinslot mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.educgengap##c.HDPwinslot IYIwinslot mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// For the table
eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPwinslot AKPwinslot mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPwinslot CHPwinslot mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPwinslot MHPwinslot mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPwinslot IYIwinslot mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

esttab using "table55.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("AKP CHP MHP IYI") title("HDP Vote Share in 2018")


***************************END OF CODE FOR MODELS****************************

*******************************MARGINS ANALYSES*******************************
// Note: The bulk of the graph code is adapted from Garcia and Wimpy (2016).
// 1. IYI Party Marginal Effects Plot

mkdir marginplots 

set graphics off

qui spregress IYI c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

// Graph one: Mean educational gender gap -1 s.d. moderated by change in CHP Vote Share

qui margins, at(chpdiff=(-18(1)0) educgengap=(0.0846255)) 


// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(5 10 15 20 25 30, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-18 -15 -12 -9 -6 -3 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(a, size(large) color(gs6)) ///
ylabel(, labsize(small) noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(5)30) title("") ytitle("") 

gr save marginplots/chp1.gph, replace

// Graph two: Mean educational gender gap moderated by change in CHP Vote Share

qui margins, at(chpdiff=(-18(1)0) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(5 10 15 20 25 30, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-18 -15 -12 -9 -6 -3 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(b, size(large) color(gs6)) ///
ylabel(, labsize(small) noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(5)30) title("") ytitle("") 

gr save marginplots/chp2.gph, replace

// Graph three: Mean educational gender gap + 1 s.d. moderated by change in CHP Vote Share 

qui margins, at(chpdiff=(-18(1)0) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(5 10 15 20 25 30, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-18 -15 -12 -9 -6 -3 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(c, size(large) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(5)30) title("") ytitle("") 

gr save marginplots/chp3.gph, replace
			

qui spregress IYI c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

// Graph four: Mean educatinal gender gap -1 s.d. moderated by change in MHP Vote Share

qui margins, at(mhpdiff=(-24(1)0) educgengap=(0.0846255)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(5 10 15 20 25 30, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-24 -20 -16 -12 -8 -4 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(d, size(large) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(5)30) xlabel(-24(4)0) title("") ytitle("") 

gr save marginplots/mhp1.gph, replace

// Graph five: Mean educational gender gap moderated by change in MHP Vote Share

qui margins, at(mhpdiff=(-24(1)0) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(5 10 15 20 25 30, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-24 -20 -16 -12 -8 -4 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(e, size(large) color(gs6)) ///
ylabel(, labsize(small) noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(5)30) xlabel(-24(4)0) title("") ytitle("") 

gr save marginplots/mhp2.gph, replace

// Graph six: Mean educational gender gap + 1 s.d. moderated by change in MHP Vote Share 

qui margins, at(mhpdiff=(-24(1)0) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(5 10 15 20 25 30, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-24 -20 -16 -12 -8 -4 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(f, size(large) color(gs6)) ///
ylabel(, labsize(small) noticks) yscale(noline) xsize(4) ysize(4) legend(off)  ylabel(0(5)30) xlabel(-24(4)0) title("") ytitle("") 

gr save marginplots/mhp3.gph, replace

// Turn the graphics back on for combined graph:

set graphics on			

// Now to combine the graphs into one image:

graph combine marginplots/chp1.gph marginplots/chp2.gph marginplots/chp3.gph ///
marginplots/mhp1.gph marginplots/mhp2.gph marginplots/mhp3.gph, ///
graphregion(fcolor(white)) l1(Predicted Vote Share of IYI Party (2018), size(medium)) b1(Change in Vote Shares of CHP(a-c) and MHP(d-f) (2015-2018), size(medium)) title(Mean Educational Gender Gap ± 1 s.d., size(medium)) 

gr save marginplots/marginsiyi.gph, replace 

gr export marginplots/marginsiyi.pdf, replace
gr export marginplots/marginsiyi.svg, replace

// 2. HDP Marginal Effect Plots
// There will be two marginal effects plots for HDP to map the interaction terms from the econometric analyses presented. The first will show the interaction terms from the three models presented in Table 6 of the paper (interaction terms are significant across the three models) while the second will show the interaction of educational gender gap and the highest position of female candidates fielded by HDP in the 2018 election. 

// 2a

set graphics off

qui spregress HDP HDP15 c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// Graph 1: Mean educational gender gap -1 s.d. moderated by change in AKP vote share

qui margins, at(akpdiff=(-36(1)0) educgengap=(0.0846255)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(8 10 12 14, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-36 -30 -24 -18 -12 -6 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(a, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(2.5) ysize(2.5) legend(off) ylabel(8(2)14) title("") ytitle("") 

gr save marginplots/hdp1.gph, replace
			

// Graph 2: Mean educational gender gap noderated by change in AKP vote share

qui margins, at(akpdiff=(-36(1)0) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(8 10 12 14, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-36 -30 -24 -18 -12 -6 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(b, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(2.5) ysize(2.5) legend(off) ylabel(8(2)14) title("") ytitle("") 

gr save marginplots/hdp2.gph, replace

// Graph 3: Mean educational gender gap +1 s.d. moderated by change in AKP vote share

qui margins, at(akpdiff=(-36(1)0) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(8 10 12 14, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-36 -30 -24 -18 -12 -6 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(c, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(2.5) ysize(2.5) legend(off) ylabel(8(2)14) title("") ytitle("") 

gr save marginplots/hdp3.gph, replace

qui spregress HDP HDP15 c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// Graph 4: Mean educational gender gap -1 s.d. moderated by change in CHP vote share

qui margins, at(chpdiff=(-18(3)0) educgengap=(0.0846255)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(8 10 12 14, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-18 -15 -12 -9 -6 -3 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(d, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(2.5) ysize(2.5) legend(off) ylabel(8(2)14) title("") ytitle("") 

gr save marginplots/hdp4.gph, replace

// Graph 5: Mean educational gender gap moderated by change in CHP vote share

qui margins, at(chpdiff=(-18(3)0) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(8 10 12 14, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-18 -15 -12 -9 -6 -3 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(e, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(2.5) ysize(2.5) legend(off) ylabel(8(2)14) title("") ytitle("") 

gr save marginplots/hdp5.gph, replace

// Graph 6: Mean educational gender gap +1 s.d. moderated by change in CHP vote share

qui margins, at(chpdiff=(-18(3)0) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(8 10 12 14, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-18 -15 -12 -9 -6 -3 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(f, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(2.5) ysize(2.5) legend(off) ylabel(8(2)14) title("") ytitle("") 

gr save marginplots/hdp6.gph, replace

qui spregress HDP HDP15 c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// Graph 7: Mean educational gender gap -1 s.d. moderated by change in MHP vote share

qui margins, at(mhpdiff=(-24(1)0) educgengap=(0.0846255)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(8 10 12 14, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-24 -20 -16 -12 -8 -4 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(g, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(2.5) ysize(2.5) legend(off) ylabel(8(2)14) title("") ytitle("") 

gr save marginplots/hdp7.gph, replace

// Graph 8: Mean educational gender gap moderated by change in MHP vote share

qui margins, at(mhpdiff=(-24(1)0) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(8 10 12 14, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-24 -20 -16 -12 -8 -4 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(h, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(2.5) ysize(2.5) legend(off) ylabel(8(2)14) title("") ytitle("") 

gr save marginplots/hdp8.gph, replace

// Graph 9: Mean educational gender gap +1 s.d. moderated by change in MHP vote share

qui margins, at(mhpdiff=(-24(1)0) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(8 10 12 14, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(-24 -20 -16 -12 -8 -4 0, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(i, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(2.5) ysize(2.5) legend(off) ylabel(8(2)14) title("") ytitle("") 

gr save marginplots/hdp9.gph, replace

// Turn the graphics back on for combined graph:

set graphics on			

// Now to combine the graphs into one image:

graph combine marginplots/hdp1.gph marginplots/hdp2.gph marginplots/hdp3.gph ///
marginplots/hdp4.gph marginplots/hdp5.gph marginplots/hdp6.gph ///  
marginplots/hdp7.gph marginplots/hdp8.gph marginplots/hdp9.gph, rows(3) ///
graphregion(fcolor(white)) l1(Predicted Vote Share of HDP (2018), size(medium)) b1(Change in Vote Shares of AKP(a-c) CHP(d-f) and MHP(g-i) (2015-2018), size(medium)) title(Mean Educational Gender Gap ± 1 s.d., size(medium)) 

gr save marginplots/marginshdp.gph, replace 
gr export marginplots/marginshdp.pdf, replace
gr export marginplots/marginshdp.svg, replace

// 2b
// The marginal effects are calculated from Models I, II, and IV from Table 8, where the interaction of educational gender gap and highest position of female candidate placed by HDP is significant.
mkdir marginplots
set graphics off

qui spregress HDP HDP15 c.educgengap##c.HDPrank IYIrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// Graph 1: From Model IV (controlling for IYI candidate placement) with mean educational gender gap -1 s.d.

qui margins, at(HDPrank=(1(0.1)5) educgengap=(0.0846255)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(1 2 3 4 5, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(a, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(3.5) ysize(3.5) legend(off) ylabel(6(2)12) title("") ytitle("") 

gr save marginplots/hdprank1.gph, replace

// Graph 2: From Model IV (controlling for IYI candidate placement) with mean educational gender gap 

qui margins, at(HDPrank=(1(0.1)5) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(1 2 3 4 5, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(b, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(3.5) ysize(3.5) legend(off) ylabel(6(2)12) title("") ytitle("") 

gr save marginplots/hdprank2.gph, replace

// Graph 3: From Model IV (controlling for IYI candidate placement) with mean educational gender gap +1 s.d.

qui margins, at(HDPrank=(1(0.1)5) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(1 2 3 4 5, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(c, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(3.5) ysize(3.5) legend(off) ylabel(6(2)12) title("") ytitle("") 

gr save marginplots/hdprank3.gph, replace

qui spregress HDP HDP15 c.educgengap##c.HDPrank AKPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// Graph 4: From Model I (controlling for AKP candidate placement) with mean educational gender gap -1 s.d.

qui margins, at(HDPrank=(1(0.1)5) educgengap=(0.0846255)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(1 2 3 4 5, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(d, size(large) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(3.5) ysize(3.5) legend(off) ylabel(6(2)12) title("") ytitle("") 

gr save marginplots/hdprank4.gph, replace

// Graph 5: From Model I (controlling for AKP candidate placement) with mean educational gender gap 

qui margins, at(HDPrank=(1(0.1)5) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(1 2 3 4 5, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(e, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(3.5) ysize(3.5) legend(off) ylabel(6(2)12) title("") ytitle("") 

gr save marginplots/hdprank5.gph, replace

// Graph 6: From Model I (controlling for AKP candidate placement) with mean educational gender gap +1 s.d.

qui margins, at(HDPrank=(1(0.1)5) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(1 2 3 4 5, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(f, size(large) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(3.5) ysize(3.5) legend(off) ylabel(6(2)12) title("") ytitle("") 

gr save marginplots/hdprank6.gph, replace

qui spregress HDP HDP15 c.educgengap##c.HDPrank CHPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// Graph 7: From Model II (controlling for CHP candidate placement) with mean educational gender gap -1 s.d.

qui margins, at(HDPrank=(1(0.1)5) educgengap=(0.0846255)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(1 2 3 4 5, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(g, size(large) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(3.5) ysize(3.5) legend(off) ylabel(6(2)12) title("") ytitle("") 

gr save marginplots/hdprank7.gph, replace

// Graph 8: From Model II (controlling for CHP candidate placement) with mean educational gender gap 

qui margins, at(HDPrank=(1(0.1)5) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(1 2 3 4 5, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(h, size(medium) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(3.5) ysize(3.5) legend(off) ylabel(6(2)12) title("") ytitle("") 

gr save marginplots/hdprank8.gph, replace

// Graph 9: From Model I (controlling for CHP candidate placement) with mean educational gender gap +1 s.d.

qui margins, at(HDPrank=(1(0.1)5) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(1 2 3 4 5, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(i, size(large) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(3.5) ysize(3.5) legend(off) ylabel(6(2)12) title("") ytitle("") 

gr save marginplots/hdprank9.gph, replace

// Turn the graphics back on for combined graph:

set graphics on			

// Now to combine the graphs into one image:

graph combine marginplots/hdprank1.gph marginplots/hdprank2.gph marginplots/hdprank3.gph marginplots/hdprank4.gph marginplots/hdprank5.gph marginplots/hdprank6.gph marginplots/hdprank7.gph marginplots/hdprank8.gph marginplots/hdprank9.gph, rows(3) ///
graphregion(fcolor(white)) l1(Predicted Vote Share of HDP (2018), size(medium)) /// 
b1(`"Highest Placement of HDP Female Candidate"' `"versus IYI (a-c), AKP(d-f) and CHP(g-i)"' , size(medium)) title(Mean Educational Gender Gap ± 1 s.d., size(medium)) 

gr save marginplots/marginshdprank.gph, replace 

gr export marginplots/marginshdprank.pdf, replace
gr export marginplots/marginshdprank.svg, replace


// End of Main Replication Files

********************************************************************************
**************************** ROBUSTNESS CHECKS *********************************
// 1. Using Educational Gender Gap as the main independent variable as in the model, but not modelling the spatial aspect of the data to demonstrate that the results are not driven by spatial modelling.

// Use the dataset without the weight matrix.
use "spatial_abus_incongruity_dataset.dta", clear

reg IYI c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

reg IYI c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

reg IYI c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

reg IYI c.hdpdiff##c.educgengap HDP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

// For presenting the robustness check results in the same format as the main model:

eststo: quietly reg IYI c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

eststo: quietly reg IYI c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

eststo: quietly reg IYI c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

eststo: quietly reg IYI c.hdpdiff##c.educgengap HDP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

esttab using "table6.tex", replace b(%9.3g) se r2 label nogap onecell ///
nonumbers mtitle("AKP CHP MHP HDP") title("IYI Party Vote Share in 2018 (Robustness Check)")

// 2. Running the model with a spatial lag specification instead of a spatial error specification. Although theoretically the spatial error is the correct one to utilize, this is still meant to show that the results are not due simply to model specification. 

spregress IYI c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

spregress IYI c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

spregress IYI c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

spregress IYI c.hdpdiff##c.educgengap HDP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

// For presenting the robustness check results in the same format as the main model:

eststo: quietly spregress IYI c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

eststo: quietly spregress IYI c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

eststo: quietly spregress IYI c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

eststo: quietly spregress IYI c.hdpdiff##c.educgengap HDP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

esttab using "table7.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("AKP CHP MHP HDP") title("IYI Party Vote Share in 2018 (Robustness Check)") 

// 3. The same robustness check as in (1) above is done for the HDP model:

// Use the dataset without the weight matrix.
use "abus_incongruity_dataset.dta", clear

reg HDP HDP15 c.educgengap##c.HDPrank AKPrank distmag mosquedens margin lnpop dideprat1 dibulge1 , robust

reg HDP HDP15 c.educgengap##c.HDPrank CHPrank distmag mosquedens margin lnpop dideprat1 dibulge1 , robust

reg HDP HDP15 c.educgengap##c.HDPrank MHPrank distmag mosquedens margin lnpop dideprat1 dibulge1 , robust 

reg HDP HDP15 c.educgengap##c.HDPrank IYIrank distmag mosquedens margin lnpop dideprat1 dibulge1 , robust

// For presenting the robustness check results in the same format as the main model:

eststo: quietly reg HDP HDP15 c.educgengap##c.HDPrank AKPrank mosquedens margin lnpop dideprat1 dibulge1 , robust

eststo: quietly reg HDP HDP15 c.educgengap##c.HDPrank CHPrank mosquedens margin lnpop dideprat1 dibulge1 , robust

eststo: quietly reg HDP HDP15 c.educgengap##c.HDPrank MHPrank mosquedens margin lnpop dideprat1 dibulge1 , robust

eststo: quietly reg HDP HDP15 c.educgengap##c.HDPrank IYIrank mosquedens margin lnpop dideprat1 dibulge1 , robust

esttab using "table8.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("AKP CHP MHP IYI") title("HDP Vote Share in 2018 (Robustness Check)") 

// 4. The same robustness check as in (2) above is done for the HDP model:

spregress HDP HDP15 c.educgengap##c.HDPrank AKPrank mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

spregress HDP HDP15 c.educgengap##c.HDPrank CHPrank mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

spregress HDP HDP15 c.educgengap##c.HDPrank MHPrank osquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

spregress HDP HDP15 c.educgengap##c.HDPrank IYIrank mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force 
estat impact

// For presenting the robustness check results in the same format as the main model:

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPrank AKPrank mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPrank CHPrank mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPrank MHPrank mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

eststo: quietly spregress HDP HDP15 c.educgengap##c.HDPrank IYIrank mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

esttab using "table9.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("AKP CHP MHP IYI") title("HDP Vote Share in 2018 (Robustness Check)") 

// 5. Robustness check using socio-economic development index instead of educational gender gap for IYI (spatial error models):

spregress IYI c.akpdiff##c.sedi AKP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress IYI c.chpdiff##c.sedi CHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress IYI c.mhpdiff##c.sedi MHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress IYI c.hdpdiff##c.sedi HDP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// For the table 
eststo: quietly spregress IYI c.akpdiff##c.sedi AKP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress IYI c.chpdiff##c.sedi CHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress IYI c.mhpdiff##c.sedi MHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress IYI c.hdpdiff##c.sedi HDP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

esttab using "table11.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("I II III IV") title("IYI Party Vote Share in 2018")

// 6. Robustness check using socio-economic development index instead of educational gender gap for HDP (spatial error models):

spregress HDP HDP15 c.akpdiff##c.sedi AKP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.chpdiff##c.sedi CHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.mhpdiff##c.sedi MHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// For the table 
eststo: quietly spregress HDP HDP15 c.akpdiff##c.sedi AKP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.chpdiff##c.sedi CHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.mhpdiff##c.sedi MHP15 mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

esttab using "table12.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("I II III") title("HDP Vote Share in 2018")

// 7. Using Educational Gender Gap as the main independent variable as in the model, but not modelling the spatial aspect of the data to demonstrate that the results are not driven by spatial modelling (HDP Model):

// Use the dataset without the weight matrix.
use "spatial_abus_incongruity_dataset.dta", clear

reg HDP HDP15 c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

reg HDP HDP15 c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

reg HDP HDP15 c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

// For presenting the robustness check results in the same format as the main model:

eststo: quietly reg HDP HDP15 c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

eststo: quietly reg HDP HDP15 c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

eststo: quietly reg HDP HDP15 c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , robust

esttab using "table13.tex", replace b(%9.3g) se r2 label nogap onecell ///
nonumbers mtitle("I II III") title("HDP Vote Share in 2018") 

// 8. Running the model with a spatial lag specification instead of a spatial error specification. Although theoretically the spatial error is the correct one to utilize, this is still meant to show that the results are not due simply to model specification (HDP models):

spregress HDP HDP15 c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

spregress HDP HDP15 c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

spregress HDP HDP15 c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force
estat impact

// For presenting the robustness check results in the same format as the main model:

eststo: quietly spregress HDP HDP15 c.akpdiff##c.educgengap AKP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

eststo: quietly spregress HDP HDP15 c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

eststo: quietly spregress HDP HDP15 c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

esttab using "table14.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("I II III") title("HDP Vote Share in 2018") 

// 9. The utilization of socio-economic development index in the additional implication derived from the theory in the paper:

spregress HDP HDP15 c.sedi##c.HDPrank AKPrank distmag mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.sedi##c.HDPrank CHPrank distmag mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.sedi##c.HDPrank MHPrank distmag mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

spregress HDP HDP15 c.sedi##c.HDPrank IYIrank distmag mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

// For the table
eststo: quietly spregress HDP HDP15 c.sedi##c.HDPrank AKPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.sedi##c.HDPrank CHPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.sedi##c.HDPrank MHPrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

eststo: quietly spregress HDP HDP15 c.sedi##c.HDPrank IYIrank mosquedens margin lnpop dideprat1 dibulge1 , gs2sls errorlag(W) force

esttab using "table15.tex", replace b(%9.3g) se pr2 label nogap onecell ///
nonumbers mtitle("I II III IV") title("HDP Vote Share in 2018")

// 10. Calculating the predicted vote share of IYI by using the vote gain of other parties (in the paper, magnitude of effects section shows the analysis using the vote loss of other parties). 

set graphics off

qui spregress IYI c.chpdiff##c.educgengap CHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

// Graph one: Mean educational gender gap -1 s.d. moderated by change in CHP Vote Share

qui margins, at(chpdiff=(0(1)6) educgengap=(0.0846255)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(2 4 6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(0 1 2 3 4 5 6, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(a, size(large) color(gs6)) ///
ylabel(, labsize(small) noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(2)12) title("") ytitle("") 

gr save marginplots/chp1robust.gph, replace

// Graph two: Mean educational gender gap moderated by change in CHP Vote Share

qui margins, at(chpdiff=(0(1)6) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(2 4 6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(0 1 2 3 4 5 6, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(b, size(large) color(gs6)) ///
ylabel(, labsize(small) noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(2)12) title("") ytitle("") 

gr save marginplots/chp2robust.gph, replace

// Graph three: Mean educational gender gap + 1 s.d. moderated by change in CHP Vote Share 

qui margins, at(chpdiff=(0(1)6) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(2 4 6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(0 1 2 3 4 5 6, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(c, size(large) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(2)12) title("") ytitle("") 

gr save marginplots/chp3robust.gph, replace
			

qui spregress IYI c.mhpdiff##c.educgengap MHP15 mosquedens margin lnpop dideprat1 dibulge1 , ml vce(robust) dvarlag(W) force

// Graph four: Mean educatinal gender gap -1 s.d. moderated by change in MHP Vote Share

qui margins, at(mhpdiff=(0(1)10) educgengap=(0.0846255)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(2 4 6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(0 2 4 6 8 10, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(d, size(large) color(gs6)) ///
ylabel(, labsize(small)noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(2)12) xlabel(0(2)10) title("") ytitle("") 

gr save marginplots/mhp1robust.gph, replace

// Graph five: Mean educational gender gap moderated by change in MHP Vote Share

qui margins, at(mhpdiff=(0(1)10) educgengap=(0.2244962)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(2 4 6 8 10 12, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(0 2 4 6 8 10, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(e, size(large) color(gs6)) ///
ylabel(, labsize(small) noticks) yscale(noline) xsize(4) ysize(4) legend(off) ylabel(0(2)12) xlabel(0(2)10) title("") ytitle("") 

gr save marginplots/mhp2robust.gph, replace

// Graph six: Mean educational gender gap + 1 s.d. moderated by change in MHP Vote Share 

qui margins, at(mhpdiff=(0(1)10) educgengap=(0.3644)) 

// At this point we can run marginsplot:

marginsplot, recast(line) recastci(rarea) plotopts(lcolor(black) /// 
lwidth(medthick) lpattern(solid)) ciopts(fcolor(gs11) lcolor(gs10)) ///
yline(0, lwidth(medthick) lpattern(solid) lcolor(white)) ///
yline(2 4 6 8 10, lwidth(medium) lpattern(solid) lcolor(white)) ///
xline(0 2 4 6 8 10, lwidth(medium) lpattern(solid) lcolor(white)) ///
graphregion(fcolor(white)) plotregion(fcolor(gs15)) xlabel(, labsize(small) ///
noticks) xscale(noline) xtitle(f, size(large) color(gs6)) ///
ylabel(, labsize(small) noticks) yscale(noline) xsize(4) ysize(4) legend(off)  ylabel(0(2)12) xlabel(0(2)10) title("") ytitle("") 

gr save marginplots/mhp3robust.gph, replace

// Turn the graphics back on for combined graph:

set graphics on			

// Now to combine the graphs into one image:

graph combine marginplots/chp1robust.gph marginplots/chp2robust.gph marginplots/chp3robust.gph ///
marginplots/mhp1robust.gph marginplots/mhp2robust.gph marginplots/mhp3robust.gph, ///
graphregion(fcolor(white)) l1(Predicted Vote Share of IYI Party (2018), size(medium)) b1(Change in Vote Shares of CHP(a-c) and MHP(d-f) (2015-2018), size(medium)) title(Mean Educational Gender Gap ± 1 s.d., size(medium)) 

gr save marginplots/marginsiyirobust.gph, replace 

gr export marginplots/marginsiyirobust.pdf, replace
gr export marginplots/marginsiyirobust.svg, replace

*************************END OF ROBUSTNESS CHECKS*****************************	

/* Appendix E Chloropleth Maps for Candidate Placement

The code to reproduce the maps offered in Appendix E of the paper is presented below: Maps are produced to show the spread of the number of female candidates, highest placement of femela candidates, and the number of female candidates in winnable slots. These maps are presented for the reader to visualize the basis of the main argument and theory presented in the paper. In order to provide a visual contrast between IYI and HDP on one hand and the other parliamentary parties on the other hand, these variables are mapped for all parliemanetary parties. Be sure to have spset the data before usıng the grmap command to reprıduce the maps. Because each variable will be mapped five times, IYI and HDP is combined into one page and other parties are grouped into one page for all variables to control the length of the appendix.       */   
// We will make a new folder for these appendix graphs 

mkdir appendixmaps

set graphics off 

//Figure 9 (Mapping Number of Female Candidates for IYI)
grmap IYIfem, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure9.gph, replace
graph export appendixmaps/figure9.pdf, replace
graph export appendixmaps/figure9.svg, replace

grmap IYIfem, clmethod(unique) fcolor(Reds) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure9color.gph, replace
graph export appendixmaps/figure9color.pdf, replace
graph export appendixmaps/figure9color.svg, replace

// Figure 10 (Mapping Number of Female Candidates for HDP)
grmap HDPfem, clmethod(unique) fcolor(Greys2) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure10.gph, replace
graph export appendixmaps/figure10.pdf, replace
graph export appendixmaps/figure10.svg, replace

grmap HDPfem, clmethod(unique) fcolor(Reds2) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure10color.gph, replace
graph export appendixmaps/figure10color.pdf, replace
graph export appendixmaps/figure10color.svg, replace

// Combining the 2 graphs into one image:

set graphics on 

graph combine appendixmaps/figure9.gph appendixmaps/figure10.gph, rows(2) altshrink

gr save appendixmaps/combined1.gph, replace 
gr export appendixmaps/combined1.pdf, replace
gr export appendixmaps/combined1.svg, replace

graph combine appendixmaps/figure9color.gph appendixmaps/figure10color.gph, rows(2) altshrink

gr save appendixmaps/combinedcolor1.gph, replace 
gr export appendixmaps/combinedcolor1.pdf, replace
gr export appendixmaps/combinedcolor1.svg, replace

set graphics off

// Figure 11 (Mapping Number of Female Candidates for AKP)
grmap AKPfem, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure11.gph, replace
graph export appendixmaps/figure11.pdf, replace
graph export appendixmaps/figure11.svg, replace

grmap AKPfem, clmethod(unique) fcolor(Reds) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure11color.gph, replace
graph export appendixmaps/figure11color.pdf, replace
graph export appendixmaps/figure11color.svg, replace

// Figure 12 (Mapping Number of Female Candidates for CHP)
grmap CHPfem, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure12.gph, replace
graph export appendixmaps/figure12.pdf, replace
graph export appendixmaps/figure12.svg, replace

grmap CHPfem, clmethod(unique) fcolor(Reds) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure12color.gph, replace
graph export appendixmaps/figure12color.pdf, replace
graph export appendixmaps/figure12color.svg, replace

// Figure 13 (Mapping Number of Female Candidates for MHP)
grmap MHPfem, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure13.gph, replace
graph export appendixmaps/figure13.pdf, replace
graph export appendixmaps/figure13.svg, replace

grmap MHPfem, clmethod(unique) fcolor(Reds) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number of Candidates) 
graph save appendixmaps/figure13color.gph, replace
graph export appendixmaps/figure13color.pdf, replace
graph export appendixmaps/figure13color.svg, replace

// Combining the 3 graphs into one image:

set graphics on

graph combine appendixmaps/figure11.gph appendixmaps/figure12.gph appendixmaps/figure13.gph, rows(3) altshrink

gr save appendixmaps/combined2.gph, replace 
gr export appendixmaps/combined2.pdf, replace
gr export appendixmaps/combined2.svg, replace

graph combine appendixmaps/figure11color.gph appendixmaps/figure12color.gph appendixmaps/figure13color.gph, rows(3) altshrink

gr save appendixmaps/combinedcolor2.gph, replace 
gr export appendixmaps/combinedcolor2.pdf, replace
gr export appendixmaps/combinedcolor2.svg, replace

set graphics off

// Figure 14 (Mapping Highest Position of Female Candidates fpr IYI)
grmap IYIrank, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure14.gph, replace
graph export appendixmaps/figure14.pdf, replace
graph export appendixmaps/figure14.svg, replace

grmap IYIrank, clmethod(unique) fcolor(Blues) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure14color.gph, replace
graph export appendixmaps/figure14color.pdf, replace
graph export appendixmaps/figure14color.svg, replace

// Figure 15 (Mapping Highest Position of Female Candidates for HDP)
grmap HDPrank, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure15.gph, replace
graph export appendixmaps/figure15.pdf, replace
graph export appendixmaps/figure15.svg, replace

grmap HDPrank, clmethod(unique) fcolor(Blues) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure15color.gph, replace
graph export appendixmaps/figure15color.pdf, replace
graph export appendixmaps/figure15color.svg, replace

// Combining the 2 graphs into one image:

set graphics on

graph combine appendixmaps/figure14.gph appendixmaps/figure15.gph, rows(2) altshrink

gr save appendixmaps/combined3.gph, replace 
gr export appendixmaps/combined3.pdf, replace
gr export appendixmaps/combined3.svg, replace

graph combine appendixmaps/figure14color.gph appendixmaps/figure15color.gph, rows(2) altshrink

gr save appendixmaps/combinedcolor3.gph, replace 
gr export appendixmaps/combinedcolor3.pdf, replace
gr export appendixmaps/combinedcolor3.svg, replace

set graphics off

// Figure 16 (Mapping Highest Position of Female Candidates for AKP),
grmap AKPrank, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure16.gph, replace
graph export appendixmaps/figure16.pdf, replace
graph export appendixmaps/figure16.svg, replace

grmap AKPrank, clmethod(unique) fcolor(Blues) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure16color.gph, replace
graph export appendixmaps/figure16color.pdf, replace
graph export appendixmaps/figure16color.svg, replace

// Figure 17 (Mapping Highest Position of Female Candidates for CHP)
grmap CHPrank, clmethod(unique) fcolor(Greys2) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure17.gph, replace
graph export appendixmaps/figure17.pdf, replace
graph export appendixmaps/figure17.svg, replace

grmap CHPrank, clmethod(unique) fcolor(Blues2) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure17color.gph, replace
graph export appendixmaps/figure17color.pdf, replace
graph export appendixmaps/figure17color.svg, replace

// Figure 18 (Mapping Highest Position of Female Candidates for MHP)
grmap MHPrank, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure18.gph, replace
graph export appendixmaps/figure18.pdf, replace
graph export appendixmaps/figure18.svg, replace

grmap MHPrank, clmethod(unique) fcolor(Blues) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Ranking of Candidates) 
graph save appendixmaps/figure18color.gph, replace
graph export appendixmaps/figure18color.pdf, replace
graph export appendixmaps/figure18color.svg, replace

// Combining the 3 graphs into one image:

set graphics on

graph combine appendixmaps/figure16.gph appendixmaps/figure17.gph appendixmaps/figure18.gph, rows(3) altshrink

gr save appendixmaps/combined4.gph, replace 
gr export appendixmaps/combined4.pdf, replace
gr export appendixmaps/combined4.svg, replace

graph combine appendixmaps/figure16color.gph appendixmaps/figure17color.gph appendixmaps/figure18color.gph, rows(3) altshrink

gr save appendixmaps/combinedcolor4.gph, replace 
gr export appendixmaps/combinedcolor4.pdf, replace
gr export appendixmaps/combinedcolor4.svg, replace

set graphics off

// Figure 19 (Mapping Number of Winnable Slots for Women - IYI)
grmap IYIwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save appendixmaps/figure19.gph, replace
graph export appendixmaps/figure19.pdf, replace
graph export appendixmaps/figure19.svg, replace

grmap IYIwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save appendixmaps/figure19color.gph, replace
graph export appendixmaps/figure19color.pdf, replace
graph export appendixmaps/figure19color.svg, replace

// Figure 20 (Mapping Number of Winnable Slots for Women - HDP)
grmap HDPwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnig Slots) 
graph save appendixmaps/figure20.gph, replace
graph export appendixmaps/figure20.pdf, replace
graph export appendixmaps/figure20.svg, replace

grmap HDPwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save appendixmaps/figure20color.gph, replace
graph export appendixmaps/figure20color.pdf, replace
graph export appendixmaps/figure20color.svg, replace

// Combining the 2 graphs into one image:

set graphics on

graph combine appendixmaps/figure19.gph appendixmaps/figure20.gph, rows(2) altshrink

gr save appendixmaps/combined5.gph, replace 
gr export appendixmaps/combined5.pdf, replace
gr export appendixmaps/combined5.svg, replace

graph combine appendixmaps/figure19color.gph appendixmaps/figure20color.gph, rows(2) altshrink

gr save appendixmaps/combinedcolor5.gph, replace 
gr export appendixmaps/combinedcolor5.pdf, replace
gr export appendixmaps/combinedcolor5.svg, replace

set graphics off 

// Figure 21 (Mapping Number of Winnable Slots for Women - AKP)
grmap AKPwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save appendixmaps/figure21.gph, replace
graph export appendixmaps/figure21.pdf, replace
graph export appendixmaps/figure21.svg, replace

grmap AKPwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save appendixmaps/figure21color.gph, replace
graph export appendixmaps/figure21color.pdf, replace
graph export appendixmaps/figure21color.svg, replace

// Figure 22 (Mapping Number of Winnable Slots for Women - CHP)
grmap CHPwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save appendixmaps/figure22.gph, replace
graph export appendixmaps/figure22.pdf, replace
graph export appendixmaps/figure22.svg, replace

grmap CHPwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save appendixmaps/figure22color.gph, replace
graph export appendixmaps/figure22color.pdf, replace
graph export appendixmaps/figure22color.svg, replace

// Figure 23 (Mapping Number of Winnable Slots for Women - MHP)
grmap MHPwinslot, clmethod(unique) fcolor(Greys) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save appendixmaps/figure23.gph, replace
graph export appendixmaps/figure23.pdf, replace
graph export appendixmaps/figure23.svg, replace

grmap MHPwinslot, clmethod(unique) fcolor(Greens) osize(vvthin) ndfcolor(black) legenda(on) legstyle(1) legtitle(Number in Winnable Slots) 
graph save appendixmaps/figure23color.gph, replace
graph export appendixmaps/figure23color.pdf, replace
graph export appendixmaps/figure23color.svg, replace

// Combining the 3 graphs into one image:

set graphics on

graph combine appendixmaps/figure21.gph appendixmaps/figure22.gph appendixmaps/figure23.gph, rows(3) altshrink

gr save appendixmaps/combined6.gph, replace 
gr export appendixmaps/combined6.pdf, replace
gr export appendixmaps/combined6.svg, replace

graph combine appendixmaps/figure21color.gph appendixmaps/figure22color.gph appendixmaps/figure23color.gph, rows(3) altshrink

gr save appendixmaps/combinedcolor6.gph, replace 
gr export appendixmaps/combinedcolor6.pdf, replace
gr export appendixmaps/combinedcolor6.svg, replace

*************************END OF APPENDIX E MAP CODE*************************

// APPENDIX G : Gender Gap In Global Perspective 
// The following code reproduces the time-series plots found in Appendix G.

use "gendergapdata.dta", clear

describe

// The dataset contains the annual values on gender gap index and two subindices for Turkey from 2006-2018. Inorder to plot, we will define the dataset as a time series. 

mkdir appendixgraphs

tsset year

tsline gendergapind gendergapave, xtitle("")

gr save appendixgraphs/gendergap1.gph, replace 
gr export appendixgraphs/gendergap1.pdf, replace
gr export appendixgraphs/gendergap1.svg, replace

tsline econ econave, xtitle("")

gr save appendixgraphs/gendergap2.gph, replace 
gr export appendixgraphs/gendergap2.pdf, replace
gr export appendixgraphs/gendergap2.svg, replace

tsline politemp politempave, xtitle("")

gr save appendixgraphs/gendergap3.gph, replace 
gr export appendixgraphs/gendergap3.pdf, replace
gr export appendixgraphs/gendergap3.svg, replace

*************************END OF APPENDIX CODE******************************

******************************END OF DO-FILE********************************
