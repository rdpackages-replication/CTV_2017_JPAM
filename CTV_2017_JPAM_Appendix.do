********************************************************************************
** Comparing Inference Approaches for RD Designs:
** A Reexamination of the Effect of Head Start on Child Mortality
** SUPPLEMENTAL APPENDIX
** Authors: Matias D. Cattaneo, Rocio Titiunik and Gonzalo Vazquez-Bare
** Last update: 23-FEB-2021
********************************************************************************
** SOFTWARE WEBSITE: https://rdpackages.github.io/
********************************************************************************
** TO INSTALL STATA PACKAGES:
** RDROBUST: net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace
** RDDENSITY: net install rdlocrand, from(https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata) replace
** RDLOCRAND: net install rddensity, from(https://raw.githubusercontent.com/rdpackages/rddensity/master/stata) replace
********************************************************************************
** NOTE: if you are using RDROBUST version 2020 or newer, the option 
** "masspoints(off) stdvars(on)" may be needed to replicate the results in
** the paper.
********************************************************************************
** NOTE: if you are using RDDENSITY version 2020 or newer, the option 
** "nomasspoints" may be needed to replicate the results in the paper.
********************************************************************************

********************************************************************************
** Load and Setup Data
********************************************************************************
use headstart, clear

gl Y mort_age59_related_postHS
gl R povrate60

drop if $R==.

foreach var of varlist census1960_pop* {
	qui replace `var'=`var'/1000 		// express pop variables in thousands
}

gl covs60 "census1960_pop census1960_pctsch1417 census1960_pctsch534 census1960_pctsch25plus census1960_pop1417 census1960_pop534 census1960_pop25plus census1960_pcturban census1960_pctblack" 
gl covs90 "census1990_pop census1990_pop1824 census1990_pop2534 census1990_pop3554 census1990_pop55plus census1990_pcturban census1990_pctblack census1990_percapinc"
global placebos "mort_age59_injury_postHS mort_age59_all_postHS mort_age25plus_related_postHS mort_age25plus_injuries_postHS mort_age59_related_preHS mort_wh_age59_related_postHS mort_bl_age59_related_postHS"

qui {
	gen double R = $R - 59.1984
	gen double D = $R >= 59.1984
	gen double R1 = ($R-59.1984)
	gen double R2 = ($R-59.1984)^2
	gen double R3 = ($R-59.1984)^3
	gen double R4 = ($R-59.1984)^4
	gen double DR1 = D*R1 
	gen double DR2 = D*R2
	gen double DR3 = D*R3
	gen double DR4 = D*R4
}

********************************************************************************
** Table SA-1: Outcome and Score Descriptive Stats
********************************************************************************

mat T = J(10,2,.)
local count = 1
foreach var in $R $Y{
	qui sum `var', det
	mat T[1,`count'] = round(r(mean),.001)
	mat T[2,`count'] = round(r(sd),.001)
	mat T[3,`count'] = round(r(min),.001)
	mat T[4,`count'] = round(r(p10),.001)
	mat T[5,`count'] = round(r(p25),.001)
	mat T[6,`count'] = round(r(p50),.001)
	mat T[7,`count'] = round(r(p75),.001)
	mat T[8,`count'] = round(r(p90),.001)
	mat T[9,`count'] = round(r(max),.001)
	mat T[10,`count'] = r(N)
	local ++count
}

matlist T

********************************************************************************
** Figure SA-1 and SA-2: Histograms
********************************************************************************

** Histograms with full sample

hist $R, bin(60) addplot(pci 0 59.1984 .032 59.1984, lpattern(shortdash) lcolor(black)) ///
	legend(off) graphregion(color(white))

hist $Y, bin(60) graphregion(color(white))

hist $Y if $Y>0, bin(60) graphregion(color(white))

** Histograms in window

local w = 1
hist $R if abs(R)<=`w', bin(20) addplot(pci 0 59.1984 1.6 59.1984, lpattern(shortdash) lcolor(black)) ///
	legend(off) graphregion(color(white))

hist $Y if abs(R)<=`w', bin(20) graphregion(color(white))

hist $Y if $Y>0 & abs(R)<=`w', bin(12) graphregion(color(white))

local w = 3
hist $R if abs(R)<=`w', bin(50) addplot(pci 0 59.1984 .4 59.1984, lpattern(shortdash) lcolor(black)) ///
	legend(off) graphregion(color(white))

hist $Y if abs(R)<=`w', bin(60) graphregion(color(white))

hist $Y if $Y>0 & abs(R)<=`w', bin(60) graphregion(color(white))

local w = 9
hist $R if abs(R)<=`w', bin(50) addplot(pci 0 59.1984 .1 59.1984, lpattern(shortdash) lcolor(black)) ///
	legend(off) graphregion(color(white))

hist $Y if abs(R)<=`w', bin(60) graphregion(color(white))

hist $Y if $Y>0 & abs(R)<=`w', bin(60) graphregion(color(white))

********************************************************************************
** Figure SA-3: Rdplots
********************************************************************************

local nombre = 1
foreach var of varlist $Y $placebos {
	rdplot `var' $R if $R<90, c(59.1984) graph_options(graphregion(color(white)) ///
		ytitle("`lab'") title(""))
	local ++nombre
}

********************************************************************************
** Tables SA-2 to SA-5: Summary Statistics and Difference in Means
********************************************************************************

mat T1 = J(18,9,.)
mat T2 = J(18,9,.)
mat T3 = J(18,9,.)
mat T4 = J(18,9,.)

local tabla = 1
foreach w of numlist 100 9 3 1{
	local fila = 1
	foreach var of varlist $Y $R $placebos $covs60 {
		qui ttest `var' if abs(R)<=`w', by(D) unequal
		mat T`tabla'[`fila',1] = r(N_1)
		mat T`tabla'[`fila',2] = round(r(mu_1),.001)
		mat T`tabla'[`fila',3] = round(r(sd_1)/sqrt(r(N_1)),.001)
		mat T`tabla'[`fila',4] = r(N_2)
		mat T`tabla'[`fila',5] = round(r(mu_2),.001)
		mat T`tabla'[`fila',6] = round(r(sd_2)/sqrt(r(N_2)),.001)
		mat T`tabla'[`fila',7] = round(r(mu_2)-r(mu_1),.001)
		mat T`tabla'[`fila',8] = round(r(se),.001)
		mat T`tabla'[`fila',9] = round(r(p),.001)
		
		local ++fila
	}
	matlist T`tabla'

	local ++tabla
}

********************************************************************************
** Figure SA-4: Falsification Test - Rdplots
********************************************************************************

rdplot $Y $R, c(59.1984) binselect(es) graph_options(graphregion(color(white)) ///
	ytitle("`lab'") title(""))

rdplot $Y $R, c(59.1984) binselect(qs) graph_options(graphregion(color(white)) ///
	ytitle("`lab'") title(""))

rdplot $Y $R, c(59.1984) binselect(qsmv) graph_options(graphregion(color(white)) title(""))

********************************************************************************
** Tables SA-6 to SA-8: Local Polynomial, Main and Placebo Outcomes
********************************************************************************

** Constant Fit (p=0)
mat T30 = J(48,4,.)
** Linear Fit (p=1)
mat T31 = J(48,4,.)

foreach tabla of numlist 0 1 {
	local fila = 0
	foreach var of varlist $Y $placebos {
		qui rdrobust `var' R, p(`tabla') bwselect(cerrd)
		mat aux = e(b)
		mat T3`tabla'[1+`fila',1] = e(p)
		mat T3`tabla'[2+`fila',1] = round(e(h_l),.001)
		mat T3`tabla'[3+`fila',1] = round(aux[1,1],.001)
		mat T3`tabla'[4+`fila',1] = round(e(pv_rb),.001)
		mat T3`tabla'[5+`fila',1] = e(N_h_l)
		mat T3`tabla'[6+`fila',1] = e(N_h_r)

		qui rdrobust `var' R, p(`tabla') bwselect(mserd)
		mat aux = e(b)
		mat T3`tabla'[1+`fila',2] = e(p)
		mat T3`tabla'[2+`fila',2] = round(e(h_l),.001)
		mat T3`tabla'[3+`fila',2] = round(aux[1,1],.001)
		mat T3`tabla'[4+`fila',2] = round(e(pv_rb),.001)
		mat T3`tabla'[5+`fila',2] = e(N_h_l)
		mat T3`tabla'[6+`fila',2] = e(N_h_r)

		qui rdrobust `var' R, p(`tabla') h(9)
		mat aux = e(b)
		mat T3`tabla'[1+`fila',3] = e(p)
		mat T3`tabla'[2+`fila',3] = round(e(h_l),.001)
		mat T3`tabla'[3+`fila',3] = round(aux[1,1],.001)
		mat T3`tabla'[4+`fila',3] = round(e(pv_rb),.001)
		mat T3`tabla'[5+`fila',3] = e(N_h_l)
		mat T3`tabla'[6+`fila',3] = e(N_h_r)

		qui rdrobust `var' R, p(`tabla') h(18)
		mat aux = e(b)
		mat T3`tabla'[1+`fila',4] = e(p)
		mat T3`tabla'[2+`fila',4] = round(e(h_l),.001)
		mat T3`tabla'[3+`fila',4] = round(aux[1,1],.001)
		mat T3`tabla'[4+`fila',4] = round(e(pv_rb),.001)
		mat T3`tabla'[5+`fila',4] = e(N_h_l)
		mat T3`tabla'[6+`fila',4] = e(N_h_r)
		
		local fila = `fila'+6
	}

	matlist T3`tabla'

}

** Constant Fit (p=0)
mat T40 = J(13,4,.)
** Linear Fit (p=1)
mat T41 = J(13,4,.)
** Quartic Fit (p=4)
mat T44 = J(13,4,.)


foreach tabla of numlist 0 1 4{
	qui rdrobust $Y R, p(`tabla') bwselect(cerrd)
	mat aux = e(b)
    mat T4`tabla'[1,1] = e(p)
	mat T4`tabla'[2,1] = round(e(h_l),.001)
	mat T4`tabla'[3,1] = round(aux[1,1],.001)
	mat T4`tabla'[4,1] = round(e(pv_rb),.001)
	mat T4`tabla'[5,1] = e(N_h_l)
	mat T4`tabla'[6,1] = e(N_h_r)

	qui rdrobust $Y R, p(`tabla') bwselect(mserd)
	mat aux = e(b)
    mat T4`tabla'[1,2] = e(p)
	mat T4`tabla'[2,2] = round(e(h_l),.001)
	mat T4`tabla'[3,2] = round(aux[1,1],.001)
	mat T4`tabla'[4,2] = round(e(pv_rb),.001)
	mat T4`tabla'[5,2] = e(N_h_l)
	mat T4`tabla'[6,2] = e(N_h_r)

	qui rdrobust $Y R, p(`tabla') h(9)
	mat aux = e(b)
    mat T4`tabla'[1,3] = e(p)
	mat T4`tabla'[2,3] = round(e(h_l),.001)
	mat T4`tabla'[3,3] = round(aux[1,1],.001)
	mat T4`tabla'[4,3] = round(e(pv_rb),.001)
	mat T4`tabla'[5,3] = e(N_h_l)
	mat T4`tabla'[6,3] = e(N_h_r)

	qui rdrobust $Y R, p(`tabla') h(18)
	mat aux = e(b)
    mat T4`tabla'[1,4] = e(p)
	mat T4`tabla'[2,4] = round(e(h_l),.001)
	mat T4`tabla'[3,4] = round(aux[1,1],.001)
	mat T4`tabla'[4,4] = round(e(pv_rb),.001)
	mat T4`tabla'[5,4] = e(N_h_l)
	mat T4`tabla'[6,4] = e(N_h_r)

	local fila = 7
	foreach y of varlist mort_age59_injury_postHS mort_age59_all_postHS mort_age25plus_related_postHS ///
						 mort_age25plus_injuries_postHS mort_age59_related_preHS mort_wh_age59_related_postHS mort_bl_age59_related_postHS {
	   qui rdrobust `y' R, p(`tabla') bwselect(cerrd)
	   mat T4`tabla'[`fila',1] = round(e(pv_rb),.001)
	   qui rdrobust `y' R, p(`tabla') bwselect(mserd)
	   mat T4`tabla'[`fila',2] = round(e(pv_rb),.001)
	   qui rdrobust `y' R, p(`tabla') h(9)
	   mat T4`tabla'[`fila',3] = round(e(pv_rb),.001)
	   qui rdrobust `y' R, p(`tabla') h(18)
	   mat T4`tabla'[`fila',4] = round(e(pv_rb),.001)
	   local ++fila
	}
}

matlist T40
matlist T41
matlist T44

********************************************************************************
** Figure SA-5: Local Randomization Methods - Window Selection
********************************************************************************

gl wreps = 1000

** WINDOW SELECTION
** NOTE: the plots are drawn using the asymptotic p-value to speed up the process.
** Remove the "approx" option to use randinf and replicate the results in the paper.

rdwinselect R mort_age59_related_preHS $covs60, reps($wreps) stat(ksmirnov) wmin(.3) wstep(.2) nwin(40) level(.2) approx
mat Res = r(results)
preserve
svmat Res
rename Res1 pvalues 
rename Res6 w
gen red = pval if Res3==43
twoway(scatter pval w)(scatter red w, msize(vlarge) msymbol(circle_hollow) mlwidth(medthick)), ///
	xline(1.1,lpattern(shortdash)) ytitle(p-values) xtitle(bandwidth) ///
	xlabel(0.3(.4)8.1, labsize(small)) legend(off) graphregion(color(white))
restore

rdwinselect R mort_age59_related_preHS $covs60, reps($wreps) stat(ttest) wmin(.3) wstep(.2) nwin(40) level(.2) approx
mat Res = r(results)
preserve
svmat Res
rename Res1 pvalues 
rename Res6 w
gen red = pval if Res3==53
twoway(scatter pval w)(scatter red w, msize(vlarge) msymbol(circle_hollow) mlwidth(medthick)), ///
	xline(1.5,lpattern(shortdash)) ytitle(p-values) xtitle(bandwidth) ///
	xlabel(0.3(.4)8.1, labsize(small)) legend(off) graphregion(color(white))
restore

rdwinselect R mort_age59_related_preHS $covs60, reps($wreps) stat(ranksum) wmin(.3) wstep(.2) nwin(40) level(.2) approx
mat Res = r(results)
preserve
svmat Res
rename Res1 pvalues 
rename Res6 w
gen red = pval if Res3==51
twoway(scatter pval w)(scatter red w, msize(vlarge) msymbol(circle_hollow) mlwidth(medthick)), ///
	xline(1.3,lpattern(shortdash)) ytitle(p-values) xtitle(bandwidth) ///
	xlabel(0.3(.4)8.1, labsize(small)) legend(off) graphregion(color(white))
restore

rdwinselect R mort_age59_related_preHS $covs60, reps($wreps) stat(hotelling) wmin(.3) wstep(.2) nwin(40) level(.2) approx
mat Res = r(results)
preserve
svmat Res
rename Res1 pvalues 
rename Res6 w
gen red = pval if Res3==81
twoway(scatter pval w)(scatter red w, msize(vlarge) msymbol(circle_hollow) mlwidth(medthick)), ///
	xline(2.7,lpattern(shortdash)) ytitle(p-values) xtitle(bandwidth) ///
	xlabel(0.3(.4)8.1, labsize(small)) legend(off) graphregion(color(white))
restore

********************************************************************************
** Tables SA-9 and SA-10: Local Randomization Methods
********************************************************************************

gl rreps = 1000

** No Adjustment (p=0)

mat T50 = J(13,5,.)
mat T60 = J(13,5,.)

** Linear Adjustment (p=1)

mat T51 = J(13,5,.)
mat T61 = J(13,5,.)

foreach tabla of numlist 0 1 {
	qui rdrandinf $Y R, wl(-0.9) wr(0.9) reps($rreps) p(`tabla')
	mat T5`tabla'[1,1] = r(p)
	mat T5`tabla'[2,1] = round(r(wr),.001)
	mat T5`tabla'[3,1] = round(r(obs_stat),.001)
	mat T5`tabla'[4,1] = round(r(randpval),.001)
	mat T5`tabla'[5,1] = r(N_left)
	mat T5`tabla'[6,1] = r(N_right)
	
	mat T6`tabla'[1,1] = r(p)
	mat T6`tabla'[2,1] = round(r(wr),.001)
	mat T6`tabla'[3,1] = round(r(obs_stat),.001)
	mat T6`tabla'[4,1] = round(r(asy_pval),.001)
	mat T6`tabla'[5,1] = r(N_left)
	mat T6`tabla'[6,1] = r(N_right)

	qui rdrandinf $Y R, wl(-1.1) wr(1.1) reps($rreps) p(`tabla')
	mat T5`tabla'[1,2] = r(p)
	mat T5`tabla'[2,2] = round(r(wr),.001)
	mat T5`tabla'[3,2] = round(r(obs_stat),.001)
	mat T5`tabla'[4,2] = round(r(randpval),.001)
	mat T5`tabla'[5,2] = r(N_left)
	mat T5`tabla'[6,2] = r(N_right)
	
	mat T6`tabla'[1,2] = r(p)
	mat T6`tabla'[2,2] = round(r(wr),.001)
	mat T6`tabla'[3,2] = round(r(obs_stat),.001)
	mat T6`tabla'[4,2] = round(r(asy_pval),.001)
	mat T6`tabla'[5,2] = r(N_left)
	mat T6`tabla'[6,2] = r(N_right)

	qui rdrandinf $Y R, wl(-1.3) wr(1.3) reps($rreps) p(`tabla')
	mat T5`tabla'[1,3] = r(p)
	mat T5`tabla'[2,3] = round(r(wr),.001)
	mat T5`tabla'[3,3] = round(r(obs_stat),.001)
	mat T5`tabla'[4,3] = round(r(randpval),.001)
	mat T5`tabla'[5,3] = r(N_left)
	mat T5`tabla'[6,3] = r(N_right)
	
	mat T6`tabla'[1,3] = r(p)
	mat T6`tabla'[2,3] = round(r(wr),.001)
	mat T6`tabla'[3,3] = round(r(obs_stat),.001)
	mat T6`tabla'[4,3] = round(r(asy_pval),.001)
	mat T6`tabla'[5,3] = r(N_left)
	mat T6`tabla'[6,3] = r(N_right)

	qui rdrandinf $Y R, wl(-1.5) wr(1.5) reps($rreps) p(`tabla')
	mat T5`tabla'[1,4] = r(p)
	mat T5`tabla'[2,4] = round(r(wr),.001)
	mat T5`tabla'[3,4] = round(r(obs_stat),.001)
	mat T5`tabla'[4,4] = round(r(randpval),.001)
	mat T5`tabla'[5,4] = r(N_left)
	mat T5`tabla'[6,4] = r(N_right)
	
	mat T6`tabla'[1,4] = r(p)
	mat T6`tabla'[2,4] = round(r(wr),.001)
	mat T6`tabla'[3,4] = round(r(obs_stat),.001)
	mat T6`tabla'[4,4] = round(r(asy_pval),.001)
	mat T6`tabla'[5,4] = r(N_left)
	mat T6`tabla'[6,4] = r(N_right)

	qui rdrandinf $Y R, wl(-2.7) wr(2.7) reps($rreps) p(`tabla')
	mat T5`tabla'[1,5] = r(p)
	mat T5`tabla'[2,5] = round(r(wr),.001)
	mat T5`tabla'[3,5] = round(r(obs_stat),.001)
	mat T5`tabla'[4,5] = round(r(randpval),.001)
	mat T5`tabla'[5,5] = r(N_left)
	mat T5`tabla'[6,5] = r(N_right)

	mat T6`tabla'[1,5] = r(p)
	mat T6`tabla'[2,5] = round(r(wr),.001)
	mat T6`tabla'[3,5] = round(r(obs_stat),.001)
	mat T6`tabla'[4,5] = round(r(asy_pval),.001)
	mat T6`tabla'[5,5] = r(N_left)
	mat T6`tabla'[6,5] = r(N_right)
	
	local fila = 7
	foreach y of varlist mort_age59_injury_postHS mort_age59_all_postHS mort_age25plus_related_postHS ///
						 mort_age25plus_injuries_postHS mort_age59_related_preHS mort_wh_age59_related_postHS mort_bl_age59_related_postHS {
	   qui rdrandinf `y' R, wl(-0.9) wr(0.9) reps($rreps) p(`tabla')
	   mat T5`tabla'[`fila',1] = round(r(randpval),.001)
	   mat T6`tabla'[`fila',1] = round(r(asy_pval),.001)
	   qui rdrandinf `y' R, wl(-1.1) wr(1.1) reps($rreps) p(`tabla')
	   mat T5`tabla'[`fila',2] = round(r(randpval),.001)
	   mat T6`tabla'[`fila',2] = round(r(asy_pval),.001)
	   qui rdrandinf `y' R, wl(-1.3) wr(1.3) reps($rreps) p(`tabla')
	   mat T5`tabla'[`fila',3] = round(r(randpval),.001)
	   mat T6`tabla'[`fila',3] = round(r(asy_pval),.001)
	   qui rdrandinf `y' R, wl(-1.5) wr(1.5) reps($rreps) p(`tabla')
	   mat T5`tabla'[`fila',4] = round(r(randpval),.001)
	   mat T6`tabla'[`fila',4] = round(r(asy_pval),.001)
	   qui rdrandinf `y' R, wl(-2.7) wr(2.7) reps($rreps) p(`tabla')
	   mat T5`tabla'[`fila',5] = round(r(randpval),.001)
	   mat T6`tabla'[`fila',5] = round(r(asy_pval),.001)
	   local ++fila
	}
}

matlist T50
matlist T60
matlist T51
matlist T61

********************************************************************************
** Figure SA-6: Sensitivity of Linear Adjustment Model
********************************************************************************

twoway (scatter $Y R if abs(R)<=1.1)(lfit $Y R if abs(R)<=1.1 & D==1, lcolor(black) lpattern(shortdash)) ///
	(lfit $Y R if abs(R)<=1.1 & D==0, lcolor(black) lpattern(shortdash)), ///
	graphregion(color(white)) title("") legend(off) ylabel(0(20)80) xlabel(-2(1)2) ///
	ytitle("`lab'")

twoway (scatter $Y R if abs(R)<=1.3)(lfit $Y R if abs(R)<=1.3 & D==1, lcolor(black) lpattern(shortdash)) ///
	(lfit $Y R if abs(R)<=1.3 & D==0, lcolor(black) lpattern(shortdash)), ///
	graphregion(color(white)) title("") legend(off) ylabel(0(20)80) xlabel(-2(1)2) ///
	ytitle("`lab'")

********************************************************************************
** Table SA-11: Comparison of Inference Approaches
********************************************************************************

** NOTE: this table is generated at the end of the main paper replication do file.
** See Figure 4: Summary of Results


