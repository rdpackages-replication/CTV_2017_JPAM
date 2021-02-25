********************************************************************************
** Comparing Inference Approaches for RD Designs:
** A Reexamination of the Effect of Head Start on Child Mortality
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
** NOTE: If you are using RDROBUST version 2020 or newer, the option 
** "masspoints(off) stdvars(on)" may be needed to replicate the results in 
** the paper. For example, line 210:
**
**     rdrobust $Y R, p(0)
**
** should be replaced by:
**
**     rdrobust $Y R, p(0) masspoints(off) stdvars(on)
********************************************************************************
** NOTE: If you are using RDDENSITY version 2020 or newer, the option 
** "nomasspoints" may be needed to replicate the results in the paper.
** For example, line 77:
**
**     rddensity R
**
** should be replaced by:
**
**     rddensity R, nomasspoints
********************************************************************************

********************************************************************************
** Load Data
********************************************************************************

use headstart, clear
gl Y mort_age59_related_postHS
gl R povrate60
gl c = 59.1984

gen double R = $R - $c

********************************************************************************
** Figure 1: Scatter and RD plot
********************************************************************************

rdplot $Y $R if $Y <= 20, c($c) nbins(3000)	///
	graph_options(graphregion(color(white)) title("") /// 
	ytitle("`lab'") legend(off))

rdplot $Y $R, c($c) graph_options(graphregion(color(white)) ///
	title("") ytitle("`lab'") legend(off))

********************************************************************************
** Table 1: Binomial Tests
********************************************************************************

rdwinselect R if $Y!=., wmin(.3) wstep(.2) nwin(6)
mat tmp = r(results)
mat T = (tmp[1..6,5],tmp[1..6,3..4],tmp[1..6,2])
forvalues row=1/6{
	mat T[`row',1] = round(T[`row',1],.001)
	mat T[`row',4] = round(T[`row',4],.001)
}

matlist T

********************************************************************************
** Table 2: Nonparametric Density Continuity Test
********************************************************************************

mat T = J(3,5,.)

rddensity R
mat T[1,1] = round(e(h_l),.001)
mat T[1,2] = round(e(h_r),.001)
mat T[1,3] = e(N_h_l)
mat T[1,4] = e(N_h_r)
mat T[1,5] = round(e(pv_q),.001)

rddensity R, bwselect(diff) 
mat T[2,1] = round(e(h_l),.001)
mat T[2,2] = round(e(h_r),.001)
mat T[2,3] = e(N_h_l)
mat T[2,4] = e(N_h_r)
mat T[2,5] = round(e(pv_q),.001)

rddensity R, fitselect(restricted)
mat T[3,1] = round(e(h_l),.001)
mat T[3,2] = round(e(h_r),.001)
mat T[3,3] = e(N_h_l)
mat T[3,4] = e(N_h_r)
mat T[3,5] = round(e(pv_q),.001)

matlist T

********************************************************************************
** Table 3: Flexible Parametric RD Methods
********************************************************************************

mat T = J(10,4,.)
mat CI = J(7,3,.)

gen D = $R >= $c

gen R1 = ($R-$c)
gen R2 = ($R-$c)^2
gen R3 = ($R-$c)^3
gen R4 = ($R-$c)^4

gen DR1 = D*R1 
gen DR2 = D*R2
gen DR3 = D*R3
gen DR4 = D*R4

mat T[1,1] = 1
mat T[1,2] = 1
mat T[1,3] = 4
mat T[1,4] = 4
mat T[2,1] = 9
mat T[2,2] = 18
mat T[2,3] = 20
mat T[2,4] = 100

** Outcome Variable
qui count if abs(R1)<=9 & D==0 & $Y!=. & $R!=.
mat T[7,1] = r(N)
qui count if abs(R1)<=9 & D==1 & $Y!=. & $R!=.
mat T[8,1] = r(N)
reg $Y D R1 DR1 if abs(R1)<=9 , robust
mat aux = r(table)
mat T[3,1] = round(aux[1,1],.001)
mat T[4,1] = round(aux[5,1],.001)
mat T[5,1] = round(aux[6,1],.001)
mat T[6,1] = round(aux[4,1],.001)
mat CI[5,1] = round(aux[5,1],.001)
mat CI[5,2] = round(aux[1,1],.001)
mat CI[5,3] = round(aux[6,1],.001)

qui count if abs(R1)<=18 & D==0 & $Y!=. & $R!=.
mat T[7,2] = r(N)
qui count if abs(R1)<=18 & D==1 & $Y!=. & $R!=.
mat T[8,2] = r(N)
reg $Y D R1 DR1 if abs(R1)<=18, robust
mat aux = r(table)
mat T[3,2] = round(aux[1,1],.001)
mat T[4,2] = round(aux[5,1],.001)
mat T[5,2] = round(aux[6,1],.001)
mat T[6,2] = round(aux[4,1],.001)
mat CI[6,1] = round(aux[5,1],.001)
mat CI[6,2] = round(aux[1,1],.001)
mat CI[6,3] = round(aux[6,1],.001)

qui count if abs(R1)<=20 & D==0 & $Y!=. & $R!=.
mat T[7,3] = r(N)
qui count if abs(R1)<=20 & D==1 & $Y!=. & $R!=.
mat T[8,3] = r(N)
reg $Y D R1-R4 DR1-DR4 if abs(R1)<=20 , robust
mat aux = r(table)
mat T[3,3] = round(aux[1,1],.001)
mat T[4,3] = round(aux[5,1],.001)
mat T[5,3] = round(aux[6,1],.001)
mat T[6,3] = round(aux[4,1],.001)


qui count if abs(R1)<=100 & D==0 & $Y!=. & $R!=.
mat T[7,4] = r(N)
qui count if abs(R1)<=100 & D==1 & $Y!=. & $R!=.
mat T[8,4] = r(N)
reg $Y D R1-R4 DR1-DR4 if abs(R1)<=100, robust
mat aux = r(table)
mat T[3,4] = round(aux[1,1],.001)
mat T[4,4] = round(aux[5,1],.001)
mat T[5,4] = round(aux[6,1],.001)
mat T[6,4] = round(aux[4,1],.001)
mat CI[7,1] = round(aux[5,1],.001)
mat CI[7,2] = round(aux[1,1],.001)
mat CI[7,3] = round(aux[6,1],.001)

** Placebo Outcome Variables
local count=9
foreach y of varlist mort_age59_injury_postHS mort_age59_related_preHS {
	reg `y' D R1 DR1 if abs(R1)<=9 , robust
	mat aux = r(table)
	mat T[`count',1] = round(aux[4,1],.001)
	reg `y' D R1 DR1 if abs(R1)<=18, robust
	mat aux = r(table)
	mat T[`count',2] = round(aux[4,1],.001)
	reg `y' D R1-R4 DR1-DR4 if abs(R1)<=20 , robust
	mat aux = r(table)
	mat T[`count',3] = round(aux[4,1],.001)
	reg `y' D R1-R4 DR1-DR4 if abs(R1)<=100, robust
	mat aux = r(table)
	mat T[`count',4] = round(aux[4,1],.001)
	local ++count
}

matlist T

********************************************************************************
** Table 4: Robust Nonparametric Local Polynomial Methods
********************************************************************************

mat T = J(10,4,.)

** rdrobust: MSE data-driven bandwidth - local constant regression
rdrobust $Y R, p(0)
mat aux = e(b)
mat T[1,1] = e(p)
mat T[2,1] = round(e(h_l),.001)
mat T[3,1] = round(aux[1,1],.001)
local lb = e(ci_l_rb)
local ub = e(ci_r_rb)
mat T[4,1] = round(`lb',.001)
mat T[5,1] = round(`ub',.001)
mat T[6,1] = round(e(pv_rb),.001)
mat T[7,1] = e(N_h_l)
mat T[8,1] = e(N_h_r)
mat CI[3,1] = round(`lb',.001)
mat CI[3,2] = round(aux[1,1],.001)
mat CI[3,3] = round(`ub',.001)

** rdrobust: parametric bandwidth - local constant regression
rdrobust $Y R, p(0) h(9)
mat aux = e(b)
mat T[1,2] = e(p)
mat T[2,2] = round(e(h_l),.001)
mat T[3,2] = round(aux[1,1],.001)
local lb = e(ci_l_rb)
local ub = e(ci_r_rb)
mat T[4,2] = round(`lb',.001)
mat T[5,2] = round(`ub',.001)
mat T[6,2] = round(e(pv_rb),.001)
mat T[7,2] = e(N_h_l)
mat T[8,2] = e(N_h_r)

** rdrobust: MSE data-driven bandwidth - local linear regression
rdrobust $Y R, p(1)
mat aux = e(b)
mat T[1,3] = e(p)
mat T[2,3] = round(e(h_l),.001)
mat T[3,3] = round(aux[1,1],.001)
local lb = e(ci_l_rb)
local ub = e(ci_r_rb)
mat T[4,3] = round(`lb',.001)
mat T[5,3] = round(`ub',.001)
mat T[6,3] = round(e(pv_rb),.001)
mat T[7,3] = e(N_h_l)
mat T[8,3] = e(N_h_r)
local lb = e(ci_l_rb)
local ub = e(ci_r_rb)
mat CI[4,1] = round(`lb',.001)
mat CI[4,2] = round(aux[1,1],.001)
mat CI[4,3] = round(`ub',.001)

** rdrobust: parametric bandwidth - local linear regression
rdrobust $Y R, p(1) h(9)
mat aux = e(b)
mat T[1,4] = e(p)
mat T[2,4] = round(e(h_l),.001)
mat T[3,4] = round(aux[1,1],.001)
local lb = e(ci_l_rb)
local ub = e(ci_r_rb)
mat T[4,4] = round(`lb',.001)
mat T[5,4] = round(`ub',.001)
mat T[6,4] = round(e(pv_rb),.001)
mat T[7,4] = e(N_h_l)
mat T[8,4] = e(N_h_r)

** Placebo Outcome Variables
local count = 9
foreach y of varlist mort_age59_injury_postHS mort_age59_related_preHS {
   rdrobust `y' R, p(0)
   mat T[`count',1] = round(e(pv_rb),.001)
   rdrobust `y' R, p(0) h(9)
   mat T[`count',2] = round(e(pv_rb),.001)
   rdrobust `y' R, p(1)
   mat T[`count',3] = round(e(pv_rb),.001)
   rdrobust `y' R, p(1) h(9)
   mat T[`count',4] = round(e(pv_rb),.001)
   local ++count
}

matlist T

********************************************************************************
** Figure 2: Window Selection
********************************************************************************

gl covs60 "census1960_pop census1960_pctsch1417 census1960_pctsch534 census1960_pctsch25plus census1960_pop1417 census1960_pop534 census1960_pop25plus census1960_pcturban census1960_pctblack" 
gl covs90 "census1990_pop census1990_pop1824 census1990_pop2534 census1990_pop3554 census1990_pop55plus census1990_pcturban census1990_pctblack census1990_percapinc"
gl wreps = 1000
gl rreps = 5000

** Window Selection
rdwinselect R mort_age59_related_preHS $covs60, reps($wreps) stat(ksmirnov) wmin(.3) wstep(.2) level(.2)

** Generate p-values plot 
** NOTE: the plot is drawn using the asymptotic p-value to speed up the process.
** Remove the "approx" option to use randinf and replicate the results in the paper.

rdwinselect R mort_age59_related_preHS $covs60, reps($wreps) stat(ksmirnov) nwin(40) wmin(.3) wstep(.2) level(.2) plot approx
mat Res = r(results)
preserve
svmat Res
rename Res1 pvalues 
rename Res6 w
gen red=pval if Res3==43
twoway(scatter pval w)(scatter red w, msize(vlarge) msymbol(circle_hollow) mlwidth(medthick)), ///
	xline(1.1,lpattern(shortdash)) ytitle(p-values) xtitle(bandwidth) ///
	xlabel(0.3(.4)8.1, labsize(small)) legend(off) graphregion(color(white))
restore

gl w0 = 1.1

** Scatter Plot with Means 
tempvar mt mc
qui sum $Y if abs(R)<=$w0 & D==1
gen `mt'=r(mean) if D==1
qui sum $Y if abs(R)<=$w0 & D==0
gen `mc'=r(mean) if D==0
local lab: variable label $Y
twoway (scatter $Y R if abs(R)<=$w0) (line `mt' R if abs(R)<=$w0, sort lcolor(black) lpattern(shortdash)) ///
	(line `mc' R if abs(R)<=$w0, sort lcolor(black) lpattern(shortdash)),  ///
	graphregion(color(white)) title("") legend(off) ytitle("`lab'")

********************************************************************************
** Table 5: Local Randomization Methods
********************************************************************************

mat T = J(8,6,.)
gl w0 = 1.1

rdrandinf $Y R, wl(-$w0) wr($w0) reps($rreps)
mat T[1,1] = r(p)
mat T[2,1] = round(r(wr),.001)
mat T[3,1] = round(r(obs_stat),.001)
mat T[4,1] = round(r(randpval),.001)
mat T[5,1] = r(N_left)
mat T[6,1] = r(N_right)

rdrandinf $Y R, wl(-3.235) wr(3.235) reps($rreps)
mat T[1,2] = r(p)
mat T[2,2] = round(r(wr),.001)
mat T[3,2] = round(r(obs_stat),.001)
mat T[4,2] = round(r(randpval),.001)
mat T[5,2] = r(N_left)
mat T[6,2] = r(N_right)

rdrandinf $Y R, wl(-9) wr(9) reps($rreps)
mat T[1,3] = r(p)
mat T[2,3] = round(r(wr),.001)
mat T[3,3] = round(r(obs_stat),.001)
mat T[4,3] = round(r(randpval),.001)
mat T[5,3] = r(N_left)
mat T[6,3] = r(N_right)

rdrandinf $Y R, wl(-$w0) wr($w0) reps($rreps) p(1)
mat T[1,4] = r(p)
mat T[2,4] = round(r(wr),.001)
mat T[3,4] = round(r(obs_stat),.001)
mat T[4,4] = round(r(randpval),.001)
mat T[5,4] = r(N_left)
mat T[6,4] = r(N_right)

rdrandinf $Y R, wl(-3.235) wr(3.235) reps($rreps) p(1)
mat T[1,5] = r(p)
mat T[2,5] = round(r(wr),.001)
mat T[3,5] = round(r(obs_stat),.001)
mat T[4,5] = round(r(randpval),.001)
mat T[5,5] = r(N_left)
mat T[6,5] = r(N_right)

rdrandinf $Y R, wl(-9) wr(9) reps($rreps) p(1)
mat T[1,6] = r(p)
mat T[2,6] = round(r(wr),.001)
mat T[3,6] = round(r(obs_stat),.001)
mat T[4,6] = round(r(randpval),.001)
mat T[5,6] = r(N_left)
mat T[6,6] = r(N_right)

** Placebo outcomes
local i=7
foreach var of varlist mort_age59_injury_postHS mort_age59_related_preHS {
	rdrandinf `var' R, wl(-$w0) wr($w0) reps($rreps)
	mat T[`i',1] = round(r(randpval),.001)
	rdrandinf `var' R, wl(-3.235) wr(3.235) reps($rreps)
	mat T[`i',2] = round(r(randpval),.001)
	rdrandinf `var' R, wl(-9) wr(9) reps($rreps)
	mat T[`i',3] = round(r(randpval),.001)

	rdrandinf `var' R, wl(-$w0) wr($w0) reps($rreps) p(1)
	mat T[`i',4] = round(r(randpval),.001)
	rdrandinf `var' R, wl(-3.235) wr(3.235) reps($rreps) p(1)
	mat T[`i',5] = round(r(randpval),.001)
	rdrandinf `var' R, wl(-9) wr(9) reps($rreps) p(1)
	mat T[`i',6] = round(r(randpval),.001)
	
	local ++i
}

matlist T

********************************************************************************
** Figure 3: Sensitivity to bandwidth choice
********************************************************************************

rdsensitivity $Y R, wlist(0.3(0.2)10.1) tlist(-7(.25)2) saving(graphdata)
preserve
use graphdata, clear
twoway contour pvalue t w, ccuts(0(0.05)1) ccolors(gray*0.01 gray*0.05 ///
	gray*0.1 gray*0.15 gray*0.2 gray*0.25 gray*0.3 gray*0.35 ///
	gray*0.4 gray*0.5 gray*0.6 gray*0.7 gray*0.8 gray*0.9 gray ///
	black*0.5  black*0.6 black*0.7 black*0.8 black*0.9 black) ///
	xlabel(.3(1)10.1, labsize(small)) ylabel(-10(1.5)6, nogrid labsize(small)) ///
	graphregion(color(white)) ytitle("null hypothesis") xtitle(bandwidth)
restore

rdsensitivity $Y R, wlist(0.3(0.2)10.1) tlist(-7(.25)2) p(1) saving(graphdata_p1)
preserve
use graphdata_p1, clear
twoway contour pvalue t w, ccuts(0(0.05)1) ccolors(gray*0.01 gray*0.05 ///
	gray*0.1 gray*0.15 gray*0.2 gray*0.25 gray*0.3 gray*0.35 ///
	gray*0.4 gray*0.5 gray*0.6 gray*0.7 gray*0.8 gray*0.9 gray ///
	black*0.5  black*0.6 black*0.7 black*0.8 black*0.9 black) ///
	xlabel(.3(1)10.1, labsize(small)) ylabel(-10(1.5)6, nogrid labsize(small)) ///
	graphregion(color(white)) ytitle("null hypothesis") xtitle(bandwidth)
restore

********************************************************************************
** Table 6: Local Randomization Methods -- CI and Interference
********************************************************************************

mat T = J(10,2,.)

rdrandinf $Y R, wl(-$w0) wr($w0) reps($rreps) interfci(.05) p(0)
mat T[1,1] = r(p)
mat T[2,1] = round(r(wr),.001)
mat T[3,1] = round(r(obs_stat),.001)
mat T[4,1] = round(r(randpval),.001)
mat T[7,1] = round(r(int_lb),.001)
mat T[8,1] = round(r(int_ub),.001)
mat T[9,1] = r(N_left)
mat T[10,1] = r(N_right)
mat CI[1,2] = round(r(obs_stat),.001)
rdsensitivity $Y R, wlist($w0) wlist_left(-$w0) tlist(-5(.025)0) ci(-$w0 $w0) reps($rreps) nodraw 
mat T[5,1] = round(r(ci_lb),.001)
mat T[6,1] = round(r(ci_ub),.001)
mat CI[1,1] = round(r(ci_lb),.001)
mat CI[1,3] = round(r(ci_ub),.001)

rdrandinf $Y R, wl(-$w0) wr($w0) reps($rreps) interfci(.05) p(1)
mat T[1,2] = r(p)
mat T[2,2] = round(r(wr),.001)
mat T[3,2] = round(r(obs_stat),.001)
mat T[4,2] = round(r(randpval),.001)
mat T[7,2] = round(r(int_lb),.001)
mat T[8,2] = round(r(int_ub),.001)
mat T[9,2] = r(N_left)
mat T[10,2] = r(N_right)
mat CI[2,2] = round(r(obs_stat),.001)
rdsensitivity $Y R, wlist($w0) wlist_left(-$w0) tlist(-5(.025)0) ci(-$w0 $w0) reps($rreps) nodraw p(1)
mat T[5,2] = round(r(ci_lb),.001)
mat T[6,2] = round(r(ci_ub),.001)
mat CI[2,1] = round(r(ci_lb),.001)
mat CI[2,3] = round(r(ci_ub),.001)

matlist T

********************************************************************************
** Table 7: Local Randomization Methods -- Rosenbaum Bounds
********************************************************************************

rdrbounds $Y R, expgamma(1.1 1.2 1.3 1.4) wlist(.3 .5 .7 .9 1.1 1.3 1.5) ///
				statistic(ttest) fmp bound(upper) reps(5000)
mat Ta = r(pvals)
mat Tb = r(ubound)
mat Tc = (.09,.18,.26,.34 \ 1.1,1.2,1.3,1.4)'
mat T = (.,.,.3,.5,.7,.9,1.1,1.3,1.5 \ J(2,2,.),Ta \ Tc,Tb)

matlist T

********************************************************************************
** Figure 4: Summary of Results
********************************************************************************

mat H = (1.0,1.2,3.235,6.811,9,18,25)'
mat T = (H,CI[1...,2],CI[1...,1],CI[1...,3])

matlist T

** Plot all results
preserve
svmat CI
svmat H
rename CI1 lb
rename CI2 pe
rename CI3 ub
rename H1 h

twoway (rcap lb ub h)(scatter pe h), yline(0, lcolor(black) lpattern(dash)) ///
	xlabel(1.1 "LR" 3.235 "MSE_0" 6.811 "MSE_1" 9 "LM9" 18 "LM18" 25 "GP", labsize(small)) ///
	xscale(range(-2 30)) xline(2.2285 7.699 21.5, lpattern(shortdash) lcolor(black)) ///
	legend(off) graphregion(color(white)) yscale(range(-10 2)) ylabel(-6(2)2, labsize(small)) ///
	xtitle(bandwidth) ytitle("RD inference results") ///
	text(-8 0 "Local" "Rand.") ///
	text(-8 5 "Non" "Parametric") ///
	text(-8 14 "Flexible" "Parametric") ///
	text(-8 25 "Global" "Parametric")
restore

********************************************************************************
** Additional empirical analysis
********************************************************************************

** NOTE: this analysis is not reported in the paper.

** Robust Nonparametric Methods with Covariates

rdrobust $Y R, covs($covs60)

** Robust Nonparametric Methods: Different Bandwdiths at Each Side

rdrobust $Y R, bwselect(msetwo)
