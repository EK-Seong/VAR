***** KOSPI *****
clear
drop _all

import excel "kospi.xlsx", sheet("Sheet1") firstrow
rename 종가 kospi
gen time = _n
tsset time
gen R = 100*(kospi-L.kospi)/kospi

twoway (tsline R, lwidth(vthin)), ytitle("KOSPI returns rate") tlabel()
graph save "Graph" "KOSPI수익률.gph",replace

** model selection procedure **
arch R, arch(1) garch(1)
arch R, arch(1) garch(1/2)
arch R, arch(1/2) garch(1)

arch R, arch(1) garch(1)
est store GARCH_result
predict var_garch, variance

twoway (tsline var_garch, lwidth(vthin)), ytitle("Conditional Variance (GARCH(1,1)") 
graph save "Graph" "KOSPI조건부분산garch.gph" ,replace


** model selection procedure **
arch R, earch(1) egarch(1)
arch R, earch(1/2) egarch(1)
arch R, earch(1) egarch(1/2)

arch R, earch(1) egarch(1)
est store EGARCH_result
predict var_egarch, variance
twoway (tsline var_egarch) , title("Conditional Variance EGARCH(1,1)")
graph save "Graph" "KOSPI조건부분산egarch.gph" ,replace

** News Responce Function **
generate et = (_n-1850)/1850*8
predict sigma2, variance at(et 1)
line sigma2 et , m(i) c(l) title(News response function)

** GARCH-in-mean
arch R, earch(1) egarch(1) archm

*************************
***** Exchange Rate *****
*************************

clear
drop _all

import excel "ExchangeRate.xlsx", sheet("Sheet1") firstrow clear
gen t = _n
tsset t
gen d_e = 100*(erate-L.erate)/erate

twoway (tsline d_e, lwidth(vthin)), title("Exchange Rate")
graph save "Graph" "ExchangeRate.gph" , replace

** model selection procedure **
arch d_e, arch(1) garch(1)
arch d_e, arch(1) garch(1) ar(1)
arch d_e, arch(1/2) garch(1)
arch d_e, arch(1/3) garch(1)
arch d_e, arch(1/2) garch(1/2)
arch d_e, arch(1/2) garch(1/3)
arch d_e, arch(1/2) garch(1/4)
arch d_e, arch(1/2) garch(1/5)

arch d_e, arch(1/2) garch(1/4)
est store GARCH_result
predict var_garch, variance


twoway (tsline var_garch , lwidth(vthin)), title("Conditional Variance GARCH(2,4)")
graph save "Graph" "ExchangeRate_garch.gph" , replace

** model selection procedure **
arch d_e, earch(1) egarch(1)
arch d_e, earch(1) egarch(1) ar(1)
arch d_e, earch(1/2) egarch(1)
arch d_e, earch(1) egarch(1/2)
arch d_e, earch(1) egarch(1/3)

arch d_e, earch(1) egarch(1/2)
est store EGARCH_result
predict var_egarch, variance
twoway (tsline var_egarch,lwidth(vthin)) , title("Conditional Variance EGARCH(1,2)")
graph save "Graph" "ExchangeRate_egarch.gph" , replace

** News Responce Function **
generate et = (_n-2550)/2550*8
predict sigma2, variance at(et 1)
line sigma2 et , m(i) c(l) title(News response function)

** GARCH-in-mean
arch d_e, earch(1) egarch(1) archm
predict var_egarchm, variance
twoway (tsline var_egarchm)



*******************************
***** Inflation ***************
*******************************

clear
drop _all

import excel "cpi.xlsx", sheet("데이터") firstrow clear

gen t = _n
tsset t

twoway(tsline CPI)

gen inflation=400*(CPI-L.CPI)/CPI

twoway(tsline inflation)
graph save "Graph" "CPI.gph" , replace

/*
** model selection procedure **
arch inflation, arch(1) garch(1)
arch inflation, arch(1/2) garch(1)
arch inflation, arch(1/3) garch(1)
arch inflation, arch(1/2) garch(1/2)
arch inflation, arch(1/2) garch(1)
arch inflation, arch(1/2) garch(1) ar(1)
arch inflation, arch(1/2) garch(1) ar(1/2)
arch inflation, arch(1/2) garch(1) ar(1/3)
*/

/*
gen d = t<201
arch inflation, arch(1/2) garch(1) ar(1/2)
arch inflation, arch(1/2) garch(1) ar(1/2) archm
arch inflation d, arch(1/2) garch(1) ar(1/2)
drop var_garch
predict var_garch,variance
twoway (tsline var_garch) , tlabel(none)
graph save "Graph" "CPI_garch.gph" , replace
*/

** model selection procedure **
arch inflation, earch(1) egarch(1)
arch inflation, earch(1/2) egarch(1)
arch inflation, earch(1/3) egarch(1)
arch inflation, earch(1/2) egarch(1/2)
arch inflation, earch(1/2) egarch(1) ar(1)

arch inflation, earch(1/2) egarch(1)
predict var_egarch,variance
twoway (tsline var_egarch), title("Inlation Conditional Variance EGARCH(2,1)")
twoway (tsline var_egarch if t>200, lwidth(vthin)), title("Inflation") subtitle("Conditional Variance EGARCH(2,1)")

** News Responce Function **
generate et = (_n-350)/350*4
predict sigma2, variance at(et 1)
line sigma2 et , m(i) c(l) title(News response function)

** GARCH-in-mean
arch inflation, arch(1/2) garch(1) archm

predict sigmag,variance
twoway (tsline sigmag)

/*
* 1981.09부터분석함
drop if t < 201 
arch inflation, earch(1) egarch(1) ar(1/2)
predict sigmaeg,variance
twoway (tsline sigmaeg)
*/



