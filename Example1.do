
** Load data and estimate quantile regression 
webuse womenwk,clear

** Estimation using AB method
global wage_eqn wage educ age
global seleqn married children educ age
qregsel $wage_eqn, select($seleqn) quantile(.1 .5 .9) 
ereturn list
svmat e(grid), name(col)
qui gen lvalue = log10(value)
twoway connected lvalue spearman

** Prediction
set seed 1
predict wage_hat participation
_pctile wage_hat, nq(20)
mat qs = J(19,3,.)
forvalues i=1/19 {
	mat qs[`i',1] = r(r`i')
}
_pctile wage, nq(20)
forvalues i=1/19 {
	mat qs[`i',2] = r(r`i')
	mat qs[`i',3] = `i'
}
svmat qs, name(quantiles)
twoway connected quantiles1 quantiles2 quantiles3, ///
 xtitle("Ventile") ytitle("Wage") legend(order(1 "Corrected" 2 "Uncorrected"))

** Inference
bootstrap rho=e(rho) _b, reps(100) seed(2) notable: qregsel $wage_eqn, ///
	select($seleqn) quantile(.1 .5 .9)  
estat bootstrap, percentile
