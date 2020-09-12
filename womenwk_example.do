********************************************************************************
** Example 1: Wages of women
********************************************************************************
webuse womenwk,clear

** Estimate the model
global wage_eqn wage education age
global seleqn married children education age
qregsel $wage_eqn, select($seleqn) quantile(.1 .5 .9)
ereturn list
svmat e(grid), name(col)
gen lvalue = log10(value)
twoway connected lvalue spearman

** Create counterfactual outcomes
set seed 1
predict wage_hat participation
_pctile wage_hat, nq(20)
mat qs = J(19,3,.)
forvalues i=1/19{
mat qs[`i',1] = r(r`i')
}
_pctile wage, nq(20)
forvalues i=1/19{
mat qs[`i',2] = r(r`i')
mat qs[`i',3] = `i'
}
svmat qs, name(quantiles)
twoway connected quantiles1 quantiles2 quantiles3, ///
 xtitle("Ventile") ytitle("Wage") legend(order(1 "Corrected" 2 "Uncorrected"))

** Inference with bootstrap
set seed 2
webuse womenwk,clear
capture program drop myqregsel
program myqregsel, eclass
    version 16
    tempname bb
    quietly qregsel $wage_eqn, select($seleqn) quantile(.1 .5 .9)  
	local colnames : colfullnames e(coefs)
	local rownames : rowfullnames e(coefs)
	foreach lname1 of local colnames   {
		foreach lname2 of local rownames   {
			local names = "`names' `lname1':`lname2'"
		}
	}
	mata: st_matrix("`bb'", vec(st_matrix("e(coefs)"))')
    matrix `bb' = `bb',e(rho)
	mat colnames `bb' = `names' rho:rho
    ereturn post `bb'
	ereturn local cmd="bootstrap"
end
bootstrap _b, reps(99) nowarn: myqregsel 






