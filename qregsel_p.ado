* Authors: Ercio Munoz & Mariel Siravegna 
*! version 1.0.0 20April2020

*** predict command for qregsel ***
cap program drop qregsel_p
#delimit ;

program define qregsel_p;
	version 16.0;
	
	syntax newvarname [if] [in];
	
	marksample touse, novarlist;
		
	local copula  "`e(copula)'";
	local selection_eq = "`e(selection_eq)'";
	local outcome_eq = "`e(outcome_eq)'";	
	local indepvars = "`e(indepvars)'";	
	local py "`py'";
	local rho = e(rho);
	
	tempname coefs;
	tempvar q;

quietly {;
qregsel `outcome_eq', select(`selection_eq') quantile(
 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 
 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 
 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 
 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99) `py'
 copula(`copula') grid_min(`rho') grid_max(`rho') grid_length(0);

gen `q' = (int(uniform()*99 + 1));
mat `coefs' = e(coefs)';

quietly {;
preserve;
clear;
svmat `coefs', names(a);
gen `q' = _n;
tempfile temp1;
save "`temp1'";
restore;
merge m:1 `q' using "`temp1'", nogenerate;
local vars `e(indepvars)';
local i = 0;
foreach variable of local vars {;
if "`variable'"!= "_cons" {;
local i = 1+`i';
local y_hat "`variable'*a`i'+`y_hat'";
};
};
local i = `i'+1;
local y_hat "`y_hat'a`i'";
gen `typlist' `varlist' = `y_hat';
};
};	
end;

