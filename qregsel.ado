*! version 1.1 24July2020
* Authors: Ercio Munoz & Mariel Siravegna 

/* Wrap file to implement copula-based sample selection model in quantile 
regression as suggested by Arellano and Bonhomme (2017). */

* qregselection wage educ age, quantile(.1 .5 .9) select(married children educ age) py grid_min(-.99) grid_max(.99) grid_length(.02)

cap program drop qregsel
program define qregsel, eclass sortpreserve
    version 16.0

    syntax varlist(numeric) [if] [in], SELect(string) quantile(string) ///
	grid_min(real) grid_max(real) grid_length(real) ///
	[ copula(string) noCONStant plot ]
    
    gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'

    fvexpand `indepvars' 
    local cnames `r(varlist)'
	
	tokenize `select', parse("=")
	if "`2'" != "=" {
		local x_s `select'
		tempvar y_s
		qui gen `y_s' = (`depvar'!=.)
	}
	else {
		local y_s `1'
		local x_s `3'
	}		
	capture unab x_s : `x_s'
	
********************************************************************************	
** Marking the sample to use (selected observations)
********************************************************************************	
	marksample touse
	markout `touse' `y_s' `x_s'
	tempvar touse1
	mark `touse1' if `y_s'==1
	markout `touse1' `depvar' `indepvars'
	qui replace `touse' = 0 if `touse1'==0 & `y_s'==1 

	
********************************************************************************	
** Checking errors in the selection indicator
********************************************************************************	
	qui tab `y_s' 
	if r(r) != 2 {
	dis as error "`y_s' should be binary."
	exit 198 
	}
	qui sum `y_s' 
	if r(max) != 1 & r(min) != 0 {
	dis as error "`y_s' should be either 0 or 1."
	exit 198
	}
	
	
********************************************************************************	
** Checking errors with copula specification (default is gaussian)
********************************************************************************	
	if "`copula'" == "" {
	local copula "gaussian"
	} 			
	if  "`copula'" != "gaussian" & "`copula'" != "fgm" & /// 
		"`copula'" != "amh" & "`copula'" != "frank"	  { 
		dis as error "`copula' is not an available copula."
		exit 198
	}
	

********************************************************************************	
** Checking errors in the grid specification
********************************************************************************	
	if `grid_max'<`grid_min' {
	dis as error "The maximum value of the grid must be greater or equal than the minimum"
	exit 198 
	} 
	if (`grid_max'==`grid_min' & `grid_length'==0) {
	local grid_rho = 1
    }
	else {
	local grid_rho = (`grid_max'-`grid_min')/`grid_length' + 1 
	}
	if mod(`grid_rho',1) !=0 {
	dis as error "Specify a length that allows an evenly spaced grid"
	exit 198 
	}
	if (`grid_max'>=1 | `grid_min'<=-1) & ///
	( "`copula'" == "gaussian" | "`copula'" == "fgm" & /// 
		"`copula'" == "amh" ) {
	dis as error "Grid values are out of range"
	exit 198 
	} 

	
********************************************************************************
** Generate the propensity score and the instrument	
********************************************************************************	
	tempname pZ varphi1 b N rank df_r it obj_j index object betas_taus rhos x_grid
	tempvar copula_cdf
	qui: probit `y_s' `x_s'
	qui: predict `pZ'
	qui: gen `varphi1' = `pZ' if `y_s'==1


********************************************************************************
** Estimate rho looping over quantiles and the grid for rho	
********************************************************************************		
	mat `object' = J(`grid_rho',1,0)
	mat `x_grid' = J(`grid_rho',1,0)

forvalues j = 1(1)`grid_rho' {

local rhoa = `grid_min' + (`j'-1)*`grid_length'
mat `x_grid'[`j',1] = `rhoa'
local obj = 0

* Loop over 9 quantiles to compute values of the grid to be minimized
forvalues k = 1(1)9 {
local tau = `k'/10

** Obtain copula
qui:	mata: copulafn("`pZ'",`rhoa',`tau',"`touse'","`copula_cdf'","`copula'")

** Rotated quantile regression
qui:	mata: mywork("`depvar'", "`cnames'", "`touse'", "`constant'", ///
	"`b'", "`N'", "`rank'", "`df_r'","`it'","`copula_cdf'")

** Obtain value for the grid to be minimized
qui:	mata: objective("`pZ'",`rhoa',`tau',"`touse'", ///
 	"`depvar'", "`cnames'", "`constant'","`obj_j'","`b'","`copula_cdf'")	
qui: cap drop `copula_cdf'
local obj = `obj'+`obj_j'
}

mat `object'[`j',1] = `obj'*`obj'
}
numlist "`grid_min'(`grid_length')`grid_max'"
local rownames = r(numlist)
mat rownames `object' = `rownames'
mat colnames `object' = grid_rho

** Minimize the objective function
qui: mata: minmatrix("`object'","`index'")	
local index_min = `index'[1,1]


********************************************************************************	
** Plot the grid used for minimization
********************************************************************************
if "`plot'" != "" {
tempvar rho_eval1 rho_evali1
svmat `object', names(`rho_eval1')
svmat `x_grid', names(`rho_evali1')
twoway line `rho_eval1' `rho_evali1' if _n<=`grid_rho', name(rho_plot,replace) scheme(s1color) xtitle(rho) ytitle(objective)
}


********************************************************************************	
** Checking errors in the quantiles requested
********************************************************************************
tokenize `quantile'
local orig `1'
macro shift
if "`orig'" == "" {
	di in red "option quantile() required"
	exit 198
}
capture confirm number `orig'
if _rc {
	di "`orig' not a number"
	exit 198
}
if `orig' >= 1 {
	local orig = `orig'/100
}
if `orig'<=0 | `orig' >=1 {
	local orig = 100*`orig'
	di "`orig' out of range"
	exit 198
}
local quants = `orig'

while "`1'" != "" {
local orig `1'
macro shift
if "`orig'" == "" {
	di in red "option quantile() required"
	exit 198
}
capture confirm number `orig'
if _rc {
	di "`orig' not a number"
	exit 198
}
if `orig' >= 1 {
	local orig = `orig'/100
}
if  `orig' >=1 {
	local orig = 100*`orig'
	di "`orig' out of range"
	exit 198
}
if `orig'<=0 {
	di "`orig' out of range"
	exit 198
}
local quants = "`quants' `orig'"
}


********************************************************************************	
** Estimate rotated quantile regression using the selected rho
********************************************************************************
local count: word count `indepvars' 
mat `betas_taus' = J(`count'+1,1,.)
foreach tau of local quants {
local rho = `grid_min' + (`index_min'-1)*`grid_length'

** Obtain copula
qui:	mata: copulafn("`pZ'",`rho',`tau',"`touse'","`copula_cdf'","`copula'")

** Rotated quantile regression
qui:	mata: mywork("`depvar'", "`cnames'", "`touse'", "`constant'", ///
	"`b'", "`N'", "`rank'", "`df_r'","`it'","`copula_cdf'")
qui: cap drop `copula_cdf'

local qtau = 100*`tau'
mat colnames `b' = "q`qtau'"
mat `betas_taus' = `betas_taus',`b'
}

mat `betas_taus' = `betas_taus'[1...,2...]
    if "`constant'" == "" {
    	local cnames `cnames' _cons
    }
matrix rownames `betas_taus' = `cnames'	


********************************************************************************	
** Generating the output
********************************************************************************	
	dis " "
	dis in green "Quantile selection model" ///
		_column (50) "Number of obs" _column(69) "=" _column(71) %8.0f in yellow `N'

    ereturn post , esample(`touse') buildfvinfo
	ereturn matrix grid    = `object'
    ereturn matrix coefs   = `betas_taus'
    ereturn scalar N       = `N'
	ereturn scalar rank    = `rank'
    ereturn scalar df_r    = `df_r'
    ereturn scalar rho     = `rho'
	ereturn scalar it      = `it'
	ereturn local title   "Quantile selection model"
	ereturn local predict "qregsel_p"
    ereturn local cmd     "qregsel"
	ereturn local selection_eq "`select'"	
	ereturn local outcome_eq "`depvar' `indepvars'"
	ereturn local cmdline "qregsel `depvar' `indepvars', select(`select')"
	ereturn local indepvars "`cnames'"
	ereturn local depvar  "`depvar'"
	ereturn local copula "`copula'"

    ereturn display
	matlist e(coefs)
	
end

********************************************************************************	
** Auxiliary functions needed for the estimation
********************************************************************************
mata:

void copulafn( string scalar pscore, numeric vector rho,
			   numeric vector tau,	 string scalar touse,   
			   string scalar cdf,    string scalar name) 
{

    real matrix pZ1, G, vs, v1

    pZ1  = st_data(., pscore, touse)
	
	if (name=="gaussian") {
	vs = J(rows(pZ1),1,invnormal(tau))
	v1 = invnormal(pZ1)
	st_view(G, ., st_addvar("float", cdf),touse)
	G[.,.] = binormal(vs,v1,rho) :/ pZ1
	}
	else if (name=="fgm") {
	st_view(G, ., st_addvar("float", cdf),touse)
	G[.,.] = (tau*pZ1) :* ( 1:+rho*(1-tau):*(1:-pZ1) ) :/ pZ1
	}
	else if (name=="frank") {
	st_view(G, ., st_addvar("float", cdf),touse)
	G[.,.] = -ln(1:+(exp(-rho*tau):-1):*(exp(-rho:*pZ1):-1):/(exp(-rho)-1)):/(rho:*pZ1)
	}
	else {
	st_view(G, ., st_addvar("float", cdf),touse)
	G[.,.] = ( tau*pZ1:/(1:-rho*(1:-pZ1):*(1-tau)) ) :/ pZ1
	}
		
}

void mywork( string scalar depvar,  	string scalar indepvars, 
             string scalar touse,   	string scalar constant,  
			 string scalar bname,       string scalar nname,   
			 string scalar rname,       string scalar dfrname, 
			 string scalar itname, 		string scalar cdf) 
{
    real vector y, p
    real matrix X, u, a, b
    real scalar m, n, k, it

    y    = st_data(., depvar, touse)
    X    = st_data(., indepvars, touse)
 	p    = st_data(., cdf, touse) 
	
	m    = rows(X)
	n    = cols(X)
	
    if (constant == "") {
    X    = X,J(m,1,1)
    }
	k    = cols(X) 
	
	u    = J(m, 1, 1)
	a    = (1:-p):*u
	it=0
	b    = -lp_fnm(X',-y',X'*a,u,a,it)'
	
	st_matrix(bname, b)
    st_numscalar(itname, it)
    st_numscalar(nname, m)
    st_numscalar(rname, k)
    st_numscalar(dfrname, m-k)
	
}

void objective( string scalar pscore, 	numeric vector rhoa,
				numeric vector tau, 	string scalar touse,   
				string scalar depvar,   string scalar indepvars, 
				string scalar constant, string scalar G2, 
				string scalar betas   , string scalar G1						) 
{

    real matrix varphi1, pZ1, G, y, X, b
	real scalar n
	
	y    = st_data(., depvar, touse)
    X    = st_data(., indepvars, touse)
	varphi1  = st_data(., pscore, touse)
	pZ1  = st_data(., pscore, touse)'
    n    = rows(X)
	b = st_matrix(betas)'
	copula_p = st_data(., G1, touse)

    if (constant == "") {
        X    = X,J(n,1,1)
    }
	
	st_numscalar(G2, mean( varphi1 :* ((y:<=X*b'):-copula_p) ))
	
}

void minmatrix(string scalar obj, string scalar G2) 
{
    real matrix values, i, w
	values = st_matrix(obj)
	i = J(0,0,.)
	w = J(0,0,.)
	minindex(values,1,i,w)
	st_matrix(G2, i)
	
}

end

cap mata mata drop lp_fnm()
mata:
real matrix function lp_fnm(real matrix A,
							real matrix c,
							real matrix b,
							real matrix u,
							real matrix x,
							real scalar it)
{ 

  beta = 0.9995
  small = 1e-5
  max_it = 50
  m = rows(A)
  n = cols(A)

// Generate initial feasible point 
  s = u - x  
  y = svsolve(A',c')'
  r = c - y * A
  r = mm_cond(r:==0,r:+0.001,r)  
  z = mm_cond(r:>0,r,0)
  w = z - r
  gap = c * x - y * b + w * u

// Start iterations
 it = 0
while (gap > small & it < max_it) {
    it++

// Compute affine step
    q = 1 :/ (z' :/ x + w' :/ s)
    r = z - w
    Q = SPMATbandedmake(diag(sqrt(q)),0,0)
    AQ = SPMATbandedmultfull(Q,0,0,A')'
	rhs = SPMATbandedmultfull(Q,0,0,r')
    dy = (svsolve(AQ',rhs))'
    dx = q :* (dy * A - r)'
    ds = -dx
    dz = -z :* (1 :+ dx :/ x)'
    dw = -w :* (1 :+ ds :/ s)'

// Compute maximum allowable step lengths
    fx = bound(x, dx)
    fs = bound(s, ds)
    fw = bound(w, dw)
    fz = bound(z, dz)
    fp = mm_cond(fx:<fs,fx,fs) 
	fd = mm_cond(fw:<fz,fw,fz)
	fp = mm_cond(min(beta * fp):<1,min(beta * fp),1) 
	fd = mm_cond(min(beta * fd):<1,min(beta * fd),1) 
	
if (mm_cond(fp:<fd,fp,fd) < 1) {
    
// Update mu
      mu = z * x + w * s
      g = (z + fd * dz) * (x + fp * dx) + (w + fd * dw) * (s + fp * ds)
      mu = mu * (g / mu) ^3 / ( 2 * n)

// Compute modified step
      dxdz = dx :* dz'
      dsdw = ds :* dw'
      xinv = 1 :/ x
      sinv = 1 :/ s
      xi = mu * (xinv - sinv)
	  rhs = rhs + SPMATbandedmultfull(Q,0,0,( dxdz - dsdw -xi ))
	  dy = (svsolve(AQ',rhs))'
      dx = q :* (A' * dy' + xi - r' -dxdz + dsdw)
      ds = -dx
      dz = mu * xinv' - z - xinv' :* z :* dx' - dxdz'
      dw = mu * sinv' - w - sinv' :* w :* ds' - dsdw'

// Compute maximum allowable step lengths
      fx = bound(x, dx)
      fs = bound(s, ds)
      fw = bound(w, dw)
      fz = bound(z, dz)
	  fp = mm_cond(fx:<fs,fx,fs) 
	  fd = mm_cond(fw:<fz,fw,fz)
	  fp = mm_cond(min(beta * fp):<1,min(beta * fp),1) 
	  fd = mm_cond(min(beta * fd):<1,min(beta * fd),1) 
}

// Take the step
    x = x + fp * dx
    s = s + fp * ds
    y = y + fd * dy
    w = w + fd * dw
    z = z + fd * dz
    gap = c * x - y * b + w * u
}
return(y)
}
end
********************************************************************************
// Bound function
********************************************************************************
cap mata: mata drop bound()
mata:
real matrix function bound(real matrix x, real matrix dx)
{
	return(mm_cond(dx:<0,-x:/dx,1e20 :+ 0 :* x))
}
end




