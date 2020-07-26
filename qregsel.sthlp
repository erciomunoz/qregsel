{smcl}
{* *! version 1.1 July 2020}{...}
{cmd: help qregsel}
{hline}

{title:Title}

{phang}
{bf:Quantile Regression Corrected for Sample Selection}

{title:Syntax}

{p 8 17 2}
{cmd:qregsel}
{it:depvar} {it:varlist}
{ifin} 
{cmd:,}
{cmdab:sel:ect(}[{it:depvar_s} {cmd:=}] {it:varlist_s}{cmd:)}
{cmd:quantile(}{it:#}{cmd:)}
{cmd:grid_min(}{it:#}{cmd:)}
{cmd:grid_max(}{it:#}{cmd:)}
{cmd:grid_length(}{it:#}{cmd:)}
[
{cmd:copula(}{it:copula}{cmd:)}
{cmdab:nocons:tant}
{cmdab:plot}
]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt sel:ect()}}specifies a selection equation.{p_end}

{synopt:{opt quantile:(#)}}specifies the quantiles to be estimated.{p_end}

{synopt:{opt grid_min:(#)}}specifies the minimum value to be considered in the grid search.{p_end}

{synopt:{opt grid_max:(#)}}specifies the maximum value to be considered in the grid search.{p_end}

{synopt:{opt grid_length:(#)}}specifies the length of the (evenly spaced) grid to be considered in the grid search.{p_end}

{synopt:{opt copula:(copula)}}specifies a copula;
	default is gaussian.{p_end}

{synopt:{opt nocons:tant}}suppresses a constant term in the outcome equation.{p_end}

{synopt:{opt plot:}}generates a line graph of the value of the objective function over rho's grid.{p_end}



{title:Description}

{pstd}
{cmd:qregsel} estimates a copula-based sample selection model for quantile regression.	
Users can specify a copula from the lists below. 

{p 4 4 2}
Available copulas are 

{p 8 8 2} {it: gaussian, fgm, amh}, and {it:frank}.

{p 4 4 2} 
Notes: The name of the copula is case-sensitive. 


{title:Options}

{dlgtab:Main}

{phang}
{opt sel:ect()} is required. It specifies a selection equation. 
If {it:depvar_s} is specified, it should be coded as 0 or 1, which 0 indicating {it:depvar} not observed for an observation 
and 1 indicating {it:depvar} observed for an observation. 

{phang} 
{opt quantile(#)} is required. It specifies a set of quantiles to be estimated.

{phang} 
{opt grid_min:(#)} specifies the minimum value to be considered in the grid search. grid_min() is required. It must be greater than -1 when copula is gaussian, fgm, or amh.

{phang} 
{opt grid_max:(#)} specifies the maximum value to be considered in the grid search. grid_max() is required and cannot be smaller than grid_min(). It must be smaller than 1 when copula is gaussian, fgm, or amh.

{phang} 
{opt grid_length:(#)} specifies the length of the (evenly spaced) grid to be considered in the grid search. grid_length() is required. It can be set to 0 if grid_min() and grid_max() are equal.

{phang}
{opt copula(copula)} specifies a copula function for the dependence between outcome and selection equation.
See above for the list of available copulas. 
Default is {bf:gaussian}. 

{phang}
{opt noncons:tant} suppresses a constant term of the outcome equation.

{phang}
{opt plot} generates a line graph of the value of the objective function over rho's grid. By default the grid is composed by 199 evenly space values. 



{title:Saved results}

{pstd}
{cmd:qregsel} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(rank)}}number of parameters{p_end}
{synopt:{cmd:e(df_r)}}degrees of freedom{p_end}
{synopt:{cmd:e(rho)}}copula parameter{p_end}
{synopt:{cmd:e(it)}}number of iterations{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(copula)}}specified {cmd:copula}{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(indepvar)}}independent variables{p_end}
{synopt:{cmd:e(cmdline)}}command line{p_end}
{synopt:{cmd:e(outcome_eq)}}outcome equation{p_end}
{synopt:{cmd:e(select_eq)}}selection equation{p_end}
{synopt:{cmd:e(predict)}}predict command name{p_end}
{synopt:{cmd:e(cmd)}}{cmd:qregsel}{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(coefs)}}coefficient matrix. Each column corresponds to the coefficients for a quantile{p_end}
{synopt:{cmd:e(grid)}}value of the objective function minimized over the grid{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{title:References}

{pstd}
Arellano, M., and S. Bonhomme. 2017.
"Quantile  Selection  Models  With  an  Application  to  Understanding Changes in Wage Inequality." Econometrica 85(1): 1–28.

{title:Authors}
{p}
{p_end}

{pstd}
Ercio Munoz, CUNY Graduate Center, New York, US.

{pstd}
Email: {browse "mailto:emunozsaavedra@gc.cuny.edu":emunozsaavedra@gc.cuny.edu}

{pstd}
Mariel Siravegna, Georgetown University, Washington DC, US.

{pstd}
Email: {browse "mailto:mcs92@georgetown.edu":mcs92@georgetown.edu}

{title: Also see}

{psee}
Online: {help heckman}, {help qreg}, {help heckmancopula}





