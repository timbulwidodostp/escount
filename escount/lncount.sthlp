{smcl}
{* *! version 1.1.0 Aug2019}{...}
{cmd:help lncount}{right: ({browse "https://doi.org/10.1177/1536867X20953573":SJ20-3: st0612})}
{hline}

{title:Title}

{p2colset 5 16 18 2}{...}
{p2col:{cmd:lncount} {hline 2}}Maximum likelihood estimation of lognormal count-data regression model{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 15 2}
{cmd:lncount}
{depvar} [{cmd:=}] {indepvars}
{ifin} 
[{it:{help lncount##weight:weight}}]
[{cmd:,} {it:options}]

{synoptset 30}{...}
{synopthdr}
{synoptline}
{synopt:{opt model(model)}}specify distribution of count-data outcomes{p_end}
{synopt:{opt nocons:tant}}suppress constant term{p_end}
{synopt:{opt exp:osure(varname_e)}}include ln({it:varname_e}) in model with coefficient constrained to 1{p_end}
{synopt:{opt off:set(varname_o)}}include {it:varname_o} in model with coefficient constrained to 1{p_end}
{synopt:{opt const:raints(constraints)}}apply specified linear constraints{p_end}
{synopt:{opt intp:oints(#)}}use {it:#} Gauss-Hermite quadrature points; default is {cmd:intpoints(24)}{p_end}
{synopt:{opt vce:(vcetype)}}{it:vcetype} may be {cmd:robust}, {cmd:cluster}
{it:clustvar}, {cmd:oim}, or {cmd:opg}{p_end}
{synopt:{help lncount##maximize:{it:maximize_options}}}control the maximization process{p_end}
{synoptline}
{marker weight}{...}
{p 4 6 2}
{cmd:aweight}s, {cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s are
allowed; see {help weight}.{p_end}


{title:Description}

{pstd}
{cmd:lncount} estimates the parameters of a lognormal count-data model.


{title:Options}
 
{phang}
{opt model(model)} specifies the distribution for count-data outcomes.  The
default is {cmd:model(poisson)}.  The alternative specification is {cmd:nb}.

{phang}
{cmd:noconstant} suppresses a constant term of the outcome equation.

{phang}
{opt exposure(varname_e)} includes ln({it:varname_e}) in the model
with the coefficient constrained to 1.

{phang}
{opt offset(varname_o)} includes {it:varname_o} in the model with
the coefficient constrained to 1.

{phang}
{opt constraints(constraints)}; see
{manhelp estimation_options R:Estimation options}.

{phang}
{opt intpoints(#)} specifies the number of integration points to use for
integration by quadrature.  The default is {cmd:intpoints(24)}; the maximum is
{cmd:intpoints(512)}.

{phang}
{opt vce(vcetype)} specifies the type of standard error reported, which
includes types that are robust to some kinds of misspecification
({cmd:robust}), that allow for intragroup correlation ({cmd:cluster}
{it:clustvar}), and that are derived from asymptotic theory ({cmd:oim},
{cmd:opg}); see {manhelpi vce_option R}.

{marker maximize}{...}
{phang}
{it:maximize_options}:  
{opt dif:ficult},
{opt tech:nique(algorithm_spec)},
{opt iter:ate(#)},
[{cmd:no}]{opt log},
{opt tr:ace},
{opt grad:ient},
{opt showstep},
{opt hess:ian},
{opt showtol:erance},
{opt tol:erance(#)},
{opt ltol:erance(#)},
{opt nrtol:erance(#)},
{opt nonrtol:erance}, and
{opt from(init_specs)}; see {helpb maximize:[R] Maximize}.
These options are seldom used.


{title:Remarks}

{p 4 4 2}
The following predictions are available.


{title:Syntax for predict}

{p 8 17 2}
{cmd:predict} [{it:type}] 
{it:newvar} 
{ifin} 
[{cmd:,} {cmd:n} {cmd:xb}]


{title:Options for predict}

{phang}
{cmd:n}, the default, calculates the predicted number of events, which is
exp(xb)*exp(0.5*sigma^2) if neither {opt offset(varname_o)} nor 
{opt exposure(varname_e)} is specified; exp(xb + offset)*exp(0.5*sigma^2) if
{cmd:offset()} is specified; or exp(xb)exposure*exp(0.5*sigma^2) if
{cmd:exposure()} is specified.

{phang}
{cmd:xb} computes the linear prediction of the dependent variable for each
observation.

{title:Stored results}

{pstd}
{cmd:lncount} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(k_eq)}}number of equations in {cmd:e(b)}{p_end}
{synopt:{cmd:e(k_eq_model)}}number of equations in overall model test{p_end}
{synopt:{cmd:e(k_aux)}}number of auxiliary parameters{p_end}
{synopt:{cmd:e(k_dv)}}number of dependent variables{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(chi2)}}chi-squared{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(p)}}significance of model test{p_end}
{synopt:{cmd:e(n_quad)}}number of quadrature points{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:lncount}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(wtype)}}weight type{p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(offset)}}offset for outcome equation{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared
	test{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(opt)}}type of optimization{p_end}
{synopt:{cmd:e(which)}}{cmd:max} or {cmd:min}; whether optimizer is to perform maximization or minimization{p_end}
{synopt:{cmd:e(ml_method)}}type of {cmd:ml} method{p_end}
{synopt:{cmd:e(user)}}name of likelihood-evaluator program{p_end}
{synopt:{cmd:e(technique)}}maximization technique{p_end}
{synopt:{cmd:e(crittype)}}optimization criterion{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(gradient)}}gradient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{title:Author}

{pstd}
Takuya Hasebe{break}
Sophia University{break}
Tokyo, Japan{break}
thasebe@sophia.ac.jp


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 20, number 3: {browse "https://doi.org/10.1177/1536867X20953573":st0612}{p_end}
