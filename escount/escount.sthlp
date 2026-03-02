{smcl}
{* *! version 1.1.0 Aug2019}{...}
{cmd:help escount}{right: ({browse "https://doi.org/10.1177/1536867X20953573":SJ20-3: st0612})}
{hline}

{title:Title}

{p2colset 5 16 18 2}{...}
{p2col:{cmd:escount} {hline 2}}Maximum likelihood estimation of endogenous
switching count-data regression model{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 15 2}
{cmd:escount}
{depvar} [{cmd:=}] {indepvars}
{ifin} 
[{it:{help escount##weight:weight}}]{cmd:,}
{cmdab:sel:ect(}{it:depvar_s} {cmd:=} {it:varlist_s} [{cmd:,}
{cmdab:nocons:tant} {cmdab:off:set(}{it:varname_o}{cmd:)}]{cmd:)} 
[{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {cmdab:sel:ect()}}specify selection equation{p_end}
{synopt:{opt model(model)}}specify distribution of count-data outcomes{p_end}
{synopt:{opt nocons:tant}}suppress constant term{p_end}
{synopt:{opt exp:osure(varname_e)}}include ln({it:varname_e}) in model with coefficient constrained to 1{p_end}
{synopt:{opt off:set(varname_o)}}include {it:varname_o} in model with coefficient constrained to 1{p_end}
{synopt:{opt const:raints(constraints)}}apply specified linear constraints{p_end}
{synopt:{opt intp:oints(#)}}use {it:#} Gauss-Hermite quadrature points; default is {cmd:intpoints(24)}{p_end}
{synopt:{opt vce:(vcetype)}}{it:vcetype} may be {cmd:robust}, {cmd:cluster}
{it:clustvar}, {cmd:oim}, or {cmd:opg}{p_end}
{synopt:{help escount##maximize:{it:maximize_options}}}control the maximization process{p_end}
{synoptline}
{pstd}
* {cmd:select()} is required.  The full specification is{p_end}
{p 10 10 2}
{cmdab:sel:ect(}{it:depvar_s} {cmd:=} {it:varlist_s}[{cmd:,} {opt nocons:tant} {opt off:set(varname_o)}]{cmd:)}{p_end}
{marker weight}{...}
{p 4 6 2}
{cmd:aweight}s, {cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s are
allowed; see {help weight}.{p_end}


{title:Description}

{pstd}
{cmd:escount} estimates the parameters of a switching regression model with
count-data outcomes, where the possible outcome differs across two alternate
states.  The model assumes lognormal latent heterogeneity.


{title:Options}

{phang}
{cmd:select(}{depvar:_s} {cmd:=} {indepvars:_s}[{cmd:,} {opt noconstant}
{cmd:offset(}{it:{help varname:varname_o}}{cmd:)}]{cmd:)} specifies the
variables and options for the selection equation.  {cmd:select()} is required.
 
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
[{cmd:,}
{it:options}]


{title:Options for predict}

{phang}
{cmd:psel}, the default, calculates the probability of receiving the treatment
status d=1 for each observation.

{phang}
{cmd:xb0} calculates the linear prediction of the outcome variable of the
treatment status d=0 for each observation.

{phang}
{cmd:xb1} calculates the linear prediction of the outcome variable of the
treatment status d=1 for each observation.

{phang}
{cmd:xbsel} calculates the linear prediction for the selection equation for
each observation.

{phang}
{cmd:mu0} calculates the expected value of the potential-outcome variable of
the treatment status d=0 for each observation: E(y0) =
exp(xb0)*exp(0.5*sigma0^2).

{phang}
{cmd:mu1} calculates the expected value of the potential-outcome variable of
the treatment status d=1 for each observation: E(y1) =
exp(xb1)*exp(0.5*sigma1^2).

{phang}
{cmd:mu0d0} calculates the expected value of the potential-outcome variable of
the treatment status d=0 conditional on receiving the treatment status d=0 for
each observation:
E(y0|d=0)=exp(xb0)*exp(0.5*sigma0^2)(normal(-(rho0*sigma0+xbsel))/normal(-xbsel)).

{phang}
{cmd:mu0d1} calculates the expected value of the potential-outcome variable of
the treatment status d=0 conditional on receiving the treatment status d=1 for
each observation:
E(y0|d=1)=exp(xb0)*exp(0.5*sigma0^2)*(normal(rho0*sigma0+xbsel)/normal(zg)).

{phang}
{cmd:mu1d0} calculates the expected value of the potential-outcome variable of
the treatment status d=1 conditional on receiving the treatment status d=0 for
each observation:
E(y1|d=0)=exp(xb1)*exp(0.5*sigma1^2)*(normal(-(rho1*sigma1+xbsel))/normal(-xbsel)).

{phang}
{cmd:mu1d1} calculates the expected value of the potential-outcome variable of
the treatment status d=1 conditional on receiving the treatment status d=1 for
each observation:
E(y1|d=1)=exp(xb1)*exp(0.5*sigma1^2)*(normal(rho1*sigma1+xbsel)/normal(xbsel)).


{title:Stored results}

{pstd}
{cmd:escount} stores the following in {cmd:e()}:

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
{synopt:{cmd:e(p)}}significance of comparison test{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(chi2)}}chi-squared{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(n_quad)}}number of quadrature points{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:escount}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(wtype)}}weight type{p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(offset1)}}offset for outcome equation{p_end}
{synopt:{cmd:e(offset2)}}offset for selection equation{p_end}
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
