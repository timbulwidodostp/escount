{smcl}
{* *! version 1.1.0 Aug2019}{...}
{cmd:help teescount}{right: ({browse "https://doi.org/10.1177/1536867X20953573":SJ20-3: st0612})}
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{cmd:teescount} {hline 2}}Treatment-effect estimator based on
endogenous switching count-data regression model{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 17 2}
{cmd:teescount} [{cmd:,} {it:options}]

{synoptset 28}{...}
{synopthdr}
{synoptline}
{synopt:{opt ate}}average treatment effect{p_end}
{synopt:{opt att}}average treatment effect on the treated{p_end}
{synopt:{opt atu}}average treatment effect on the untreated{p_end}
{synopt:{opt late(varname lower upper)}}local average treatment effect{p_end}
{synopt:{opt mte(nu)}}marginal treatment effect{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:teescount} estimates the treatment effects using the estimated parameters
of the switching regression model with count-data outcomes.  This command can
be executed only after the execution of {cmd:escount}.


{title:Options}

{phang}
{cmd:ate} estimates average treatment effect.  This is the default.

{phang}
{cmd:att} estimates average treatment effect on the treated.

{phang}
{cmd:atu} estimates average treatment effect on the untreated.

{phang}
{cmd:late(}{varname} {it:lower} {it:upper}{cmd:)} estimates the local average
treatment effect when the value of {it:varname} changes from {it:lower} to
{it:upper}.  {it:lower} and {it:upper} should be numerical.

{phang}
{cmd:mte(}{it:nu}{cmd:)} estimates marginal treatment effect evaluated at
{it:nu}.  {it:nu} should be numerical.


{title:Stored results}

{pstd}
{cmd:teescount} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}table of result{p_end}


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
