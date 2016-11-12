{smcl}
{* *! version 1.0  01dec2012}{...}
{cmd:help sls}
{hline}


{title:Title}

{phang}
{bf:sls} {hline 2} Semiparametric Least Squares from Ichimura, 1993.


{title:Syntax}

{p 8 17 2}
{cmd:sls}
{depvar} {indepvars}
{ifin}
[{cmd:,} {it:options}] 

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt tr:im(#,#)}}lower and upper trimming centiles. Default is trim(2,98){p_end}
{synopt:{cmdab:init(}{it:{help sls##init_options:init_options}}{cmd:)}}specify moptimize_init options for minimization{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}
{cmd:sls} performs semi-parametric estimation as described in Ichimura's 1993 paper. 
Kernel density estimates are calculated using a twicing kernel, as described in 
Newey, Hsieh, and Robins, 2004. The bandwidth parameter is estimated using a two-step
procedure. First, the index parameters are estimated using a plug-in optimal bandwidth 
estimate. The bandwidth and index parameter estimates from this preliminary estimation are
used to construct bounds on the final bandwidth estimate. The final optimal bandwidth 
parameter is estimated simultaneously with the index parameters. The bandwidth parameter 
is chosen as the minimizer of the squared error objective function, 
as described in Hardle, Hall, & Ichimura, 1993. 

{pstd}
Because no functional form is assumed, parameters are not identified in location or scale. Only 
ratios of coefficients are identified. To achieve identification, a normalization assumption
is required. This estimator does not estimate a constant term and restricts the coefficient 
on the first independent variable to 1. At least two independent variables are required
to run the estimator. At least one independent variable must be continuous.


{title:Options}

{dlgtab:Main}

{phang}
{opt tr:im(#,#)} sets the percentile bounds that identify boundary observations. These
observations are excluded from the value function in each iteration of the minimization 
procedure. They are also excluded from variance estimates. Trimmed observations are included 
in the the estimation of conditional expectation for non-trimmed observation. Users who 
wish to completely exclude certain observations should use {ifin} statements. The default
trimming centiles are (2,98). Bounds of 0 or 100 may be specified if no trimming is desired
on the lower or upper tail, respectively. 

{marker init_options}{...}
{phang}
{opt init(init_options)} passes {help mf_moptimize##syn_step2:moptimize_init} options directly to moptimize problem.

{pmore}
Passing {help mf_moptimize##syn_step2:moptimize_init} options directly provides total control over the optimization procedure. 
These options take precedence over all pre-set minimization options.  
These options may be used to specify initial values, add linear coefficient constraints, or change the optimization technique. 
Because the bandwidth parameter is estimated simultaneously, the moptimize problem has two equations. The first equation 
estimates index parameters and the second estimates the bandwidth. Equation number should be specified with setting
equation specific {help mf_moptimize##syn_step2:moptimize_init} options. 

{pmore}
{it:init_options} may contain one or more {help mf_moptimize##syn_step2:moptimize_init} options, separated by spaces or commas.  
In each option, the name of the moptimize problem, along with the subsequent comma, should be omitted. 
Similarly, the characters "moptimize_init_" should be omitted from each {it:init_option}. 
For example, the moptimize option, {...} 
{cmd:moptimize_init_search(M,} {c -(}{cmd:"on"}|{cmd:"off"}{c )-}{cmd:)} {...}
would be specified as {cmd:search("on")} or {cmd:search("off")}. 


{title:Examples}

{phang}{cmd:. sls y1 x1 x2 x3}

{phang}{cmd:. sls y1 x1 x2 x3, init(tracelevel("step"), eq_coefs(1, (1,0.75)))}


{marker references}{...}
{title:References}

{p 0 4}Hardle, W., Hall, P., & Ichimura, H. (1993). 
Optimal Smoothing in Single-Index Models. 
The Annals of Statistics, 21(1), 157–178.

{p 0 4}Ichimura, H. (1993). 
Semiparametric least squares (SLS) and weighted SLS estimation of single-index models. 
Journal of Econometrics, 58(1-2), 71–120. 

{p 0 4}Newey, W. K., Hsieh, F., & Robins, J. M. (2004). 
Twicing Kernels and a Small Bias Property of Semiparametric Estimators. 
Econometrica, 72(3), 947–962. 


{title:Author}

{pstd} Michael Barker {p_end}
{pstd} Georgetown University {p_end}
{pstd} mdb96@georgetown.edu {p_end}


{title:Also see}

{pstd}
{help sls_postestimation:sls postestimation}



