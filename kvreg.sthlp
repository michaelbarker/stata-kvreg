{smcl}
{* *! version 1.0  01dec2012}{...}
{cmd:help kvreg}
{hline}


{title:Title}

{phang}
{bf:kvreg} {hline 2} Covariance-dependent control function estimator from Klein Vella, 2010


{title:Syntax}

{p 8 17 2}
{cmd:kvreg}
{depvar:1} {depvar:2} {indepvars}
{ifin}
[{cmd:,} {it:options}] 

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt detail}}display intermediate estimates; reported standard errors are not correct{p_end}
{synopt:{opt query}}display moptimize settings after options have been parsed{p_end}
{synopt:{opt nogls}}skip GLS step for secondary equation{p_end}
{synopt:{opt tr:im(#,#)}}lower and upper trimming centiles. Default is trim(2,98){p_end}
{synopt:{opt single}}optimize the single-index portion of the value function only{p_end}
{synopt:{opt double}}optimize the double-index portion of the value function only{p_end}
{synopt:{opt joint}}optimize both portions of the value function jointly; default{p_end}

{syntab:Minimization}
{synopt:{cmdab:prim:init(}{it:{help kvreg##init_options:init_options}}{cmd:)}}specify moptimize_init options for estimation of primary equation{p_end}
{synopt:{cmdab:sec:init(}{it:{help kvreg##init_options:init_options}}{cmd:)}}specify moptimize_init options for estimation of secondary equation{p_end}
{synopt:{it:{help kvreg##minimize_options:min_options}}}control the minimization process{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:kvreg} fits a triangular system of two simultaneous linear equations. Models of this form
are generally analyzed with two-stage least squares or IV methods, which require one 
or more exclusion restriction. The KV estimator does not require an exclusion restriction, 
so the same set of independent variables may appear in both equations. To account for endogeneity
in the triangular system, the estimator constructs a control function using information 
from the conditional distribution of the error terms.  

{pstd}
The error structure must satisfy two requirements. 
First, one or both error terms must be heteroskedastic; ther variance should 
depend on X. Second, the correlation coefficient between the homoskedastic 
components of each error term must be constant and independent of X.  

{pstd}
Following the language in KV 2010, the equation with the endogenous regressor is termed 
the primary equation. In the primary equation, {depvar:1} is a linear function 
of {depvar:2} and {indepvars}. In the secondary equation, {depvar:2} is a linear 
function of {indepvars}. 


{title:Options}

{dlgtab:Main}

{phang}
{opt detail} display estimation results from each step of the estimation procedure. 
Intermediate estimation steps treat estimated components as true values, so reported
standard errors from these steps are not correct. 

{phang} {opt query} reports initial moptimize parameters for each minimization problem. 

{phang}
{opt nogls} skips the optional gls estimation of the secondary equation and subsequent 
re-estimation of conditional variance. 

{phang}
{opt tr:im(#,#)} sets the percentile bounds that identify boundary observations. These
observations are excluded from the value function in each minimization problem. Trimmed
observations are included in every other aspect of the estimation procedure, including  
the estimation of conditional expectation for non-trimmed observation. Users who wish to 
completely exclude certain observations should use {ifin} statements. 

{phang}
{opt single} estimates the primary equation with the single-index component of the
value function only. The single index value function reduces the computational 
requirements of the estimator, and performs well in simulations. Consistency of 
the single-index estimator has not been proven, however, so this option should be
used for exploratory purposes only. 

{phang}
{opt double} estimates the primary equation with the double-index component of the
value function only. Results are not proven consistent, and the computational 
requirements are not significantly reduced from joint value function. 

{phang}
{opt joint} minimizes the single and double-index components of the value function
jointly. This is the default.   

{marker init_options}{...}
{marker minimize_options}{...}
{dlgtab:Minimization}

{phang}
{opt secinit(init_options)} passes {help mf_moptimize##syn_step2:moptimize_init} options directly to moptimize problem for the secondary equation. 

{phang}
{opt priminit(init_options)} passes {help mf_moptimize##syn_step2:moptimize_init} options directly to moptimize problem for the primary equation.

{pmore}
Passing {help mf_moptimize##syn_step2:moptimize_init} options directly provides total control over the optimization procedure. 
These options take precedence over all {it:minimize_options}, and apply only to the secondary or primary equation problem, respectively. 
These options may be used to specify initial values, add linear coefficient constraints, or specify simplex deltas for the Nelder-Mead technique.

{phang}
{it:init_options} may contain one or more {help mf_moptimize##syn_step2:moptimize_init} options, separated by spaces or commas.  
In each option, the name of the moptimize problem, along with the subsequent comma, should be omitted. 
Similarly, the characters "moptimize_init_" should be omitted from each {it:init_option}. 
For example, the moptimize option, {...} 
{cmd:moptimize_init_search(M,} {c -(}{cmd:"on"}|{cmd:"off"}{c )-}{cmd:)} {...}
would be specified as 
{cmd:search(}{c -(}{cmd:"on"}|{cmd:"off"}{c )-}{cmd:)}.

{phang}
{it:minimize_options}: 
{opt dif:ficult},
{opt tech:nique(algorithm_spec)},
{opt iter:ate(#)},
[{cmdab:no:}]{opt lo:g}, 
{opt tr:ace}, 
{opt grad:ient},
{opt showstep},
{opt hess:ian},
{opt showtol:erance},
{opt tol:erance(#)}, {opt ltol:erance(#)}, 
{opt nrtol:erance(#)}, {opt nonrtol:erance}.

{pmore}
Minimize options perform the same functions as maximize options in 
maximum likelihood estimation; see {manhelp maximize R} for full descriptions. 

{pmore}
Minimize options apply to estimation of the primary equation only.


{title:Remarks}

{pstd}
This program estimates conditional variance functions semi-parametrically. 
Index parameters are estimated within a non-parametric density calculation.

{pstd}
The final reported standard errors include corrections for the first-stage estimation
and the semiparametric nature of the estimator

{pstd}
Kernel density estimation is done using the Gaussian kernel only. Bandwidths are chosen
using the standard plug-in estimate for the Gaussian kernel. Pilot bandwidth parameters are 
updated for the first 5 iterations of the optimization process. After the 5th iteration,
bandwidth parameters are fixed for the remainder of the optimization.

{pstd}
First-stage estimation sometimes results in negative estimates for conditional variance. 
When this happens, the control function is missing, and that observation is dropped
from the estimation. These observations are not included in e(sample).


{title:Examples}

{phang}{cmd:. kvreg y1 y2 x1 x2}


{marker references}{...}
{title:References}

{p 0 4}Klein, R., & Vella, F. (2010). 
Estimating a class of triangular simultaneous equations models without exclusion restrictions. 
Journal of Econometrics, 154(2), 154–164. 

{p 0 4}Ichimura, H. (1993). 
Semiparametric least squares (SLS) and weighted SLS estimation of single-index models. 
Journal of Econometrics, 58(1-2), 71–120. 

{p 0 4}Newey, W. K., Hsieh, F., & Robins, J. M. (2004). 
Twicing Kernels and a Small Bias Property of Semiparametric Estimators. 
Econometrica, 72(3), 947–962. 

{p 0 4} Stuetzle, W., & Mittal, Y. (1979). 
Some Comments on the Asymptotic Behavior of Robust Smoothers. 
In T. Gasser & M. Rosenblatt (Eds.), 
Smoothing Techniques for Curve Estimation, Lecture Notes, 757 (pp. 191–195). 
New York: Springer-Verlag.


