/*******************************************************************************
 Michael Barker: Sept 26, 2012
 Klein and Vella 2010 estimator

X   todo: save, name, and format estimation results from each step
X   todo: allow user specification of np-options (kernel, bandwidth, trim)
X   todo: allow optimization options to be specified for both equations or individually
X   todo: create simplified user options for a restricted set of commonly used.
          - see: help mlopts , help opts_exclusive
X   todo: stop estimating density at points that are going to be trimmed.  
          - no good way to do this with "leave one out estimator"  
          X will have to write a new density calculation fxn with data pts and eval pts.
O   todo: allow specification of different varlists for each step. Currently,
            it can be done directly in optimization options, but then any new
            variables are not considered in the marksample command
O   todo: research on the inclusion of discrete and binary data in this estimator
O   todo: depending on results of above, allow for factor variables 
X   todo: consider the inclusion of bootstrap option. Should it be included as 
          an option, or just let people do it themselves if they want to.
          - Do not include bootstrap option at this time. There are too many 
            suboptions and special cases for that command. 
O   todo: research se options for this estimator in general. Is there any place
          for robust or clustered standard errors? 
X   todo: choose a name and update function names - kvreg
X   todo: organize total project into a single package (one or two files).  
           X Use an mlib for mata files
O   todo: update multivariate kernel function.
X   todo: include default options for trimming percentile in syntax statement.
          - cannot do with numlist option, only int or real.
X   todo: review data types: make sure all variables stored in Stata are doubles
X   todo: change output for trimmed observation tabulation.
*******************************************************************************/


/*******************************************************************************
 Stata Function.   
 Input: Y1 Y2 X1 X2 additional X variables

*******************************************************************************/

program kvreg
    version 11.0

    if replay() {
        if ("`e(cmd)'" != "kvreg") error 301
        Replay `0'
    }
    else Estimate `0'
end

program Replay 
    syntax 
    display as text _n "Final Coefficient Estimates" _n "---------------------------"
    ereturn display
end

program Estimate, eclass
    version 11
# delimit ;
    syntax varlist(min=4 numeric) [if] [in] , 
        [Generate(name) noGLS detail peer single double joint startsingle workaround
         TRim(numlist min=2 max=2 >=0)
         Kernel(passthru) BWidth(passthru) EVALuator(passthru) 
         SECinit(string asis) PRIMinit(string asis)

    	 random repeat(integer 10)
         DIFficult TECHnique(string) 
         ITERate(numlist max=1 >0 integer) 
         TOLerance(numlist max=1 >0) LTOLerance(numlist max=1 >0) 
         NRTOLerance(numlist max=1 >0) NONRTOLerance 
         TRace GRADient showstep HESSian SHOWTOLerance noLOg 
         query]
		 ;
 
	local minoptions `random' repeat(`repeat') `difficult' technique(`technique') 
         iterate(`iterate') tolerance(`tolerance') ltolerance(`ltolerance')
         nrtolerance(`nrtolerance') `nonrtolerance' `trace' `gradient' `showstep'
		 `hessian' `showtolerance' `log' `query'
		;   

# delimit cr
    
    marksample touse
    
    _nobs `touse' [`weight'`exp'], min(2)
    local N = r(N)

/*******************************************************************************
    Verify Arguments, Create Tempvars, Pass Options    
*******************************************************************************/
    opts_exclusive "`single' `double' `joint'"
    local objfxn "`single'`double'`joint'"
    if "`objfxn'" == "" local objfxn joint

    * Fill in default pctiles 
    if ("`trim'"=="") local trim "2,98"

    * Define mopimize problem and parse user options 
	* Parse early to prevent error later in routine. 
	* Reparse options later so user can override default optins. 
	tempname M
	mata: `M' = moptimize_init()
	np_moptimize_options `M' , `minoptions' initoptions(`priminit')

   
/*******************************************************************************
    Parse Varlist
    Note: Conditional Variance estimates should not use the same offset
          variable for primary and secondary equation estimates. Reverse the
          order of the first two x variables for the secondary equation 
		  conditional variance estimate. 
          Use original order for primary equation index.
*******************************************************************************/

    gettoken y1 xvars : varlist
    gettoken y2 xvars : xvars 

    gettoken x1 xsls : xvars 
    gettoken x2 xsls : xsls
    local xsls "`x2' `x1' `xsls'"

/*******************************************************************************
    Estimate Secondary Regression
*******************************************************************************/

    tempvar yhat v vsq vhat vsqhat index ctrl sv tx
    * Estimate variance
    quietly: regress `y2' `xvars' if `touse'
    quietly: predict `v' if `touse' , resid
    quietly: gen double `vsq' = `v'^2
    
    * Store estimates
    matrix b1 = e(b)
    matrix coleq b1 = "`y2'"
    tempvar y2reg 
    ereturn local title "Linear Regression of Secondary Equation"
    _estimates hold `y2reg'

    * Estimate Conditional Variance
	sls `vsq' `xsls' if `touse' , trim(`trim') initoptions(iterid("Secondary") , eq_name(1, "Variance") , `secinit')

    if "`gls'" != "nogls" {
             
		quietly: predict `vsqhat' if `touse' , ey
        * Repeat variance estimate with GLS
        quietly: gen double `vhat' = sqrt(`vsqhat')
		quietly: replace `touse'=0 if missing(`vhat')
        drop `v' `vsq' `vsqhat' 

        quietly: vwls `y2' `xvars' if `touse' , sd(`vhat')
        quietly: predict `yhat' if `touse' , xb 
        quietly: gen double `v' = `y2'-`yhat' if `touse'
        quietly: gen double `vsq' = `v'^2 if `touse'

        * Store estimates
        matrix b1 = e(b)
        matrix coleq b1 = "`y2'"
        tempname y2gls 
        ereturn local title "GLS Estimation of Secondary Equation"
        _estimates hold `y2gls'

        * Re-Estimate Conditional Variance
		sls `vsq' `xsls' if `touse' , trim(`trim') initoptions(iterid("Secondary") , eq_name(1, "Variance") , `secinit') 
    }

	* Calculate predicted values from sls
	quietly: predict `index' if `touse' , xb
	quietly: predict `vsqhat' if `touse' , ey
	quietly: predict `tx' if `touse' , tr 

    * Store Conditional Variance Estimates
    matrix b2 = e(b) 
    matrix coleq b2 = "Index_2"
    tempvar cvar2 
    ereturn local title "SLS Estimation of Conditional Variance from Secondary Equation"
    _estimates hold `cvar2'


*** Save vstar - partial control function 
	quietly: gen `sv' = sqrt(`vsqhat') if `touse'
    gen double `ctrl' = `v' / `sv' if `touse'
	* Drop observations that have negative conditional variance estimates.
	replace `touse'=0 if missing(`ctrl')


/*******************************************************************************
    Estimate Primary Regression
*******************************************************************************/

    *** Initialize Optimization Problem and Set Default Options
    mata: moptimize_init_evaluator(`M', &kv1`objfxn'fxn())
	mata: moptimize_init_evaluatortype(`M', "gf1")
	mata: moptimize_init_which(`M', "min")
    mata: moptimize_init_iterid(`M', "Primary")
    mata: moptimize_init_valueid(`M', "SSq(b)")
	mata: moptimize_init_touse(`M', "`touse'")

    * Dependent and non-parameter variables
    mata: moptimize_init_depvar(`M', 1, "`y1'" )
    mata: moptimize_init_depvar(`M', 2, "`tx'" )
    mata: moptimize_init_depvar(`M', 3, "`index'" )
	* Additional variables for variance calculation
    mata: moptimize_init_depvar(`M', 4, "`v'" )
    mata: moptimize_init_depvar(`M', 5, "`sv'" )
    
    * Primary Linear Equation
    mata: moptimize_init_eq_name(`M', 1, "`y1'")
    mata: moptimize_init_eq_indepvars(`M', 1, "`y2' `xvars'")

    * Control Function Term 
    mata: moptimize_init_eq_cons(`M', 2, "off")
    mata: moptimize_init_eq_indepvars(`M', 2, "`ctrl'")
    mata: moptimize_init_eq_name(`M', 2, "Ctrl_Fxn")
    mata: moptimize_init_eq_colnames(`M', 2, ("rho"))

    * Conditional Variance Equation
    gettoken x1 X2 : xvars 
    mata: moptimize_init_eq_name(`M', 3, "Index_1")
    mata: moptimize_init_eq_indepvars(`M', 3, "`X2'")
    mata: moptimize_init_eq_cons(`M', 3, "off")
    mata: moptimize_init_eq_offset(`M', 3, "`x1'")

    * Optimization Control Options
	* Reparse options in case user options were changed above
    np_moptimize_options `M' , `minoptions' initoptions(`priminit')
    
	*** Run Optimization Problem from Mata
    mata: kv1(`M')   
	* Calculate covariance matrix
	mata: kv1_cov(`M')

*** save estimation results
	tempname b V
    matrix `b' = e(b)
	matrix `V' = e(Vkvreg)	

	* Save scalar return values
	tempname SSE iterations converged 
	scalar `SSE' = e(SSE)
	scalar `iterations' = e(iterations)
	scalar `converged' = e(converged)

    tempvar kvreg
    ereturn local title "Simultaneous Estimation of Primary Equation"
    _estimates hold `kvreg'

*** Display estimates from each stage 

    if "`detail'" == "detail" {
        foreach estname in `y2reg' `y2gls' `cvar2' `kvreg' { 
        _estimates unhold `estname' 
        display as txt _n _n
        display as txt "`e(title)'"
        display as txt "-------------------------------------------------------------"
        ereturn display
        }
    }
	
	local colnames:colfullnames `b'
	matrix colnames `V' = `colnames'
	matrix rownames `V' = `colnames'

    quietly: count if `touse'
    ereturn post `b' `V' , obs(`r(N)') esample(`touse') depname("`y1'")

	
	* Additional return values
	* locals
	ereturn local cmd "kvreg"
	ereturn local cmdline "kvreg `0'"
	ereturn local depvar "`y1'"
	ereturn local indepvars = ltrim("`y2' `xvars'")
	ereturn local properties "b V"

	* scalars
	ereturn scalar SSE = `SSE' 
	ereturn scalar iterations = `iterations' 
	ereturn scalar converged = `converged' 

	* matrices

	* Add return values

    Replay    
end

/*******************************************************************************
 Begin Mata fxn definitions
*******************************************************************************/

version 11 
mata

/*******************************************************************************
 Estimate Primary Equation
*******************************************************************************/

void kv1(transmorphic scalar M) { 

	struct cexp_parameters scalar CE
    CE = cexp_define()
	CE.kernel = &kernel_gaussian()
	CE.dkernel = &dkernel_gaussian()
    moptimize_init_userinfo(M , 1, CE)

    string scalar touse
    touse = moptimize_init_touse(M)

    kvpilot(M)
    // db moptimize_query(M)
    moptimize(M)
    moptimize_result_post(M) 

}

void kv1_cov(transmorphic scalar M) {
	/*******************************************************************************
 	Calculate variance estimate 
	Component Definitions:
	Z = (y2,X,1)
	Ztheta = Z*theta  
	X = (X,1)
	X2 = (x2,x3,...,xk)

	alpha = (theta , rho, bu)

	*******************************************************************************/

	real scalar rho
	real colvector y1, Ztheta, tx, v, Sv, I2, I1, usq, Ssq, Su1, Su2, M1, M2, dphi 

	real matrix w1, w2, Z, X, X2, EX2, dSu1, dSu2, Gab, Gc, G

    // Components from primary equation (second stage estimation)
    // dependent var
    y1 = moptimize_util_depvar(M,1) 
	N = rows(y1)

	// Z = (y2,X,1) : variables in linear component of primary equation 
	Z = (moptimize_init_eq_indepvars(M, 1) , J(N,1,1))
	// X = (X,1) : variables in linear component of secondary equation
	X  = Z[.,(2..cols(Z))]
	y2 = Z[.,1]

    // Linear component of primary equation
    Ztheta = moptimize_util_xb(M,moptimize_result_coefs(M),1)

	// Indicator trimming vector
	tx = moptimize_util_depvar(M,2) 
	// Trimming vector of all ones, for no trimming 
	tx0 = J(N,1,1)

	// Parameters for conditional expectation
    struct cexp_parameters scalar CE 
    CE = moptimize_util_userinfo(M,1)

	// Control Function components
	// Su1, Su2, v, Sv, I1, I2
	// rho: control function coefficient
	rho = moptimize_result_eq_coefs(M,2)
	// v: estimated residual from secondary equation 
	v  = moptimize_init_depvar(M, 4) 	
	Sv = moptimize_init_depvar(M, 5) 	
	// Conditional Variance index from secondary equation (first stage estimation)
	I2 = moptimize_util_depvar(M,3) 
    // Conditional Variance Index 
    I1 = moptimize_util_xb(M,moptimize_result_coefs(M),3)
    // Square of uhat at estimated parameter values
    usq = (y1-Ztheta):^2

	// Bandwidth
	h1 = (*CE.bwidth)(usq,I1,CE,tx)
	h2 = (*CE.bwidth)(usq,(I1,I2),CE,tx)

	// Estimate variance, conditional on one index
    Ssq = cexp(usq,I1,CE,tx0,h1)
    // Ssq = Ssq :* strim(Ssq)
    Su1 = sqrt(Ssq) 

    // Estimate variance, conditional on two indices
	Ssq = cexp(usq,(I1,I2),CE,tx0,h2)
    // Ssq = Ssq :* strim(Ssq)
    Su2 = sqrt(Ssq) 
	

	// 1. Calculate Gradient Weights
	// W1 = dM1dalpha and W2=dM2dalpha ; alpha = (Theta, rho, bu) 

	// M1 and M2: predicted Y for single and double index control functions.
	M1 = Ztheta + (rho:*(Su1:/Sv):*v)
	M2 = Ztheta + (rho:*(Su2:/Sv):*v)
	
	// Residuals for each predicted Y
	R1 = y1 - M1
	R2 = y1 - M2

	// dM/dtheta : linear coefficients
	// Theta = (Th_y2, Th_X, Th_cons)
	dM1dTheta = (y2-(rho*(v:/(Sv:*Su1)):*cexp(y2:*R1,I1,CE,tx,h1)) , X)
	dM2dTheta = (y2-(rho*(v:/(Sv:*Su2)):*cexp(y2:*R2,(I1,I2),CE,tx,h2)) , X)

	// dM/dbu : conditional variance coefficients
	// bu has no coefficient for x1 or constant term 
	dSu1  = dcexpdb_vec(usq, I1, X[.,(2..cols(X)-1)] , CE, tx, h1) 
	dM1dBu = 0.5*rho*(v:/(Sv:*Su1)) :* dSu1

	dSu2  = dcexpdb_vec(usq, I1, X[.,(2..cols(X)-1)] , CE, tx, h2, I2) 
	dM2dBu = 0.5*rho*(v:/(Sv:*Su2)) :* dSu2

	W1 = tx :* (dM1dTheta , (v:*Su1:*Sv) , dM1dBu)
	W2 = tx :* (dM2dTheta , (v:*Su2:*Sv) , dM2dBu)


	// 2. GradA 
	// GradA Component
	Ga = -tx :* (R1:*W1 + R2:*W2) 


	// 3. GradB 

	// au: derivative of control fxn w.r.t. Su1sq
	// av: derivative of control fxn w.r.t. Svsq
	
	// Single Index Components
	au = 0.5*rho*(v:/(Sv:*Su1)) 
	av = 0.5*rho*(v:*Su1:/(Sv:^3))
	Gb1 = (usq - Su1:^2) :* cexp(tx:*W1:*au, I1, CE, tx, h1) - (v:^2 - Sv:^2) :* cexp(W1:*av, I1, CE, tx, h1)

	// Double Index Components
	au = 0.5*rho*(v:/(Sv:*Su2)) 
	av = 0.5*rho*(v:*Su2:/(Sv:^3))
	Gb2 = (usq - Su2:^2) :* cexp(tx:*W2:*au, (I1,I2), CE, tx, h2) - (v:^2 - Sv:^2) :* cexp(W2:*av, (I1,I2), CE, tx, h2)
	
	// GradB Component
	Gb = tx :* (Gb1 + Gb2)


	// 4. GradC 
	// Correct gradient for uncertainty about first stage coefficient estimates, eta 
	// Two parts: dG/deta , deta 	

	// 4.1 deta
	// Construct variance components from 1st stage estimation.

	// GLS terms
	vsq	= v:^2
	Svsq= Sv:^2 
	
	Hgls = cross(X,(tx:/Svsq),X) 
	// NxK+1 * K+1xK+1 -> NxK+1
	Ggls  = (tx :* (X) :* v :/ Svsq) * invsym(Hgls) 

	// SLS terms	
	// Swap 1st and 2nd X variables 
	Xsls = (X[.,2],X[.,1],X[.,(3..cols(X))])
	// Remove first X and constant 
	Xsls = Xsls[.,(2..cols(X)-1)]

	dSvsq= dcexpdb_vec(vsq, I2, Xsls, CE, tx, h1)
	Hsls = cross(dSvsq,dSvsq) 
	// NxK-1 * K-1xK-1 -> NxK-1
	Gsls = ((vsq-Svsq) :* (dSvsq - cexp(dSvsq,I2,CE,tx,h1))) * invsym(Hsls)

	// (NxK+1 , NxK-1) -> (N, 2K)
	// Collect GLS and SLS terms
	Geta = (Ggls,Gsls)
	
	// 4.2 dG/deta
	dR1gls = rho*(Su1:/Sv):*X
	dR2gls = rho*(Su2:/Sv):*X

	dR1sls = 0.5*rho* v:*(Su1:/Sv:^3) :* dSvsq
	dR2sls = 0.5*rho* v:*(Su2:/Sv:^3) :* dSvsq - 
				0.5*rho*(v:/(Sv:*Su2)) :* dcexpdb_vec(usq,I2,Xsls,CE,tx,h2,I1)

	dR1 = (dR1gls , dR1sls)
	dR2 = (dR2gls , dR2sls) 

	// Hc = dG/deta
	// K1xN * N*K2 -> K1xK2
	Hc = (cross(W1,tx,dR1) + cross(W2,tx,dR2)) 
	// Final Gc component
	// NxK2 * K2xK1 -> NxK1
	Gc = Geta*Hc' 

	// 5. Final Covariance Matrix Calculation
	// NxK1
	G = tx:*(Ga + Gb + Gc)
	Hinv = invsym((cross(W1,W1) + cross(W2,W2))) 

	// Klein demeans gradient. I'm not sure why, but it does not make much difference
	meanG = mean(G)	
	Var = Hinv * (crossdev(G,meanG,G,meanG)) * Hinv 

	// Return Values
	st_matrix("e(Vkvreg)" , Var)

	st_numscalar("e(SSE)" , moptimize_result_value(M))
	st_numscalar("e(iterations)" , moptimize_result_iterations(M))
	st_numscalar("e(converged)" , moptimize_result_converged(M))

}




/*******************************************************************************
Objective Function with Simultaneous Bandwidth Estimation
*******************************************************************************/

function kv1jointfxn(transmorphic M, real scalar todo, real rowvector b, 
						real colvector fv, real matrix S, real matrix H) {

    // Declare data variables 
    real colvector y, tx, I2, ZB, I1, C
    struct cexp_parameters scalar CE 
    
    // Dependent and non-parameter variables 
    y  = moptimize_util_depvar(M,1) 
    // Indicator trimming vector
    tx = moptimize_util_depvar(M,2) 
    // Index from secondary equation conditional variance estimation 
    I2 = moptimize_util_depvar(M,3) 

    // Linear Equation: Y2*gamma + X*beta 
    ZB = moptimize_util_xb(M,b,1)
    // Control Function
    C  = moptimize_util_xb(M,b,2)
   // Conditional Variance Index 
    I1  = moptimize_util_xb(M,b,3)

    // Conditional expectation parameters
    CE = moptimize_util_userinfo(M,1)

	// Bandwidth parameter
    h1 = moptimize_util_userinfo(M,2)
    h2 = moptimize_util_userinfo(M,3)

    // Residual at current parameter values
	u = y-ZB
    usq = u:^2

    // Estimate variance, conditional on one index
    Ssq1 = cexp(usq,I1,CE,tx,h1)
    Ssq1 = Ssq1 :* strim(Ssq1)
    S1 = sqrt(Ssq1) 
    
    /*** 
        For Large Negative Numbers, the trimming function returns system missing. 
        The return value is too small to be represented in a double. 
        In this case, strim() returns -smallestdouble().
        See fxn strim() in lnonparam.do
    ***/

    // Estimate variance, conditional on two indeces
    Ssq2 = cexp(usq,(I1,I2),CE,tx,h2)
    Ssq2 = Ssq2 :* strim(Ssq2)
    S2 = sqrt(Ssq2) 

    // Calculate Objective Function 
	diff1 = u-C:*S1
	diff2 = u-C:*S2
    Q1 = (diff1):^2  
    Q2 = (diff2):^2  
    Q  = (Q1 + Q2)

    // Check for missing caused by zero in denominator of nonparametric  density estimation.
    trim = missing(Q)
    if (trim>0) {
        printf("%f observations trimmed due to zero denominator in cond. exp.\n", trim)
    }
    _editmissing(Q, 0)
    
    // Return objective fxn
    fv = tx :* Q

	if (todo==1) {
		// S = -2 :* diff :* Z 			
		X =  moptimize_init_eq_indepvars(M, 3)
		Z = (moptimize_init_eq_indepvars(M, 1) , J(rows(X),1,1)) 

		// Look for missing values in conditional expectation derivative
		// that were not caused by trimming
		dctrl1 = C:/S1
		miss = sum(tx)-nonmissing(dctrl1)
    	if (miss>0) {
			printf("%f observations trimmed due to zero denominator of derivative\n", miss)
    	}
		_editmissing(dctrl1, 0)

		dQ1dB		= Z   - (dctrl1 :* cexp(Z:*u,I1,CE,tx,h1))
		dQ1drho		= S1  :* moptimize_init_eq_indepvars(M, 2)
		dQ1ddelta	= 0.5 :* (dctrl1) :* dcexpdb_vec(usq , I1 , X , CE , tx , h1) 

		S1 = -2 :* tx:* diff1 :* (dQ1dB,dQ1drho,dQ1ddelta)

		dctrl2 = C:/S2
		miss = sum(tx)-nonmissing(dctrl2)
    	if (miss>0) {
			printf("%f observations trimmed due to zero denominator of derivative\n", miss)
    	}
		_editmissing(dctrl2, 0)

		dQ2dB		= Z   - (dctrl2 :* cexp(Z:*u,(I1,I2),CE,tx,h2))
		dQ2drho		= S2  :* moptimize_init_eq_indepvars(M, 2)
		dQ2ddelta	= 0.5 :* (dctrl2) :* dcexpdb_vec(usq , I1 , X , CE , tx , h2, I2) 

		S2 = -2 :* tx:* diff2 :* (dQ2dB,dQ2drho,dQ2ddelta)

		S = S1 + S2

	}			

}

function kv1fxn_pilot(transmorphic M, real scalar todo, real rowvector b, 
						real colvector fv, real matrix S, real matrix H) {

    // Declare data variables 
    real colvector y, tx, I2, ZB, I1, C
    struct cexp_parameters scalar CE 
    
    // Dependent and non-parameter variables 
    y  = moptimize_util_depvar(M,1) 
    // Indicator trimming vector
    tx = moptimize_util_depvar(M,2) 
    // Index from secondary equation conditional variance estimation 
    I2 = moptimize_util_depvar(M,3) 

    // Linear Equation: Y2*gamma + X*beta 
    ZB = moptimize_util_xb(M,b,1)
    // Control Function
    C  = moptimize_util_xb(M,b,2)
   // Conditional Variance Index 
    I1  = moptimize_util_xb(M,b,3)

    // Conditional expectation parameters
    CE = moptimize_util_userinfo(M,1)
   
    // Residual at current parameter values
	u = y-ZB
    usq = u:^2

    // Estimate variance, conditional on one index
    Ssq1 = cexp(usq,I1,CE,tx)
    // Ssq1 = Ssq1 :* strim(Ssq1)
    S1 = sqrt(Ssq1) 
    
    /*** 
        For Large Negative Numbers, the trimming function returns system missing. 
        The return value is too small to be represented in a double. 
        In this case, strim() returns -smallestdouble().
        See fxn strim() in kvutil.do
    ***/

    // Estimate variance, conditional on two indeces
    Ssq2 = cexp(usq,(I1,I2),CE,tx)
    /// Ssq2 = Ssq2 :* strim(Ssq2)
    S2 = sqrt(Ssq2) 

    // Calculate Objective Function 
	diff1 = u-C:*S1
	diff2 = u-C:*S2
    Q1 = (diff1):^2  
    Q2 = (diff2):^2  
    Q  = (Q1 + Q2)

    // Check for missing caused by zero in denominator of nonparametric  density estimation.
    trim = missing(Q)
    if (trim>0) {
        printf("%f observations trimmed due to zero denominator in cond. exp.\n", trim)
    }
    _editmissing(Q, 0)
    
    // Return objective fxn
    fv = tx :* Q

	if (todo==1) {

		h1 = (*CE.bwidth)(usq,I1,CE,tx)
		h2 = (*CE.bwidth)(usq,(I1,I2),CE,tx)

		// S = -2 :* diff :* Z 			
		X =  moptimize_init_eq_indepvars(M, 3)
		Z = (moptimize_init_eq_indepvars(M, 1) , J(rows(X),1,1)) 

		// Look for missing values in conditional expectation derivative
		// that were not caused by trimming
		dctrl1 = C:/S1
		miss = sum(tx)-nonmissing(dctrl1)
    	if (miss>0) {
			printf("%f observations trimmed due to zero denominator of derivative\n", miss)
    	}
		_editmissing(dctrl1, 0)

		dQ1dB		= Z   - (dctrl1 :* cexp(Z:*u,I1,CE,tx,h1))
		dQ1drho		= S1  :* moptimize_init_eq_indepvars(M, 2)
		dQ1ddelta	= 0.5 :* (dctrl1) :* dcexpdb_vec(usq , I1 , X , CE , tx , h1) 

		S1 = -2 :* tx:* diff1 :* (dQ1dB,dQ1drho,dQ1ddelta)

		dctrl2 = C:/S2
		miss = sum(tx)-nonmissing(dctrl2)
    	if (miss>0) {
			printf("%f observations trimmed due to zero denominator of derivative\n", miss)
    	}
		_editmissing(dctrl2, 0)

		dQ2dB		= Z   - (dctrl2 :* cexp(Z:*u,(I1,I2),CE,tx,h2))
		dQ2drho		= S2  :* moptimize_init_eq_indepvars(M, 2)
		dQ2ddelta	= 0.5 :* (dctrl2) :* dcexpdb_vec(usq , I1 , X , CE , tx , h2, I2) 

		S2 = -2 :* tx:* diff2 :* (dQ2dB,dQ2drho,dQ2ddelta)

		S = S1 + S2

	}			

}

/*******************************************************************************
 kvpilot: Estimate starting values with pilot bandwidth. 
*******************************************************************************/

void kvpilot(transmorphic scalar M) {

    // Make a copy of current moptimize problem to estimate starting values 
    transmorphic scalar Mp
    Mp = M


    moptimize_init_evaluator(Mp, &kv1fxn_pilot())
	moptimize_init_evaluatortype(Mp, "gf1")
    moptimize_init_iterid(Mp, "Pilot")
    moptimize_init_conv_maxiter(Mp, 5) 
	moptimize_init_conv_warning(Mp, "off")
    moptimize(Mp)


	// Compute pilot bandwidth estimates and save in M
    struct cexp_parameters scalar CE 
    CE = moptimize_util_userinfo(Mp,1)
	
	b 	= moptimize_result_coefs(Mp)

    usq = (moptimize_util_depvar(Mp,1)-moptimize_util_xb(Mp,b,1)):^2 
    I1  = moptimize_util_xb(Mp,b,3)
    I2  = moptimize_util_depvar(Mp,3) 
    tx  = moptimize_util_depvar(Mp,2) 

	h1 = (*CE.bwidth)(usq,I1,CE,tx)
	h2 = (*CE.bwidth)(usq,(I1,I2),CE,tx)
    moptimize_init_userinfo(M , 2, h1)
    moptimize_init_userinfo(M , 3, h2)

	"Bandwidth Values:"
	"single index: " 
	h1
	"double index: " 
	h2

    // Store starting values as initial values in original problem, M.
    k = moptimize_init_eq_n(Mp)
    for (j=1 ; j<=k ; j++) {
        moptimize_init_eq_coefs(M, j, moptimize_result_eq_coefs(Mp, j))
    }

    // Do not search for new initial values.
    moptimize_init_search(M, "off")
}


/*******************************************************************************
 Objective Functions for single index, double index, and joint minimization
 problems. Default is joint.
*******************************************************************************/

function kv1singlefxn(transmorphic M, real scalar todo, real rowvector b, 
						real colvector fv, real matrix S, real matrix H) {

    // Declare data variables 
    real colvector y, tx, I2, ZB, I, C
    struct cexp_parameters scalar CE 
    
    // Dependent and non-parameter variables 
    y  = moptimize_util_depvar(M,1) 
    // Indicator trimming vector
    tx = moptimize_util_depvar(M,2) 
    // Index from secondary equation conditional variance estimation 
    I2 = moptimize_util_depvar(M,3) 

    // Linear Equation: Y2*gamma + X*beta 
    ZB = moptimize_util_xb(M,b,1)
    // Control Function
    C  = moptimize_util_xb(M,b,2)
    // Conditional Variance Index 
    I  = moptimize_util_xb(M,b,3)

    // Conditional expectation parameters
    CE = moptimize_util_userinfo(M,1)
   
    // Residual at current parameter values
	u = y-ZB
    usq = u:^2

	// Bandwidth
    h = moptimize_util_userinfo(M,2)

    // Estimate variance, conditional on one index
    Ssq1 = cexp(usq,I,CE,tx,h)
    // Ssq1 = Ssq1 :* strim(Ssq1)
	// Smooth trimming throws off analytical derivative.
	// Use Gaussian kernel to avoid negative predictions and don't trim.
	// Consider estimating a log or exponential transformation to avoid negative values.
    S1 = sqrt(Ssq1) 

    // Calculate Objective Function 
	// diff = (y - (ZB + C:*S1))
	diff1 = (u - C:*S1)
    Q1 = (diff1):^2  

    // Check for missing caused by zero in denominator of nonparametric  density estimation.
    trim = missing(Q1)
    if (trim>0) {
        printf("%f observations trimmed due to zero denominator in cond. exp.\n", trim)
    }
    _editmissing(Q1, 0)
    
    // Return objective fxn
    fv = tx :* Q1

	if (todo==1) {
		// S = -2 :* diff :* Z 			
		X    	= moptimize_init_eq_indepvars(M, 3)
		Z 		= (moptimize_init_eq_indepvars(M, 1) , J(rows(X),1,1)) 
    	// h = CE.optbwidth 

		// Look for missing values in conditional expectation derivative
		// that were not caused by trimming
		dctrl = C:/S1
		miss = sum(tx)-nonmissing(dctrl)
    	if (miss>0) {
			printf("%f observations trimmed due to zero denominator of derivative\n", miss)
    	}
		_editmissing(dctrl, 0)

		dQdB		= Z   - (dctrl :* cexp(Z:*u,I,CE,tx,h))
		dQdrho		= S1  :* moptimize_init_eq_indepvars(M, 2)
		dQddelta	= 0.5 :* (dctrl) :* dcexpdb_vec(usq , I , X , CE , tx , h) 

		S = -2 :* tx:* diff1 :* (dQdB,dQdrho,dQddelta)
	}			
}


function kv1doublefxn(transmorphic M, real scalar todo, real rowvector b, 
						real colvector fv, real matrix S, real matrix H) {

    // Declare data variables 
    real colvector y, tx, I2, ZB, I, C
    struct cexp_parameters scalar CE 
    
    // Dependent and non-parameter variables 
    y  = moptimize_util_depvar(M,1) 
    // Indicator trimming vector
    tx = moptimize_util_depvar(M,2) 
    // Index from secondary equation conditional variance estimation 
    I2 = moptimize_util_depvar(M,3) 

    // Linear Equation: Y2*gamma + X*beta 
    ZB = moptimize_util_xb(M,b,1)
    // Control Function
    C  = moptimize_util_xb(M,b,2)
    // Conditional Variance Index 
    I  = moptimize_util_xb(M,b,3)

    // Conditional expectation parameters
    CE = moptimize_util_userinfo(M,1)
   

    // Residual at current parameter values
	u = y-ZB
    usq = u:^2

	// Bandwidth
    h = moptimize_util_userinfo(M,3)

	// Estimate variance, conditional on two indeces
    Ssq2 = cexp(usq,(I,I2),CE,tx,h)
    // Ssq2 = Ssq2 :* strim(Ssq2)
    S2 = sqrt(Ssq2) 

    // Calculate Objective Function 
	diff2 = (u - C:*S2)
    Q2 = (diff2):^2  

    // Check for missing caused by zero in denominator of nonparametric  density estimation.
    trim = missing(Q2)
    if (trim>0) {
        printf("%f observations trimmed due to zero denominator in cond. exp.\n", trim)
    }
    _editmissing(Q2, 0)
    
    // Return objective fx
    fv = tx :* Q2


	if (todo==1) {
		// S = -2 :* diff :* Z 			
		X    	= moptimize_init_eq_indepvars(M, 3)
		Z 		=(moptimize_init_eq_indepvars(M, 1) , J(rows(X),1,1)) 
    	// h = CE.optbwidth 

		// Look for missing values in conditional expectation derivative
		// that were not caused by trimming
		dctrl2 = C:/S2
		miss = sum(tx)-nonmissing(dctrl2)
    	if (miss>0) {
			printf("%f observations trimmed due to zero denominator of derivative\n", miss)
    	}
		_editmissing(dctrl2, 0)

		dQ2dB		= Z   - (dctrl2 :* cexp(Z:*u,(I,I2),CE,tx,h))
		dQ2drho		= S2  :* moptimize_init_eq_indepvars(M, 2)
		dQ2ddelta	= 0.5 :* (dctrl2) :* dcexpdb_vec(usq , I , X , CE , tx , h, I2) 

		S = -2 :* tx:* diff2 :* (dQ2dB,dQ2drho,dQ2ddelta)
	}			


}


void d(X) {
	"dimensions"
	rows(X)
	cols(X)
}

void compare(real matrix X, real matrix Y) {

    for (i=1 ; i<=cols(X) ; i++) {
		"Compare column: "; i;
		compare_vec(X[.,i],Y[.,i])
	}

}

void compare_vec(real vector X, real vector Y) {

	"reldif"
	(X,Y,reldif(X,Y))[(1::10),.]
	"maxreldif"
	mreldif(X,Y)
	"mean covariance"
	meanvariance((X,Y))
	"correlation"
	correlation((X,Y))

}


/*******************************************************************************
 kvpeer: construct and declare contstraint vector for peer effects estimation.
 *** Cannot get equation indices without first initializing the problem.
   * Initialize problem with moptimize_query()
*******************************************************************************/

void kvpeer(transmorphic scalar M) {
    moptimize_query(M)
    i = moptimize_util_eq_indices(M,2,3)
    C = J(1,i[2,2],0)
    C[1]      = 1
    C[i[1,1]] = 1
    Cc = (C,1) 
    moptimize_init_constraints(M,(Cc))
}


end

