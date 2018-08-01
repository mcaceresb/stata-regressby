*===============================================================================
* PROGRAM: regressby.ado
* PURPOSE: Performs fast grouped univariate OLS regressions.
*          The following commands are equivalent:
*			 regressby y x, by(byvars)
*			 statsby, by(byvars) clear: reg y x
*		   Except regressby will run 10-100x faster.
*          Also computes standard errors in a variety of flavors: usual
*          asymptotic standard errors, robust standard errors, and clustered
*          standard errors.
* AUTHORS: Michael Stepner, Michael Droste, Wilbur Townsend
*===============================================================================


*-------------------------------------------------------------------------------
* Stata wrapper
*-------------------------------------------------------------------------------

program define regressby

    * timer off 99
    * timer clear 99
    * timer on 99

	version 12.0
	syntax varlist(min=2 numeric) [aweight], by(varlist) [ ///
        vce(string)                                        ///
        save(string)                                       ///
        replace                                            ///
        NOCOVariance                                       ///
        NOConstant                                         ///
    ]

	* Preserve dataset in case we crash
	preserve

    qui ds `by'
    local by `r(varlist)'

    * Restrict sample with if/in conditions. This also flags obs
    * missing independent or dependent variables, missing by variable
    * observations, and invalid or missing weights.

	marksample touse, strok
	markout `touse' `by', strok

    gettoken y x: varlist
    qui ds `x'
    local x `r(varlist)'

    * XX Remove co-linear variables by group
    * _rmcoll `x', `noconstant' expand

	* Parse VCE option, if specified
	if ( `"`vce'"' != "" ) {
		my_vce_parse , vce(`vce')
		local vcetype    "robust"
		local clusterby  "`r(clustervar)'"
		if (`"`vcetype'"'   == "robust") local robust = "robust"
		if (`"`clusterby'"' != "")       local robust = ""
	}

	* Check to make sure save data file path is valid
	if ( ("`replace'" == "") & (`"`savegraph'"' != "") ) {
		if regexm(`"`savegraph'"',"\.[a-zA-Z0-9]+$") confirm new file `"`save'"'
		else confirm new file `"`save'.dta"'
        local savedata 1
	}
    else local savedata 0

	* Error checking: can't specify both robust and clusterby
	if ( ("`robust'" != "") & ("`clusterby'" != "") ) {
		di as error "Error: can't specify both clustered and robust standard errors at once! Choose one."
		exit
	}

	* Display type of standard error chosen
	if ( ("`robust'" == "") & ("`clusterby'" == "") ) {
		di "Running regressby with normal OLS standard errors."
	}
	if ( "`robust'" != "" ) {
		di "Running regressby with robust standard errors."
	}
	if ( "`clusterby'" != "" ) {
		di "Running regressby with cluster-robust standard errors (clustered by `clusterby')."
	}

	* Construct analytical weight variable
    tempvar tmpwt
	if ( "`weight'" != "" ) {
		local wt [`weight'`exp']
		gen `tmpwt' `exp'
		di "Using analytical weights, weight `exp'."

        * Display weighting scheme, if applicable
		foreach v in `varlist' {
			qui replace `v' = `v' * sqrt(`tmpwt')
		}
		qui replace `tmpwt' = sqrt(`tmpwt')
	}
    else if ( "`noconstant'" == "" ) {
        gen byte `tmpwt' = 1
    }

	* Sort using by-groups
    qui keep if `touse'
	sort `by' `clusterby'

    * Generate single by-variable counting by groups
    local type = cond(`=_N' < maxlong(), "long", "double")
    tempvar grp
    by `by': gen `type' `grp' = (_n == 1)
    qui replace `grp' = sum(`grp')
    * sort `grp'

    if ( "`clusterby'" != "" ) {
        tempvar cluster
        * by `grp' `clusterby': gen `type' `cluster' = (_n == 1)
        by `by' `clusterby': gen `type' `cluster' = (_n == 1)
        qui replace `cluster' = sum(`cluster')
    }

	* Also count number of variables here including constant
	local num_x: list sizeof x
	local num_g: list sizeof by
	if ( "`noconstant'" == "" ){
        local ++num_x
        local cons 1
        local x `x' `tmpwt'
    }
    else {
        local cons 0
    }
    scalar num_x = `num_x'
    scalar num_g = `num_g'

	* Drop groups if num obs < parameters (NOTE: can do this ex post)
	* qui by `grp': drop if (_N < `num_x')
	* qui by `by': drop if (_N < `num_x')

	* Whether to keep covariances
    local covs = ("`nocovariance'" == "")

    * timer off 99
    * timer list
    * timer clear 99
    * timer on 99

	* Perform regressions on each by-group, store in dataset
	mata: _regressby("`y'", "`x'", "`grp'", "`by'", "`cluster'", "`robust'", `cons', `covs')

    * timer off 99
    * timer list
    * timer clear 99

	* Optionally save out to dta and just restore with a message
	if ( "`save'" == "" ) {
		restore, not
	}
	if ( "`save'" != "" ) {
		save `save', `replace'
	}

    cap scalar drop num_x
    cap scalar drop num_g
end

*-------------------------------------------------------------------------------
* Mata program: _regressby3
* Inputs:
*  	- A y-var and x-var for an OLS regression
*  	- A group var, for which each value represents a distinct by-group.
*		This var must be in ascending order.
*	- A list of numeric by-variables, whose groups correspond to th group var.
* Outputs:
*  	- dataset of coefficients from OLS regression for each by-group
*-------------------------------------------------------------------------------

* yvar      = "`y'"
* xvars     = "`x'"
* grpvar    = "`grp'"
* byvars    = "`by'"
* clusterby = "`cluster'"
* robust    = "`robust'"
* cons      = `cons'
* covs      = `covs'

version 13.1
set matastrict on

mata:
real matrix _regressby(
    string scalar yvar,
    string scalar xvars,
    string scalar grpvar,
    string scalar byvars,
    string scalar clusterby,
    string scalar robust,
    real   scalar cons,
    real   scalar covs)
{

    // class Factor scalar F
    // real matrix XX, Xy

    real scalar N, nj, i, j, k, numgrp, group, starts, ends
    real scalar kx, hac0, hac1, cluster, nc, kvars
    real matrix info, cinfo, coefs, Vs
    real matrix X, _X, XX_inv, V, Z
    real matrix xi, ei

    real rowvector covix
    real colvector beta, e, y, _y, cvar, _cvar, _grp, nobs

    string rowvector covariates, addvars, addtypes

    // Views of y and X
    st_view(_y,    ., yvar)
    st_view(_X,    ., xvars)
    st_view(_cvar, ., clusterby)
    st_view(_grp,  ., grpvar)

    // Panel setup by numeric group variable
    // F = factor(grpvar)
    // F.panelsetup()
    // info   = F.info
    // info   = panelsetup(_grp, 1)
    N      = st_nobs()
    numgrp = _grp[N]
    info   = J(numgrp, 2, .)

    j = 1
    group = _grp[1]
    info[j, 1] = 1
    for (i = 2; i <= N; i++) {
        if ( group != _grp[i] ) {
            group = _grp[i]
            info[j++, 2] = i - 1
            info[j, 1]   = i
        }
    }
    info[j, 2] = N

    // i = 1
    // j = 1
    // group = _grp[i]
    // info[j, 1] = i
    // while ( j < numgrp ) {
    //     if ( group != _grp[++i] ) {
    //         group = _grp[i]
    //         info[j++, 2] = i - 1
    //         info[j, 1]   = i
    //     }
    // }
    // info[j, 2] = st_nobs()

    // Preallocate matrices for output
    kx      = st_numscalar("num_x")
    coefs   = J(numgrp * kx, 2,    .)
    Vs      = J(numgrp,      kx^2, .)
    nobs    = info[., 2] :- info[., 1] :+ 1
    // nobs    = J(numgrp, 1, .)
    hac0    = (robust == "") & (clusterby == "")
    hac1    = (robust != "") & (clusterby == "")
    cluster = (robust == "") & (clusterby != "")

    // stata("timer off 99")
    // stata("timer list")
    // stata("timer clear 99")
    // stata("timer on 99")

    // ------------------------------------------------------------------------
    // Iterate over groups
    // ------------------------------------------------------------------------

    // Iterate over groups 1 to Ng-1
    // ends = 1
    // group = _grp[ends]
    for (j = 1; j <= numgrp; j++) {
        // starts = ends
        // do {
        //     ends++
        // } while ( (group == _grp[ends]) & (ends < N) )

        // group = _grp[ends]
        // if ( ends == N ) {
        //     ends++
        // }

        nj = nobs[j]
        // nj      = ends - starts
        // nobs[j] = nj
        if ( nj < kx )
            continue

        st_subview(y, _y, info[j, .], .)
        st_subview(X, _X, info[j, .], .)
        // st_subview(y, _y, (starts, ends - 1), .)
        // st_subview(X, _X, (starts, ends - 1), .)

        // ------------ COMPUTE COEFFICIENTS --------------------
        XX_inv = invsym(quadcross(X, X))
        beta   = XX_inv * quadcross(X, y)
        e      = y - X * beta
        k      = kx - diag0cnt(XX_inv)

        // ------------ COMPUTE STANDARD ERRORS -----------------
        if ( hac0 ) {
            V 	= XX_inv * sum(e:^2) / (nj - k)
        }
        else if ( hac1 ) {
            V   = (nj / (nj - k)) * XX_inv * quadcross(X, e:^2, X) * XX_inv
        }
        else if ( cluster ) {
            st_subview(cvar, _cvar, info[j, .], .)
            // st_subview(cvar, _cvar, (starts, ends - 1), .)
            cinfo = panelsetup(cvar, 1)
            nc    = rows(cinfo)
            Z     = J(k, k, 0)
            // TODO: What if only one or two clusters?
            if ( nc > 2 ) {
                for (i = 1; i <= nc; i++) {
                    xi = panelsubmatrix(X, i, cinfo)
                    ei = panelsubmatrix(e, i, cinfo)
                    Z  = Z + xi' * (ei * ei') * xi
                }
                V = ((nj - 1) / (nj - k)) * (nc / (nc - 1)) * XX_inv * Z * XX_inv
            }
        }

        // TODO: Improve memory use for vcov matrix
        coefs[|1 + (j - 1) * kx, . \ j * kx, 2|] = beta, sqrt(diagonal(V))
        if ( covs ) {
            Vs[j, .] = rowshape(V, 1)
        }
    }

    // stata("timer off 99")
    // stata("timer list")
    // stata("timer clear 99")
    // stata("timer on 99")

    // ------------------------------------------------------------------------
    // Gather output and pass back into Stata
    // ------------------------------------------------------------------------

    // Add all the new variables at once: nj, b, se, cov
    covariates = tokens(xvars)
    kvars      = kx + kx + covs * (kx * (kx - 1) / 2)
    addvars    = J(1, 1 + kvars, "")
    covix      = J(1, covs * (kx * (kx - 1) / 2), .)
    addvars[1] = "N"
    addtypes   = st_local("type"), J(1, kvars, st_global("c(type)"))
    if ( cons ) {
        covariates[kx] = "cons"
    }

    for (k = 1; k <= kx; k++) {
        addvars[2 * k]     = "_b_"  + covariates[k]
        addvars[2 * k + 1] = "_se_" + covariates[k]
    }

    if ( covs ) {
        i = 2
        for (k = 1; k <= kx; k++) {
            for (j = 1; j < k; j++) {
                covix[i - 1] = (k - 1) * kx + j
                // printf("%g: %s\n", 2 * kx + i++, "_cov_" + covariates[k] + "_" + covariates[j])
                addvars[2 * kx + i++] = "_cov_" + covariates[k] + "_" + covariates[j]
            }
        }
    }

    // Keep only one obs per group
    // stata(sprintf("qui by %s: keep if (_n == 1)", grpvar))
    stata(sprintf("qui by %s: keep if (_n == 1)", byvars))
    stata("keep " + byvars)
    (void) st_addvar(addtypes, addvars)

    // Number of observations, betas, se, covs
    if ( covs ) {
        st_store(., addvars, (nobs, rowshape(coefs, numgrp), Vs[., covix]))
    }
    else {
        st_store(., addvars, (nobs, rowshape(coefs, numgrp)))
    }

    stata(sprintf("qui drop if N < %g", kx))
}
end

*-------------------------------------------------------------------------------
* Auxiliary Stata programs for parsing vce (standard error options)
* Source: https://blog.stata.com/2015/12/08/programming-an-estimation-command-in-stata-using-a-subroutine-to-parse-a-complex-option/
*-------------------------------------------------------------------------------

program define my_vce_parse, rclass
    syntax  [, vce(string) ]
    local case : word count `vce'
    if `case' > 2 {
        my_vce_error , typed(`vce')
    }
    local 0 `", `vce'"'
    syntax  [, Robust CLuster * ]
    if `case' == 2 {
        if "`robust'" == "robust" | "`cluster'" == "" {
            my_vce_error , typed(`vce')
        }
        capture confirm numeric variable `options'
        if _rc {
            my_vce_error , typed(`vce')
        }
        local clustervar "`options'"
    }
    else {    // case = 1
        if "`robust'" == "" {
            my_vce_error , typed(`vce')
        }
    }
    return clear
    return local clustervar "`clustervar'"
end

program define my_vce_error
    syntax , typed(string)
    display `"{red}{bf:vce(`typed')} invalid"'
    error 498
end
