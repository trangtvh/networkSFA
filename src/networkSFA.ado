program drop _all
program networkSFA, eclass
	//version 2
	syntax anything [if] [in]
	marksample touse
	local full_arg "`0'"
	//display "`full_arg'"
	local arg `anything'
	//display "`arg'"
	local n_open = length("`arg'") - length(subinstr("`arg'", "(", "", .))
	local n_close = length("`arg'") - length(subinstr("`arg'", ")", "", .))
	//display "The number of open parentheses is `n_open'"
	//display "The number of close parentheses is `n_close'"
	if `n_open'!= `n_close' display "{err}unmatched open/close parenthesis"
	if (`n_open'== `n_close' & `n_open'<2 ) display "{err}there must be at least 2 equations" 
	local end_cha = substr("`arg'",-1,1)
	local eq1_open = strpos("`arg'", "(")
	local eq1_close =  strpos("`arg'", ")")
	if (`eq1_open'==0 | `eq1_close'==0) display "{err}arguments of equations must be in parentheses ()"
	local eq1_length = `eq1_close' - `eq1_open'+1
	local eq1_arg = substr("`arg'",`eq1_open', `eq1_length')
	//display "the first open parenthesis `eq1_open'"
	//display "the first close parenthesis `eq1_close'"
	//display "the first equation arguments `eq1_arg'"
	local eqn_open = strrpos("`arg'", "(")
	local eqn_close =  strrpos("`arg'", ")")
	local eqn_length = `eqn_close' - `eqn_open'+1
	local eqn_arg = substr("`arg'",`eqn_open', `eqn_length')
	//display "the last open parenthesis `eqn_open'"
	//display "the last close parenthesis `eqn_close'"
	//display "the last equation arguments `eqn_arg'"
	tokenize "`arg'", p(" ")
	local case "`1'"
	//display "`case'"
	local start_sys = `eq1_open'
	local length_sys = `eqn_close' - `start_sys'+1
	local system = substr("`arg'",`start_sys', `length_sys')
	//display "case is `case'"
	//display "system is `system'"
	display as text "(Note: the dependent variables should be in a log form)"
	if (`eq1_open'>0 & `eq1_close'>0 & `n_open'==`n_close' & `n_open'>=2) {
	 	
		if "`case'" == "FO" {
			if "`end_cha'" != ")" {
				display "{error}the arguments must finish with a close parenthesis if the argument starts with FO"
				exit
			} 
			else {
				quietly reg3 `system' if `touse', small dfk
				//estimates store regFOmodel
				local n `e(k_eq)'
				local N `e(N)'
				local depvar "`e(depvar)'"
				forval i=2/`n' {
					local depvar`i' `: word `i' of `depvar''
					local depvar1`i' = " `depvar`i''"+")"
					local check_depvar`i' = strpos("`eq1_arg'", " `depvar`i'' ") + strpos("`eq1_arg'", "`depvar1`i''")
					if `check_depvar`i'' ==0 {
						display "{error}`depvar`i'' does not appear in the first equation"
						exit
					}
				} 

				tempname b V Sigma k k_eq df_r ll ic sigma_u sigma_v
				matrix `b'=e(b)
				matrix `V'=e(V)
				matrix `Sigma'=e(Sigma)
				local endog = e(endog)
				local exog = e(exog)
				local eqnames=e(eqnames)
				scalar `k'= e(k)
				scalar `k_eq'=e(k_eq)
				scalar `df_r'=e(df_r)
				scalar `ll'=e(ll)
				scalar `ic'=e(ic)
				matrix `sigma_u'=J(1,`n',.)
				matrix `sigma_v'=J(1,`n',.)
				local name_sigu =""
				local name_sigv =""
				forval i=1/`n' {
					tempvar e`i' e`i'_sq e`i'_cub 
					tempname sum_e`i'_sq sum_e`i'_cub sig_u`i' sig_v`i'_sq sig_v`i' ccons`i' df_eq`i'
					predict `e`i'' if `touse', residuals eq(#`i')
					//gen `e`i''= e`i'
					//drop e`i'
					gen `e`i'_sq'= `e`i''^2
					gen `e`i'_cub'= `e`i''^3
					quietly sum `e`i'_sq'
					scalar `sum_e`i'_sq'=r(sum)
					quietly sum `e`i'_cub'
					scalar `sum_e`i'_cub'=r(sum)
					if `sum_e`i'_cub'<=0 {
						display "{error}non-positive sigma_u`i'"
						exit
					}				
					scalar `sig_u`i''=((-1/(`N'-1))*`sum_e`i'_cub'*(_pi/2)^(1/2)*(_pi/(_pi-4)))^(1/3)
					scalar `sig_v`i'_sq'=1/(`N'-1)*`sum_e`i'_sq' - (1-2/_pi)*`sig_u`i''^2
					if `sig_v`i'_sq' <=0 {
						display "{error}non-positive sigma_v`i'"
						exit
					}
					scalar `sig_v`i''=sqrt(`sig_v`i'_sq')					
					scalar `ccons`i'' = [#`i']_b[_cons] - (2/_pi)^(1/2)*`sig_u`i''
					local eqname`i' `: word `i' of `eqnames''
					local cname`i' = "`eqname`i''"+":_cons"
					local cpost`i'=colnumb(`b',"`cname`i''")
					matrix `b'[1,`cpost`i''] = `ccons`i''
					scalar `df_eq`i'' = e(df_m`i')
					ereturn scalar df_m`i'=`df_eq`i''
					matrix `sigma_u'[1,`i']= `sig_u`i''
					matrix `sigma_v'[1,`i']= `sig_v`i''				
					local name_sigu ="`name_sigu'" +" sigma_u`i'"
					local name_sigv ="`name_sigv'" +" sigma_v`i'"
				}

				matrix colnames `sigma_u'= `name_sigu'
				matrix rownames `sigma_u'= Estimate
				matrix colnames `sigma_v'= `name_sigv'
				matrix rownames `sigma_v'= Estimate
				display as text "Network SFA model"
				display as text "Number of obs = `N'"
				ereturn post `b' `V'
				ereturn display
				ereturn local predict "networkSFA_p"
				ereturn local cmd "networkSFA"
				ereturn local method "3SLS - MM"
				ereturn local depvar "`endog'"
				ereturn local exog "`exog'"	
				ereturn local eqnames "`eqnames'"
				ereturn scalar k = `k'
				ereturn scalar k_eq = `k_eq'
				ereturn scalar df_r = `df_r'
				ereturn scalar ll = `ll'
				ereturn scalar ic = `ic'
				ereturn matrix Sigma = `Sigma'
				ereturn matrix Sigma_u = `sigma_u'
				ereturn matrix Sigma_v = `sigma_v'
				ereturn local case "FO"
				matrix list e(Sigma_u)
				matrix list e(Sigma_v)	

				}
			}
		
	
		else if "`case'" == "I" {
			if "`end_cha'" != ")" {
				display "{error}the arguments must finish with a close parenthesis if the argument starts with I"
				exit
			}
			else {
				quietly reg3 `system' if `touse', small dfk
				//estimates store regImodel
				local n=`e(k_eq)' - 1
				local N `e(N)'
				local depvar "`e(depvar)'"
				forval i=1/`n' {
					local depvar`i' `: word `i' of `depvar''
					local check_depvar`i' = strpos("`eqn_arg'", " `depvar`i'' ") + strpos("`eqn_arg'", " `depvar`i'')")
					if `check_depvar`i'' ==0 {
						display "{error}`depvar`i'' does not appear in the last equation"
						exit
					}
				} 	
				tempname b V Sigma k k_eq df_r ll ic sigma_u sigma_v
				matrix `b'=e(b)
				matrix `V'=e(V)
				matrix `Sigma'=e(Sigma)
				local endog = e(endog)
				local exog = e(exog)
				local eqnames=e(eqnames)
				scalar `k'= e(k)
				scalar `k_eq'=e(k_eq)
				scalar `df_r'=e(df_r)
				scalar `ll'=e(ll)
				scalar `ic'=e(ic)
				matrix `sigma_u'=J(1,`e(k_eq)',.)
				matrix `sigma_v'=J(1,`e(k_eq)',.)
				local name_sigu =""
				local name_sigv =""				
				
				forval i=1/`e(k_eq)' {
					tempvar e`i' e`i'_sq e`i'_cub 
					tempname sum_e`i'_sq sum_e`i'_cub sig_u`i' sig_v`i'_sq sig_v`i' ccons`i' df_eq`i'
					predict `e`i'' if `touse', residuals eq(#`i')
					//gen `e`i''= e`i'
					//drop e`i'
					gen `e`i'_sq'= `e`i''^2
					gen `e`i'_cub'= `e`i''^3
					quietly sum `e`i'_sq'
					scalar `sum_e`i'_sq'=r(sum)
					quietly sum `e`i'_cub'
					scalar `sum_e`i'_cub'=r(sum)
					//display `sum_e`i'_cub'
					if `sum_e`i'_cub'>=0 {
						display "{error} non-positive sigma_u`i'"
						exit
					} 
					scalar `sig_u`i''=((1/(`N'-1))*`sum_e`i'_cub'*(_pi/2)^(1/2)*(_pi/(_pi-4)))^(1/3)
					scalar `sig_v`i'_sq'=1/(`N'-1)*`sum_e`i'_sq' - (1-2/_pi)*`sig_u`i''^2
					//display `sig_v`i'_sq'
					if `sig_v`i'_sq' <=0 {
						display "{error} non-positive sigma_v`i'"
						exit
					} 
					scalar `sig_v`i''=sqrt(`sig_v`i'_sq')
					scalar `ccons`i'' = [#`i']_b[_cons] + (2/_pi)^(1/2)*`sig_u`i''
					//display "cons`i' is `ccons`i'' "
					local eqname`i' `: word `i' of `eqnames''
					local cname`i' = "`eqname`i''"+":_cons"
					local cpost`i'=colnumb(`b',"`cname`i''")
					matrix `b'[1,`cpost`i''] = `ccons`i''				
					scalar `df_eq`i'' = e(df_m`i')
					ereturn scalar df_m`i'=`df_eq`i''
					matrix `sigma_u'[1,`i']= `sig_u`i''
					//display `sig_u`i''
					matrix `sigma_v'[1,`i']= `sig_v`i''				
					local name_sigu ="`name_sigu'" +" sigma_u`i'"
					local name_sigv ="`name_sigv'" +" sigma_v`i'"
				}
				//display "`name_sigu'"
				matrix colnames `sigma_u'= `name_sigu'
				matrix rownames `sigma_u'= Estimate
				matrix colnames `sigma_v'= `name_sigv'
				matrix rownames `sigma_v'= Estimate
				display as text "Network SFA model"
				display as text "Number of obs = `N'"
				ereturn post `b' `V'
				ereturn display
				ereturn local predict "networkSFA_p"
				ereturn local cmd "networkSFA"
				ereturn local method "3SLS - MM"
				ereturn local depvar "`endog'"
				ereturn local exog "`exog'"	
				ereturn local eqnames "`eqnames'"
				ereturn scalar k = `k'
				ereturn scalar k_eq = `k_eq'
				ereturn scalar df_r = `df_r'
				ereturn scalar ll = `ll'
				ereturn scalar ic = `ic'
				ereturn matrix Sigma = `Sigma'
				ereturn matrix Sigma_u = `sigma_u'
				ereturn matrix Sigma_v = `sigma_v'
				ereturn local case "I"
				matrix list e(Sigma_u)
				matrix list e(Sigma_v)	
			}
		}
	
		else if "`case'" == "IO" {
			if "`end_cha'" == ")" {
				display "{error}the arguments must finish with a list of intermediate outputs if the argument starts with IO"
				exit
			} 
			else {
				quietly sureg `system' if `touse', small dfk
				//estimates store regIOmodel
				local n=`e(k_eq)'
				local N `e(N)'
				local start_IO = `eqn_close'+1
				local length_IO = length("`arg'") - `start_IO'+1
				local IOvar = substr("`arg'",`start_IO', `length_IO')
				local n_IO =wordcount("`IOvar'")
				local depvar "`e(depvar)'"
				if `n_IO' != `n'-1 {
					display "{error}the number of intermediate outputs must be equal to the number of the second stage equations"
					exit
				} 
				local depvar "`e(depvar)'"
				local n_dep=wordcount("`depvar'")
				tokenize "`system'", p(")")
				forval i=1/`n_IO' {
					local j=`i'*2+1
					local Eq2_`i'= "``j''" +")"
					//local `Eq2_`i''="`Eq2_`i''"+")"
					//display "`Eq2_`i''"
					local IO`i' `: word `i' of `IOvar''
					local checkeq1_IO`i' = strpos("`eq1_arg'", " `IO`i'' ") + strpos("`eq1_arg'", " `IO`i'')")
					if `checkeq1_IO`i'' ==0 {
						display "{error} `IO`i'' does not appear in the first equation"
						exit
					} 
					
					local checkeq2_IO`i' =  strpos("`Eq2_`i''", " `IO`i'' ") + strpos("`Eq2_`i''", " `IO`i'')")
					local k=`i'+1
					if `checkeq2_IO`i'' ==0 {
						display ("{error} `IO`i'' does not appear in equation `k'")
						exit
					} 
				}
				forval k=1/`n_dep' {
					local depvar`k' `: word `k' of `depvar''
					forval  i=1/`n_IO' {
						if "`IO`i''" == "`depvar`k''" {
							display "{error} `IO`i'' must not a dependent variable if the argument starts with IO"
							exit
						} 
					}
				}
				tempname b V Sigma k k_eq df_r ll ic sigma_u sigma_v
				matrix `b'=e(b)
				matrix `V'=e(V)
				matrix `Sigma'=e(Sigma)
				local depvar = e(depvar)
				local exog = e(exog)
				local eqnames=e(eqnames)
				scalar `k'= e(k)
				scalar `k_eq'=e(k_eq)
				scalar `df_r'=e(df_r)
				scalar `ll'=e(ll)
				scalar `ic'=e(ic)
				matrix `sigma_u'=J(1,`n',.)
				matrix `sigma_v'=J(1,`n',.)
				local name_sigu =""
				local name_sigv =""			
				forval i=1/`n' {
					tempvar e`i' e`i'_sq e`i'_cub 
					tempname sum_e`i'_sq sum_e`i'_cub sig_u`i' sig_v`i'_sq sig_v`i' ccons`i' df_eq`i' 
					predict `e`i'' if `touse', residuals eq(#`i')
					//gen `e`i''=e`i'
					//drop e`i'
					gen `e`i'_sq'= `e`i''^2
					gen `e`i'_cub'= `e`i''^3
					quietly sum `e`i'_sq'
					scalar `sum_e`i'_sq'=r(sum)
					quietly sum `e`i'_cub'
					scalar `sum_e`i'_cub'=r(sum)
					local eqname`i' `: word `i' of `eqnames''
					local cname`i' = "`eqname`i''"+":_cons"
					local cpost`i'=colnumb(`b',"`cname`i''")	
					scalar `df_eq`i'' = e(df_m`i')
					ereturn scalar df_m`i'= `df_eq`i''
					local name_sigu ="`name_sigu'" +" sigma_u`i'"
					local name_sigv ="`name_sigv'" +" sigma_v`i'"
				}
				if `sum_e1_cub'<=0 {
					display "{error}non-positive sigma_u1"
					exit
				} 
				scalar `sig_u1'=((-1/(`N'-1))*`sum_e1_cub'*(_pi/2)^(1/2)*(_pi/(_pi-4)))^(1/3)
				scalar `sig_v1_sq'= 1/(`N'-1)*`sum_e1_sq' - (1-2/_pi)*`sig_u1'^2
				if `sig_v1_sq'<=0 {
					display "{error}non-positive sigma_v1"
					exit
				} 
				scalar `sig_v1'=sqrt(`sig_v1_sq')
				scalar `ccons1' = [#1]_b[_cons] - (2/_pi)^(1/2)*`sig_u1'
				matrix `b'[1,`cpost1'] = `ccons1'
				matrix `sigma_u'[1,1]= `sig_u1'
				matrix `sigma_v'[1,1]= `sig_v1'
				forval i=2/`n' {
					if `sum_e`i'_cub'>=0 {
						display "{error} non-positive sigma_u`i'"
						exit
					} 
					scalar `sig_u`i''=((1/(`N'-1))*`sum_e`i'_cub'*(_pi/2)^(1/2)*(_pi/(_pi-4)))^(1/3)
					scalar `sig_v`i'_sq'=1/(`N'-1)*`sum_e`i'_sq' - (1-2/_pi)*`sig_u`i''^2
					if `sig_v`i'_sq' <=0{
						display "{error} non-positive sigma_v`i'"
						exit
					} 
					scalar `sig_v`i''=sqrt(`sig_v`i'_sq')
					scalar `ccons`i'' = [#`i']_b[_cons] + (2/_pi)^(1/2)*`sig_u`i''
					matrix `b'[1,`cpost`i''] = `ccons`i''
					matrix `sigma_u'[1,`i']= `sig_u`i''
					matrix `sigma_v'[1,`i']= `sig_v`i''	
				}
				matrix colnames `sigma_u'= `name_sigu'
				matrix rownames `sigma_u'= Estimate
				matrix colnames `sigma_v'= `name_sigv'
				matrix rownames `sigma_v'= Estimate
				display as text "Network SFA model"
				display as text "Number of obs = `N'"
				ereturn post `b' `V'
				ereturn display
				ereturn local predict "networkSFA_p"
				ereturn local cmd "networkSFA"
				ereturn local method "SURE - MM"
				ereturn local depvar "`depvar'"
				ereturn local exog "`exog'"	
				ereturn local eqnames "`eqnames'"
				ereturn scalar k = `k'
				ereturn scalar k_eq = `k_eq'
				ereturn scalar df_r = `df_r'
				ereturn scalar ll = `ll'
				ereturn scalar ic = `ic'
				ereturn matrix Sigma = `Sigma'
				ereturn matrix Sigma_u = `sigma_u'
				ereturn matrix Sigma_v = `sigma_v'
				ereturn local case "IO"
				matrix list e(Sigma_u)
				matrix list e(Sigma_v)
				ereturn local IOvar `IOvar'
			}
		}
	
		else {
			display "{error} the first arguments must be one of FO, IO, or I"
		}	
	}
	end

	
