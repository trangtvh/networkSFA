
program networkSFA_p
	//version 2
	syntax newvarname(max=1) [if] [in],  [ XB Residuals U M JLMS BC EQuation(string)]
	marksample touse, novarlist
	local nopts : word count `residuals' `xb' `jlms' `bc' `u' `m'
	//display "number of opt is `nopts'"
	if `nopts' > 1 {
		display "{err}only one statistic may be specified"
		exit 498
	}
	tempname k_eq
	scalar `k_eq'=e(k_eq)
	//local var`i' `: word `i' of `equation''
	//display "`equation'"
	if "`equation'"!="" {
		local eqspec = subinstr("`equation'", "#", "", .)
		local eqspec=subinstr("`eqspec'", " ", "", .)
		//display "`eqspec'"
		if `eqspec' > `k_eq' {
		display "{err}the number specified in Equation() is greater than the number of equations"
		exit
		}
	}
	tempname b
	matrix `b'=e(b)
	local k_eq=e(k_eq)
	forval i=1/`k_eq' {
		tempname xb_`i'
		matrix score double `xb_`i''=`b' if `touse', eq(#`i')
		//summarize `xb_`i''
	}
	local depvar "`e(depvar)'"
	
	if (`nopts' == 0 | "`xb'" != "") {
		display "(option xb assumed; fitted values)"
		if "`equation'"=="" gen double `varlist' =	`xb_1' if `touse'
		else if `eqspec'==0 {
			display "{err} Eq(#0) is only used for options: jlms, bc, u, m"
			exit
		}
		else gen double `varlist' = `xb_`eqspec'' if `touse'
	}	
	else if ("`residuals'"!= "") {
		if "`equation'"=="" {
			local depvar `: word 1 of `depvar''
			tempvar y
			gen `y' = `depvar'
			gen double `varlist' = `y' -`xb_1' if `touse'
		} 
		else if `eqspec'==0 {
			display "{err} Eq(#0) is only used for options: jlms, bc, u, m"
			exit
		}
		else {
			local depvar `: word `eqspec' of `depvar''
			tempvar y
			gen `y' = `depvar'
			gen double `varlist' = `y' - `xb_`eqspec'' if `touse'
		}
	}		
	else {
		tempname Sigma_u Sigma_v b
		matrix `Sigma_u'=e(Sigma_u)
		matrix `Sigma_v'=e(Sigma_v)
		matrix `b'=e(b)
		local eqname "`e(eqnames)'"
		forval i=1/`k_eq' {
			tempname sig_u`i' sig_v`i'
			tempvar mu`i' e`i' y`i' u`i' m`i' bc`i' j`i' sig`i'
			scalar `sig_u`i''=`Sigma_u'[1,`i']
			scalar `sig_v`i''=`Sigma_v'[1,`i']
			local depvar`i' `: word `i' of `depvar''
			local eqname`i'  `: word `i' of `eqname''
			gen double `y`i'' = `depvar`i''
			gen double `e`i'' = `y`i'' -`xb_`i''
		}
		
		if e(case)=="FO" {
			forval i=1/`k_eq' {
				gen double `mu`i'' = `sig_u`i''^2 * `e`i''/(`sig_u`i''^2 +`sig_v`i''^2 )
				gen double `sig`i''= `sig_u`i'' * `sig_v`i''/ ((`sig_u`i''^2 +`sig_v`i''^2 )^0.5)
				gen double `u`i'' = `mu`i'' + `sig`i''* normalden(`mu`i''/`sig`i'')/normal(`mu`i''/`sig`i'')
				gen double `j`i'' = exp(-`u`i'')
				gen double `m`i'' =`mu`i''
				replace `m`i''= 0 if `mu`i''<0
				gen double `bc`i'' = exp((-1)*`mu`i''+1/2*`sig`i''^2)*normal(`mu`i''/`sig`i''-`sig`i'')/normal(`mu`i''/`sig`i'')
			}		
			
			if ("`u'"!="") {
				if "`equation'"==""  gen double `varlist' =`u1' if `touse'
				else if `eqspec'==0 {
					tempvar u
					gen double `u'=`u1'
					forval i=2/`k_eq' {
						local cname`i' = "`eqname1'"+":"+ "`depvar`i''"
						replace `u'=`u' + `u`i'' * `b'[1,"`cname`i''"]
					}
					gen double `varlist'=`u' if `touse'
				}
				else gen double `varlist'=`u`eqspec'' if `touse'	
			}
			
			if ("`m'"!="") {
				if "`equation'"==""  gen double `varlist' =`m1' if `touse'
				else if `eqspec'== 0 {
					tempvar m
					gen double `m'=`m1'
					forval i=2/`k_eq' {
						local cname`i' = "`eqname1'"+":"+ "`depvar`i''"
						replace `m'=`m' + `m`i'' * `b'[1,"`cname`i''"]
					}
					gen double `varlist'=`m' if `touse'
				}
				else gen double `varlist'=`m`eqspec'' if `touse'
			}
			
			if ("`jlms'"!="")  {
				if "`equation'"==""  gen double `varlist' =`j1' if `touse'
				else if `eqspec'== 0 {
					tempvar j
					gen double `j'=`j1'
					forval i=2/`k_eq' {
						local cname`i' = "`eqname1'"+":"+ "`depvar`i''"
						replace `j'=`j'* (`j`i''^`b'[1,"`cname`i''"])
					}
					gen double `varlist'=`j' if `touse'
				}
				else gen double `varlist'=`j`eqspec'' if `touse'
			}
			
			if ("`bc'"!="") { 
				if "`equation'"==""  gen double `varlist' =`bc1' if `touse'
				else if `eqspec'== 0 {
					tempvar bc
					gen double `bc'=`bc1'
					forval i=2/`k_eq' {
						local cname`i' = "`eqname1'"+":"+ "`depvar`i''"
						replace `bc'=`bc'* (`bc`i''^`b'[1,"`cname`i''"])
					}
					gen double `varlist'=`bc' if `touse'
				}
				else gen double `varlist'=`bc`eqspec'' if `touse'
			}
		}
		
		if e(case)=="I" {
			forval i=1/`k_eq' {
				gen double `mu`i'' = (-1)*`sig_u`i''^2 * `e`i''/(`sig_u`i''^2 +`sig_v`i''^2 )
				gen double `sig`i''= `sig_u`i'' * `sig_v`i''/ ((`sig_u`i''^2 +`sig_v`i''^2 )^0.5)
				gen double `u`i'' = `mu`i'' + `sig`i''* normalden(`mu`i''/`sig`i'')/normal(`mu`i''/`sig`i'')
				gen double `j`i'' = exp(-`u`i'')
				gen double `m`i'' = `mu`i''
				replace `m`i'' = 0 if `mu`i''<0
				gen double `bc`i'' = exp((-1)*`mu`i''+1/2*`sig`i''^2)*normal(`mu`i''/`sig`i''-`sig`i'')/normal(`mu`i''/`sig`i'')			
			}
			if ("`u'"!="") {
				if "`equation'"==""  gen double `varlist' =`u1' if `touse'
				else if `eqspec'== 0 {
					tempvar u
					gen double `u'=`u`k_eq''
					local n = `k_eq'-1
					forval i=1/`n' {
						local cname`i' = "`eqname`k_eq''"+":"+ "`depvar`i''"
						replace `u'=`u' + `u`i'' * `b'[1,"`cname`i''"]
						//display "`cname`i''"
					}
					gen double `varlist'=`u' if `touse'
				}
				else gen double `varlist'=`u`eqspec'' if `touse'
			}
			if ("`m'"!="") { 
				if "`equation'"==""  gen double `varlist' =`m1' if `touse'
				else if `eqspec'== 0 {
					tempvar m
					gen double `m'=`m`k_eq''
					local n = `k_eq'-1
					forval i=1/`n' {
						local cname`i' = "`eqname`k_eq''"+":"+ "`depvar`i''"
						replace `m'=`m' + `m`i'' * `b'[1,"`cname`i''"]
					}
					gen double `varlist'=`m' if `touse'
				}
				else gen double `varlist'=`m`eqspec'' if `touse'		
			}
			if ("`jlms'"!="") {
				if "`equation'"==""  gen double `varlist' =`j1' if `touse'
				else if `eqspec'== 0 {
					tempvar j
					gen double `j'=`j`k_eq''
					local n = `k_eq'-1
					forval i=1/`n' {
						local cname`i' = "`eqname`k_eq''"+":"+ "`depvar`i''"
						replace `j'=`j' * (`j`i'' ^ `b'[1,"`cname`i''"])
					}
					gen double `varlist'=`j' if `touse'
				}
				else gen double `varlist'=`j`eqspec''	 if `touse'
			}
			if ("`bc'"!="") {
				if "`equation'"==""  gen double `varlist' =`bc1' if `touse'
				else if `eqspec'== 0 {
					tempvar bc
					gen double `bc'=`bc`k_eq''
					local n = `k_eq'-1
					forval i=1/`n' {
						local cname`i' = "`eqname`k_eq''"+":"+ "`depvar`i''"
						replace `bc'=`bc' * (`bc`i'' ^ `b'[1,"`cname`i''"])
					}
					gen double `varlist'=`bc' if `touse'
				}
				else gen double `varlist'=`bc`eqspec'' if `touse'					
			}		
		}
		
		if e(case)=="IO" {
			gen double `mu1' = `sig_u1'^2 * `e1'/(`sig_u1'^2 +`sig_v1'^2 )		
			forval i=2/`k_eq' {
				gen double `mu`i'' = (-1)*`sig_u`i''^2 * `e`i''/(`sig_u`i''^2 +`sig_v`i''^2 )
			} 
			local IOvar "`e(IOvar)'"
			forval i=1/`k_eq' {
				local IO`i'  `: word `i' of `IOvar''
				gen double `sig`i''= `sig_u`i'' * `sig_v`i''/ ((`sig_u`i''^2 +`sig_v`i''^2 )^0.5)
				gen double `u`i'' = `mu`i'' + `sig`i''* normalden(`mu`i''/`sig`i'')/normal(`mu`i''/`sig`i'')
				gen double `j`i'' = exp(-`u`i'')
				gen double `m`i'' = `mu`i''
				replace `m`i'' = 0 if `mu`i''<0
				gen double `bc`i'' = exp((-1)*`mu`i''+1/2*`sig`i''^2)*normal(`mu`i''/`sig`i''-`sig`i'')/normal(`mu`i''/`sig`i'')			
			}
			
			if ("`u'"!="") {
				if "`equation'"==""  gen double `varlist' =`u1' if `touse'
				else if `eqspec'== 0 {
					tempvar u
					gen double `u'=`u1'
					forval i=2/`k_eq' {
						local j=`i'-1
						local cname1`i' = "`eqname1'"+":"+ "`IO`j''"
						local cname2`i' = "`eqname`i''"+":"+ "`IO`j''"
						replace `u'=`u'+`u`i'' * `b'[1,"`cname1`i''"]/`b'[1,"`cname2`i''"]	
					}
					gen double `varlist'=`u' if `touse'
				}
				else gen double `varlist'=`u`eqspec'' if `touse'
			} 
			
			if ("`m'"!="") {
				if "`equation'"==""  gen double `varlist' =`m1' if `touse'
				else if `eqspec'== 0 {
					tempvar m
					gen double `m'=`m1'
					forval i=2/`k_eq' {
						local j=`i'-1
						local cname1`i' = "`eqname1'"+":"+ "`IO`j''"
						local cname2`i' = "`eqname`i''"+":"+ "`IO`j''"
						replace `m'=`m'+`m`i'' * `b'[1,"`cname1`i''"]/`b'[1,"`cname2`i''"]	
					}
					gen double `varlist'=`m' if `touse'
				}
				else gen double `varlist'=`m`eqspec'' if `touse'
			}
			
			if ("`jlms'"!="")  {
				if "`equation'"==""  gen double `varlist' =`j1' if `touse'
				else if `eqspec'== 0 {
					tempvar jj
					gen double `jj'=`j1'
					forval i=2/`k_eq' {
						local j=`i'-1
						local cname1`i' = "`eqname1'"+":"+ "`IO`j''"
						local cname2`i' = "`eqname`i''"+":"+ "`IO`j''"
						replace `jj'=`jj'*`j`i'' ^ (`b'[1,"`cname1`i''"]/`b'[1,"`cname2`i''"])
					}
					gen double `varlist'=`jj' if `touse'
				}
				else gen double `varlist'=`j`eqspec'' if `touse'
			}
			
			if ("`bc'"!="")  {
				if "`equation'"==""  gen double `varlist' =`bc1' if `touse'
				else if `eqspec'== 0 {
					tempvar bc
					gen double `bc'=`bc1'
					forval i=2/`k_eq' {
						local j=`i'-1
						local cname1`i' = "`eqname1'"+":"+ "`IO`j''"
						local cname2`i' = "`eqname`i''"+":"+ "`IO`j''"
						replace `bc'=`bc'*`bc`i'' ^ (`b'[1,"`cname1`i''"]/`b'[1,"`cname2`i''"])
					}
					gen double `varlist'=`bc' if `touse'
				}
				else gen double `varlist'=`bc`eqspec'' if `touse'
			}
		}
	}	
end	
