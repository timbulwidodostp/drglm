*! v.1.0.0 N.Orsini, R. Bellocco, A. Sjolander 25sep12
capture program drop drglm_momev_logit_logit_dr 
program drglm_momev_logit_logit_dr 
	version 12	
	syntax varlist if, at(name) [ ///
	mainmodel(string) ///
	mainlist(string) ///
	mainlist2(string) ///
	outcomelist(string)  ///
	exposurelist(string) ///
	y(string) ///
	a(string) ///
	derivatives(varlist)  ]  
	
	forvalues i = 1/3 {
        	local mu`i' : word `i' of `varlist'
	}
	local i 1
	
	tempvar eqY2
	quietly gen double `eqY2' = 0 `if'
	foreach var of varlist `mainlist' {
		quietly replace `eqY2' = `eqY2' + `var'*`at'[1,`i'] `if'
		local `++i'
	}
	tempvar eqY
	quietly gen double `eqY' = 0 `if'
	
	if "`outcomelist'" != "" {
	foreach var of varlist `outcomelist' {
		quietly replace `eqY' = `eqY' + `var'*`at'[1,`i'] `if'
		local `++i'
	}
	}
	quietly replace `eqY' = `eqY' + `at'[1,`i'] `if' 
	local i = `i' + 1
	tempvar eqA2
	 quietly gen double `eqA2' = 0 `if'
	foreach var of varlist `mainlist2' {
		quietly replace `eqA2' = `eqA2' + `var'*`at'[1,`i'] `if'
		local `++i'
	}
	tempvar eqA 
	quietly gen double `eqA' = 0 `if'
	if "`exposurelist'" != "" {
	foreach var of varlist `exposurelist' {
		quietly replace `eqA' = `eqA' + `var'*`at'[1,`i'] `if'
		local `++i'
	}
	}
	quietly replace `eqA' = `eqA' + `at'[1,`i'] `if' 
	local i = `i' + 1
	tempvar eqM 
	quietly gen double `eqM' = `at'[1,`i'] `if'
	local i = `i' + 1
	if ("`mainmodel'" != ""){
		foreach var of varlist `mainmodel'  {
			quietly replace `eqM' = `eqM' + `var'*`at'[1,`i'] `if'
			local `++i'
		}
	}
	
	tempvar meanY 
	quietly gen double `meanY' = invlogit(`a'*`eqM'+`eqY') `if'
	tempvar meanAstar
	quietly gen double `meanAstar' = 1/(1+((1-invlogit(`eqA'))*invlogit(`eqY'))/(invlogit(`eqA')*invlogit(`eqM'+`eqY'))) `if' 
	tempvar Uy
	quietly gen double `Uy' = `y' - invlogit(`eqY'+`eqY2') `if'
	tempvar Ua
	quietly gen double `Ua' = `a' - invlogit(`eqA'+`eqA2') `if'
	tempvar Uym
	quietly gen double `Uym' = `y'-`meanY' `if'
	tempvar Uam
	quietly gen double `Uam' = `a'-`meanAstar' `if'
	
	quietly replace `mu1' = `Uy' `if'
	quietly replace `mu2' = `Ua' `if' 
	quietly replace `mu3' = `Uam'*`Uym' `if'
	
	if "`derivatives'" == "" {
		exit
	}

	tempvar meanA meanYstar dmYdM dmAsdA scale dmAsdM dmYsdY
	quietly gen double `meanA' = invlogit(`y'*`eqM'+`eqA') `if'
	quietly gen double `meanYstar' = 1/(1+((1-invlogit(`eqY'))*invlogit(`eqA'))/(invlogit(`eqY')*invlogit(`eqM'+`eqA'))) `if' 
	quietly gen double `dmYdM' = `a'*`meanY'*(1-`meanY') `if'
	quietly gen double `dmAsdA'  = `meanAstar'*(1-`meanAstar') `if'
	quietly gen double `scale' = (1+exp(`eqY')+exp(`eqY'-`eqA'))/(1+exp(`eqY')) `if'
	quietly gen double `dmAsdM'  = `meanAstar'*(1-`meanAstar'*`scale') `if' 
	quietly gen double `dmYsdY' = `meanYstar'*(1-`meanYstar') `if'
	
	tempvar dUy
	quietly gen double `dUy' = - invlogit(`eqY'+`eqY2')*(1-invlogit(`eqY'+`eqY2')) `if'
	tempvar dUa
	quietly gen double `dUa' = - invlogit(`eqA'+`eqA2')*(1-invlogit(`eqA'+`eqA2')) `if'
	
	local nd : word count `derivatives'
	local j = 1 
	
	// Differentiate outcome model 
	foreach x of local mainlist {
		quietly replace `: word `j' of `derivatives'' = `x'*`dUy'
		local `++j'
	}
	foreach x of local outcomelist {
		quietly replace `: word `j' of `derivatives'' = `x'*`dUy'   
		local `++j'
	}
	quietly replace `: word `j' of `derivatives'' = `dUy' 
	local j  = `j' + 1
	foreach x of local mainlist2 {
		quietly replace `: word `j' of `derivatives'' = 0
		local `++j'
	}
 	foreach x of local exposurelist {
		quietly replace `: word `j' of `derivatives'' = 0
		local `++j'
	}
	quietly replace `: word `j' of `derivatives'' = 0 
	local j  = `j' + 1
	foreach x of local mainlist {
		quietly replace `: word `j' of `derivatives'' = 0
		local `++j'
	}
 
	// Differentiate exposure model 
	foreach x of local mainlist {
		quietly replace `: word `j' of `derivatives'' = 0
		local `++j'
	}
	foreach x of local outcomelist {
		quietly replace `: word `j' of `derivatives'' = 0
		local `++j'
	}
	quietly replace `: word `j' of `derivatives'' = 0 
	local j  = `j' + 1
	foreach x of local mainlist2 {
		quietly replace `: word `j' of `derivatives'' = `x'*`dUa'
		local `++j'
	}
	foreach x of local exposurelist {
		quietly replace `: word `j' of `derivatives'' = `x'*`dUa' 
		local `++j'
	}
	quietly replace `: word `j' of `derivatives'' = `dUa' 
	local j  = `j' + 1
	foreach x of local mainlist {
		quietly replace `: word `j' of `derivatives'' = 0
		local `++j'
	}

	// Differentiate main model
	foreach x of local mainlist {
		quietly replace `: word `j' of `derivatives'' = 0
		local `++j'
	}
	foreach x of local outcomelist {
		quietly replace `: word `j' of `derivatives'' = -1*`x'*`dmYsdY'*(`a'-`meanA')
		local `++j'
	}
	quietly replace `: word `j' of `derivatives'' =  -1*`dmYsdY'*(`a'-`meanA')
	local j  = `j' + 1
	foreach x of local mainlist2 {
		quietly replace `: word `j' of `derivatives'' = 0
		local `++j'
	}
 	foreach x of local exposurelist {
		quietly replace `: word `j' of `derivatives'' = -1*`x'*`dmAsdA'*`Uym'
		local `++j'
	}
	quietly replace `: word `j' of `derivatives'' = -1*`dmAsdA'*`Uym' 
	local j  = `j' + 1
	quietly replace `: word `j' of `derivatives'' =  -1*`dmAsdM'*`Uym'-1*`Uam'*`dmYdM'
	local j  = `j' + 1
	if ("`mainmodel'" != ""){
		foreach x of local mainmodel {
			quietly replace `: word `j' of `derivatives'' = -1*`x'*`dmAsdM'*`Uym'-`x'*`Uam'*`dmYdM'
			local `++j'
		}
	}

end
