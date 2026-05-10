*! v.1.0.0 N.Orsini, R. Bellocco, A. Sjolander 25sep12
capture program drop drglm_momev_ebe 
program drglm_momev_ebe 
	version 12	
	syntax varlist if, at(name) [ ///
	mainlist(string) ///
	exposurelist(string) ///
	y(string) ///
	a(string) ///
	derivatives(varlist) ///
	olink(string) ///
	elink(string) ]  
	
	forvalues i = 1/2 {
        	local mu`i' : word `i' of `varlist'
	}
	local i 1
	
	tempvar eqA 
	quietly gen double `eqA' = 0 `if'
	
	if "`exposurelist'"!="" {
	foreach var of varlist `exposurelist' {
		quietly replace `eqA' = `eqA' + `var'*`at'[1,`i'] `if'
		local `++i'
	}
	}
	quietly replace `eqA' = `eqA' + `at'[1,`i'] `if' 
	local i = `i' + 1
	
	tempvar eqM 
	quietly gen double `eqM' = 0 `if'

	foreach var of varlist `mainlist'  {
		quietly replace `eqM' = `eqM' + `var'*`at'[1,`i'] `if'
		local `++i'
	}
	
	tempvar Ua
	tempvar Uym
	tempvar Uaym
	if ("`olink'" == "identity"){
		quietly generate double `Uym' = `y' - `eqM' `if' 
	}
	if ("`olink'" == "log"){
		quietly generate double `Uym' = `y'/exp(`eqM') `if' 
	}
	if ("`elink'" == "identity"){
		quietly generate double `Ua' = `a' - `eqA' `if'
	}
	if ("`elink'" == "log"){
		quietly generate double `Ua' = `a' - exp(`eqA') `if'
	}
	if ("`elink'" == "logit"){
		quietly generate double `Ua' = `a' - invlogit(`eqA') `if'
	}
	quietly generate double `Uaym' = `Ua'*`Uym' `if'
	
	quietly replace `mu1' = `Ua' `if' 
	quietly replace `mu2' = `Uaym' `if'
	
	if "`derivatives'" == "" {
		exit
	}

	
	tempvar dUa
	tempvar dUymdm
	if ("`olink'" == "identity"){
		quietly generate double `dUymdm' = - 1 `if'
	}
	if ("`olink'" == "log"){
		quietly generate double `dUymdm' = - `y'/exp(`eqM') `if'
	}
	if ("`elink'" == "identity"){
		quietly generate double `dUa' = - 1 `if' 
	}
	if ("`elink'" == "log"){
		quietly generate double `dUa' = - exp(`eqA') `if' 
	}
	if ("`elink'" == "logit"){
		quietly generate double `dUa' = - invlogit(`eqA')*(1-invlogit(`eqA')) `if' 
	}
	
	
	local nd : word count `derivatives'
	local j = 1 
	
	// Differentiate exposure model 
	
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
 
 	foreach x of local exposurelist {
		quietly replace `: word `j' of `derivatives'' = `x'*`dUa'*`Uym'
		local `++j'
	}
	quietly replace `: word `j' of `derivatives'' = `dUa'*`Uym' 
	local j  = `j' + 1
		
	foreach x of local mainlist {
		quietly replace `: word `j' of `derivatives'' = `Ua'*`x'*`dUymdm'
		local `++j'
	}

end
