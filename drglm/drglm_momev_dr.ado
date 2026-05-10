*! v.1.0.0 N.Orsini, R. Bellocco, A. Sjolander 25sep12
capture program drop drglm_momev_dr 
program drglm_momev_dr 
	version 12	
	syntax varlist if, at(name) [ ///
	mainlist(string) ///
	outcomelist(string)  ///
	exposurelist(string) ///
	y(string) ///
	a(string) ///
	derivatives(varlist) ///
	olink(string) ///
	elink(string) ]  
	
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
	if "`outcomelist'"!="" {
	 	foreach var of varlist `outcomelist' {
		quietly replace `eqY' = `eqY' + `var'*`at'[1,`i'] `if'
		local `++i'
		}
	}
	quietly replace `eqY' = `eqY' + `at'[1,`i'] `if' 
	local i = `i' + 1
	
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
		
	tempvar Uy
	tempvar Ua
	tempvar Uym
	tempvar Uaym
	if ("`olink'" == "identity"){
		quietly generate double `Uy' = `y' - (`eqY'+`eqY2') `if'
		quietly generate double `Uym' = `y' - (`eqM'+`eqY') `if' 
	}
	if ("`olink'" == "log"){
		quietly generate double `Uy' = `y' - exp(`eqY'+`eqY2') `if'
		quietly generate double `Uym' = `y'/exp(`eqM')-exp(`eqY') `if' 
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
	
	quietly replace `mu1' =  `Uy' `if'
	quietly replace `mu2' =  `Ua' `if' 
	quietly replace `mu3' = `Uaym' `if'
	
	if "`derivatives'" == "" {
		exit
	}

	tempvar dUy
	tempvar dUa
	tempvar dUymdy
	tempvar dUymdm
	if ("`olink'" == "identity"){
		quietly generate double `dUy' = - 1 `if'
		quietly generate double `dUymdy' = - 1 `if'
		quietly generate double `dUymdm' = - 1 `if'
	}
	if ("`olink'" == "log"){
		quietly generate double `dUy' = - exp(`eqY'+`eqY2') `if'
		quietly generate double `dUymdy' = - exp(`eqY') `if'
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
		quietly replace `: word `j' of `derivatives'' = `x'*`Ua'*`dUymdy'
		local `++j'
	}
	quietly replace `: word `j' of `derivatives'' =  `Ua'*`dUymdy'
	local j  = `j' + 1
	
 
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
