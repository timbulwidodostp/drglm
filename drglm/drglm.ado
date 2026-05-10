*! v.1.0.0 N.Orsini, R. Bellocco, A. Sjolander 25sep12
capture program drop drglm
program drglm, eclass properties(mi) byable(onecall)
version 12

if _by() {
		local BY `"by `_byvars'`_byrc0':"'
}

local version : di "version " string(_caller()) ":"
`version' `BY' _vce_parserun drglm, mark(CLuster): `0'

if "`s(exit)'" != "" {
		version 10: ereturn local cmdline `"drglm `0'"'
		exit
	}

if replay() {
		if ("`e(cmd)'"!="drglm")  error 301  
		Replay `0'
	}
else `version' `BY' Estimate `0'
ereturn local cmdline `"drglm `0'"'
end

capture program drop Estimate 
version 12
program Estimate, eclass byable(recall)
syntax  varlist(min=2 max=2 numeric) [if] [in] [, ///
		Main(string) ///
		Outcome(string) ///
		Exposure(string) ///
		OLink(string) ///
		ELink(string) /// 
		EForm ///
		ebe    /// 
		obe  ///
		Level(cilevel) 			///
		VCE(passthru) 			///
        ]
 
	local cmdline : copy local 0
	marksample touse 

	local vceopt =	`:length local vce'		|	///
	   		`:length local cluster'		|	///
	   		`:length local robust'

	if `vceopt' {
		_vce_parse, argopt(CLuster) 	///
			: [`exp'], `vce'  `cluster'
		local vce
		if "`r(cluster)'" != "" {
			local clustvar `r(cluster)'
			local vce vce(cluster `r(cluster)')
		}
	}
 
 
// Get outcome and exposure names 

	gettoken y a : varlist		

// Get the estimation method

						local method "dr"  // Double   Robust Estimator (Default)
	if ("`ebe'"!="")   	local method "ebe" // Exposure Based  Estimator
	if ("`obe'"!="") 	local method "obe" // Outcome  Based  Estimator

	if  ("`ebe'"!="") & ("`obe'"!="") {
					di as err "specify just one estimation method"
					exit 198
	}

// Understand the main, outcome, and exposure model

	if ("`olink'" != "") {
	if inlist("`olink'", "identity", "logit", "log") != 1 {
		di as err "option olink can be either identity, logit, or log"
		exit 198
		}
	}

	if ("`elink'" != "") {	
	if inlist("`elink'", "identity", "logit", "log") != 1 {
		di as err "option elink() can be either identity, logit, or log"
		exit 198
	}
	}
	
	// Default links
	
	if ("`olink'" == "") | ("`olink'" == "identity") {
		local olink = "identity"
	}
	if ("`elink'" == "") | ("`elink'" == "identity") {
		local elink = "identity"
	}
	
 	if ("`olink'" == "logit") {
		local elink = "logit"
	}
	
// Check the distribution of the outcome (0/1) for the logit-logit scenario

	if ("`olink'" == "logit") & ("`elink'" == "logit") {
		
		quietly  tabulate `y' 
		if r(r) != 2 {
				di as err "with olink(`olink') and elink(`elink') `y' can take on only two values (either 0 or 1)"
				exit 198
		}
		qui levelsof `y'
		if "`r(levels)'" != "0 1" {
				di as err "with olink(`olink') and elink(`elink') `y' can take on only two values (either 0 or 1)"
				exit 198
		}
		
	}
	
// Check the distribution of the exposure (0/1) for the logit-logit scenario

	if ("`olink'" == "logit") & ("`elink'" == "logit") {
		
		quietly  tabulate `a' 
		if r(r) != 2 {
				di as err "with olink(`olink') and elink(`elink') `a' can take on only two values (either 0 or 1)"
				exit 198
		}
		qui levelsof `a'
		if "`r(levels)'" != "0 1" {
				di as err "with olink(`olink') and elink(`elink') `a' can take on only two values (either 0 or 1)"
				exit 198
		}
		
	}
	
// Understand the main model 

	_rmcoll `main'   
	local mainmodel `r(varlist)' 

	local mainlist "`a'"
	
	qui foreach v of local mainmodel {
			  capture confirm new  variable `a'`v'
				if _rc != 0 {
					di in gr "variable `a'`v' is replaced"
					quietly drop `a'`v'
				}
			  qui gen `a'`v' =  `a'*`v'
			  label var `a'`v' "`a'*`v'"
			 local mainlist "`mainlist' `a'`v'"
	}
	
// Understand the outcome model
	
	if "`outcome'" != ""  {					 
		_rmcoll `outcome'    
		local outcomelist_orig `r(varlist)' 
	}


	// Standardize each variable of the outcome list
 	quietly foreach v of local outcomelist_orig  {		
		tempvar `v's
		qui su `v'
		gen ``v's' = (`v'-r(mean))/r(sd)
		local outcomelist "`outcomelist' ``v's'"
	}
 	
// Understand the exposure model

	if "`exposure'" != ""  {			
		_rmcoll `exposure'    
		local exposurelist_orig `r(varlist)' 
	}

	
// Standardize each variable of the exposure list
	quietly foreach v of local exposurelist_orig  {		
		tempvar `v's
		qui su `v'
		gen ``v's' = (`v'-r(mean))/r(sd)
		local exposurelist "`exposurelist' ``v's'"
	}		
		
	tempname getcoefs coefs VCE getVCE 

// Select the subset of observations 

		quietly {
				preserve
				keep if `touse'
				keep `y' `a' `clustvar' `mainmodel' `mainlist' `outcomelist_orig' `exposurelist_orig' `outcomelist' `exposurelist' `mainlist2'
		}

// Get variable names for the beta vector 

	    local conams "`mainlist'"
		
// Get scalars of interest

		local p = `: word count `mainlist''
	    local nobs `c(N)'

// Create the equations to be passed to GMM

	local eqA "{xA_0}"

	foreach v of local exposurelist {
				local eqA "`eqA'+`v'*{xA_`v'}"
	}	

	local eqM "{xmain_0}"
	foreach v of local mainmodel {
				local eqM "`eqM' +`v'*{xmain_`v'}"
	}	

	local eqY "{xY_0}"

	foreach v of local outcomelist {
				local eqY "`eqY'+`v'*{xY_`v'}"
	}
	
	if "`exposure'"=="" {		
			local Astr "{xA_0}"
	}
	else {
		   local Astr "{xA:`exposurelist'}+{xA_0}"
	}
	


if ("`olink'" != "logit") {
		
		tempname bA  by ep
		
		if "`clustvar'" != ""  qui glm `y' `mainlist' `outcomelist' , link("`olink'") `vce'
		else qui glm `y' `mainlist' `outcomelist' , link("`olink'") vce(robust)
 		mat `by' = e(b)
		
		if ("`method'"=="ebe") | ("`method'"=="dr"){
			
			*Construct parameter names and intial values for the paremeters
			
			if "`clustvar'" != ""  qui glm `a' `exposurelist' , link("`elink'") `vce'
			else qui glm `a' `exposurelist' , link("`elink'") vce(robust)
			mat `bA' = e(b)
			local initA ""
			local parA ""
			local initY ""
			local initM ""
			local parY ""
			local parM ""
			mat `ep' = `bA'[1,"`a':_cons"]
			local getb = `ep'[1,1]
			local initA "xA_0 `getb'"
			foreach v of local exposurelist {
				mat `ep' = `bA'[1,"`a':`v'"]
				local getb = `ep'[1,1]
				local initA "`initA' xA_`v' `getb'"
				local parA "`parA' xA_`v'"
			}
			local parA "`parA' xA_0"
			mat `ep' = `by'[1,"`y':_cons"]
			local getb = `ep'[1,1]
			local initY "xY_0 `getb'"
			foreach v of local mainlist {
				mat `ep' = `by'[1,"`y':`v'"]
				local getb = `ep'[1,1]
				local initY "`initY' xY_`v' `getb'"
				local parY "`parY' xY_`v'"
			}
			foreach v of local outcomelist {
				mat `ep' = `by'[1,"`y':`v'"]
				local getb = `ep'[1,1]
				local initY "`initY' xY_`v' `getb'"
				local parY "`parY' xY_`v'"
			}
			local parY "`parY' xY_0"
			local fw : word 1 of `mainlist'
			mat `ep' = `by'[1,"`y':`fw'"]
			local getb = `ep'[1,1]
			local initM "xmain_0 `getb'"
			tokenize `mainmodel'
			local k = 1
			foreach v of local mainlist {
				if `k' > 1 {
					local k = `k'-1
					mat `ep' = `by'[1,"`y':`v'"]
					local getb = `ep'[1,1]
					local initM "`initM' xmain_``k'' `getb'"
					local parM "`parM' xmain_``k''"	
				}
				local k = `k'+1
			}
			local parM "xmain_0 `parM'"
		}
			
		// Double Robust Estimator 

		if  ("`method'"=="dr") {
				 
		 	qui gmm  drglm_momev_dr  , 						///
					nequations(3) 						///
					parameters(`parY' `parA' `parM') 	///
					mainlist(`mainlist') ///
					outcomelist(`outcomelist')  ///
					exposurelist(`exposurelist') ///
					y(`y') ///
					a(`a') ///
					equations(Y A main)  ///
					winitial(unadjusted, independent) ///
					instruments(Y:`mainlist' `outcomelist') ///
				    instruments(A:`exposurelist') ///
				    instruments(main: `mainmodel')  from(`initY' `initA' `initM')  ///
					onestep hasderivatives ///
					olink(`olink') elink(`elink')  `vce'	
		}
 		// Exposure model estimator 
			
		if  ("`method'"=="ebe")  {
		 
			qui gmm  drglm_momev_ebe  , 						///
					nequations(2) 						///
					parameters(`parA' `parM') 	///
					mainlist(`mainlist') ///
					exposurelist(`exposurelist') ///
					y(`y') ///
					a(`a') ///
					equations(A main)  ///
					winitial(unadjusted, independent) ///
				    instruments(A:`exposurelist') ///
				    instruments(main: `mainmodel')  from(`initA' `initM')  ///
					onestep hasderivatives ///
					olink(`olink') elink(`elink') `vce'	
		
		}
		
		

}

if ("`olink'" == "logit") {
	
		tempname bA  by ep 

		if ("`method'"=="ebe") | ("`method'"=="dr"){
		
			*Construct parameter names and intial values for the parameters
		
			local mainlist2 "`y'"
			qui foreach v of local mainmodel {
				capture confirm new  variable `y'`v'
				if _rc != 0 {
					di as res "variable `y'`v' already in the dataset. It is replaced."
					quietly drop `y'`v'
				}
				gen `y'`v' =  `y'*`v'
				local mainlist2 "`mainlist2' `y'`v'"
			}	
			if "`clustvar'" != "" qui logit `a' `exposurelist' `mainlist2' , `vce'
			else qui logit `a' `exposurelist' `mainlist2' , vce(robust)
			mat `bA' = e(b)	
			local initA ""
			local parA ""
			mat `ep' = `bA'[1,"`a':_cons"]
			local getb = `ep'[1,1]
			local initA "xA_0 `getb'"
			foreach v of local mainlist2 {
				mat `ep' = `bA'[1,"`a':`v'"]
				local getb = `ep'[1,1]
				local initA "`initA' xA_`v' `getb'"
				local parA "`parA' xA_`v'"
			}
			foreach v of local exposurelist {
				mat `ep' = `bA'[1,"`a':`v'"]
				local getb = `ep'[1,1]
				local initA "`initA' xA_`v' `getb'"
				local parA "`parA' xA_`v'"
			}
			local parA "`parA' xA_0"
		}
	
	 if ("`method'"=="obe") | ("`method'"=="dr"){
	 
		*Construct parameter names and intial values for the parameters
			
		if "`clustvar'" != "" quietly logit `y' `mainlist' `outcomelist' , `vce'
		else quietly logit `y' `mainlist' `outcomelist' , vce(robust)
		mat `by' = e(b)
		local initY ""
		local initM ""
		local parY ""
		local parM ""
	    mat `ep' = `by'[1,"`y':_cons"]
		local getb = `ep'[1,1]
		local initY "xY_0 `getb'"
		foreach v of local mainlist {
			mat `ep' = `by'[1,"`y':`v'"]
			local getb = `ep'[1,1]
			local initY "`initY' xY_`v' `getb'"
			local parY "`parY' xY_`v'"
		}
		foreach v of local outcomelist {
			mat `ep' = `by'[1,"`y':`v'"]
			local getb = `ep'[1,1]
			local initY "`initY' xY_`v' `getb'"
			local parY "`parY' xY_`v'"
		}
		local parY "`parY' xY_0"
		local fw : word 1 of `mainlist'
		mat `ep' = `by'[1,"`y':`fw'"]
		local getb = `ep'[1,1]
		local initM "xmain_0 `getb'"
		tokenize `mainmodel'
		local k = 1
		foreach v of local mainlist {
			if `k' > 1 {
				local k = `k'-1
				mat `ep' = `by'[1,"`y':`v'"]
				local getb = `ep'[1,1]
				local initM "`initM' xmain_``k'' `getb'"
				local parM "`parM' xmain_``k''"	
			}
			local k = `k'+1
		}	
		local parM "xmain_0 `parM'"
	}
	
		// Double Robust Estimator

		if  ("`method'"=="dr") {
				 
		 	 qui gmm  drglm_momev_logit_logit_dr  , 						///
					nequations(3) 						///
					parameters(`parY' `parA' `parM') ///
					mainmodel(`mainmodel') ///
					mainlist(`mainlist') ///
					mainlist2(`mainlist2') ///
					outcomelist(`outcomelist')  ///
					exposurelist(`exposurelist') ///
					winitial(unadjusted, independent) ///
					y(`y') ///
					a(`a') ///
					equations(Y A main)  ///
					instruments(Y:`mainlist' `outcomelist') ///
				    instruments(A:`mainlist2' `exposurelist') ///
				    instruments(main: `mainmodel')    ///
					from(`initY' `initA' `initM') ///
					onestep hasderivatives `vce'	
		}
		
}



// Tag a specific scenario (ebe, logit, logit)

		local tag = 0
		if (("`method'"=="ebe") &  ("`olink'" == "logit") & ("`elink'" == "logit")) local tag = 1

// Get beta vector 
			
		mat `getcoefs' = e(b)			
		if  ("`method'"=="obe")  mat `getcoefs' = `by'	
		if  ("`method'"=="ebe") &  ("`olink'" == "logit") & ("`elink'" == "logit") mat `getcoefs' = `bA'
				
// Define b and V when using GMM 

if ("`method'" != "obe")  & (`tag'!=1)    {

		if "`mainmodel'" != "" local fvmm : word 1 of `mainmodel'
		else local fvmm "0"
		
		mat `coefs' = `getcoefs'[1, "xmain_0:_cons".."xmain_`fvmm':_cons"]			
		local k : word count `mainmodel'
 
		foreach v of local mainmodel {
		local mmr "`mmr' `a'`v'"
		}
		local conams "`a' `mmr'"
		mat `getVCE' = e(V)	
		mat `VCE' = `getVCE'["xmain_0:_cons".."xmain_`fvmm':_cons",  "xmain_0:_cons".."xmain_`fvmm':_cons"]		 
} 

// Define b and V when NOT using GMM 

if ("`method'" == "obe")    {
  		local nw : word count `mainlist'
		local n1 : word 1 of `mainlist'
		local n2 : word `nw' of `mainlist'
		
		local start "`y':`n1'"
		local end "`y':`n2'"

		mat `coefs' = `getcoefs'[1, "`start'".."`end'"]	
		mat `getVCE' = e(V)	
		mat `VCE' = `getVCE'["`start'".."`end'",  "`start'".."`end'"]		 
		
		local conams "`mainlist'"
} 

 
if (`tag' == 1)   {
  		local nw : word count `mainlist2'
		local n1 : word 1 of `mainlist2'
		local n2 : word `nw' of `mainlist2'
		
		if ("`olink'" == "identity") {
		local start `n1'
		local end  `n2'
		}
	    if inlist("`olink'", "log", "logit")==1 {
		local start "`a':`n1'"
		local end "`a':`n2'"
		}
			
		mat `coefs' = `getcoefs'[1, "`start'".."`end'"]	
		mat `getVCE' = e(V)	
		mat `VCE' = `getVCE'["`start'".."`end'",  "`start'".."`end'"]		 
		
		local conams "`mainlist'"
} 
	
		
		foreach v of local conams {
		local eqnams `"`eqnams' main"'
		}	
	
		mat rownames `VCE' = `conams'
		mat colnames `VCE' = `conams'
		mat colnames `coefs' = `conams'		
		mat roweq `VCE'  = `eqnams'
		mat coleq `VCE'   =  `eqnams'
		mat coleq `coefs'   = `eqnams'
		
		ereturn post `coefs' `VCE', obs(`nobs')   depn(`y') 
		restore  
		ereturn repost, esample(`touse')
		ereturn local cmdline `"drglm `cmdline'"'
		ereturn local cmd "drglm"
		ereturn local olink "`olink'"
		ereturn local elink "`elink'"
		ereturn local estimator "`method'"
		ereturn local vcetype "Robust"
		
		if "`clustvar'" != "" {
			ereturn local clustvar = "`clustvar'"
			ereturn local vce = "cluster"
		}	

		_post_vce_rank
		
	 di _n as txt "Generalized Linear Models" _col(54) _c
	 di in gr "Number of obs =" in ye %10.0g e(N)	 
	 di in gr "Estimator: "  _c
	 if "`method'" == "dr" di in ye  "Double Robust" 
	 if "`method'" == "obe" di in ye  "Outcome Based" 
	 if "`method'" == "ebe" di in ye  "Exposure Based"

	 di in gr "Link functions: Outcome["in y "`e(olink)'" in gr ///
	          "]  Exposure[" in y "`e(elink)'" in gr "]"  

// Display  results 
 
	 Replay , level(`level') `eform'

end

capture program drop Replay
program Replay
	syntax [, Level(cilevel) eform ]
	if "`eform'" != "" ereturn display, level(`level')  eform("exp(b)")
	else ereturn display, level(`level')
end
