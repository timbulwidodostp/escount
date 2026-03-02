*! version 1.0.0 March 2016 Takuya Hasebe ;
#delimit ;
* predict command for lncount *;

program define lncount_p ;
	syntax anything [if] [in] [, xb n];
	marksample touse;
	syntax newvarname [if] [in] [, xb n];
	
	if "`e(cmd)'" != "lncount" error 301;
	local check "`xb' `n' ";
	if wordcount("`check'")>1 {;
		dis as error "Only one statistic is allowed.";
	};
	if wordcount("`check'")==0 {;
		dis as text "Option n is assumed";
		local n "n";
	};
	
	local y: word 1 of `e(depvar)';
	local outoff `e(offset)';
	
	tempname beta ;
	matrix `beta' = e(b); 
		
	local margin `e(model)';
	
	tempvar _xb ;
	tempname b sigma ;
	
	matrix `b' = `beta'[1,"`y':"];
	scalar `sigma' = exp(_b[/lnsigma]);
	
	qui matrix score double `_xb' = `b' if `touse' ;
	if ~missing("`outoff'") qui replace `_xb' = `_xb'+`outoff' if `touse' ;
		
	if ~missing("`xb'") {;
		generate `typlist' `varlist' = `_xb' if `touse' ; exit ;
	};
		
	if "`margin'"=="nb" {;
		tempname lnalpha ;
		scalar `lnalpha' = _b[/lnalpha];
	};
	
	if ~missing("`n'") {;	// 
		gen `typlist' `varlist' = exp(`_xb')*exp(`sigma'^2/2) if `touse' ; exit ;
	};
	
	
	
end; 

