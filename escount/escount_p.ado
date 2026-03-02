*! version August 2019 Takuya Hasebe ;
#delimit ;
* predict command for escount *;

program define escount_p ;
	syntax anything [if] [in] [, psel xb0 xb1 xbsel mu0 mu1 mu1d1 mu1d0 mu0d1 mu0d0 ];
	marksample touse;
	syntax newvarname [if] [in] [, psel xb0 xb1 xbsel mu0 mu1 mu1d1 mu1d0 mu0d1 mu0d0 ];
	
	if "`e(cmd)'" != "escount" error 301;
	local check "`psel' `xb0' `xb1' `xbsel' `mu0' `mu1' `cll' `mu1d1' `mu1d0' `mu0d1' `mu0d0'";
	if wordcount("`check'")>1 {;
		dis as error "Only one statistic is allowed.";
	};
	if wordcount("`check'")==0 {;
		dis as text "Option psel is assumed";
		local psel "psel";
	};
	
	local y0: word 1 of `e(depvar)';
	local y1: word 2 of `e(depvar)';
	local S: word 3 of `e(depvar)';
	local outoff `e(offset1)';
	local seloff `e(offset2)';
	
	tempname beta betas ;
	tempvar zg Fs ;
	
	matrix `beta' = e(b); 
	
	matrix `betas' = `beta'[1,"selection:"];	// coef vector of selection equation
	qui matrix score double `zg' = `betas' if `touse' ;
	if ~missing("`seloff'") qui replace	`zg' = `zg' + `seloff' if `touse' ;
	
	qui gen double `Fs' = normal(-`zg') if `touse';
	
	if ~missing("`xbsel'") {;
		generate `typlist' `varlist' = `zg' if `touse';	//linear probability;
		exit;
	};

	if ~missing("`psel'") {;
		generate `typlist' `varlist' = 1-`Fs' if `touse';
		exit ;
	};
	
		
	foreach j in 0 1 {;
		tempvar _xb`j' ;
		tempname b`j' sigma`j' rho`j';
		
		local margin`j' `e(margin`j')';
		matrix `b`j'' = `beta'[1,"status`j':"];
		scalar `sigma`j'' = exp(_b[/lnsigma`j']);
		scalar `rho`j'' = tanh(_b[/athrho`j']);
		
		
		qui matrix score double `_xb`j'' = `b`j'' if `touse' ;
		if ~missing("`outoff'") qui replace `_xb`j'' = `_xb`j''+`outoff' if `touse' ;
		
		if ~missing("`xb`j''") {;
			generate `typlist' `varlist' = `_xb`j'' if `touse' ; exit ;
		};
				
		if "`margin`j''"=="nb" {;
			tempname lnalpha`j' ;
			scalar `lnalpha`j'' = _b[/lnalpha`j'];
		};
		
		if ~missing("`mu`j''") {;	// only for Poisson and NB
			gen `typlist' `varlist' = exp(`_xb`j'')*exp(`sigma`j''^2/2) if `touse' ; exit ;
		};
		
		if ~missing("`pr`j''") {;
			if "`margin`j''"=="poisson" {;
				gen `typlist' = poissonp(`mu`j'',`pr`j'') if `touse';	exit; 
			}; 
			else if "`margin`j''"=="nb" {;
				
			};
		};
	};
	
	if ~missing("`mu0d0'") {;	// mean of status 0 conditional on being in status 0
		gen `typlist' `varlist' = exp(`_xb0')*exp(`sigma0'^2/2)*
			(normal(-(`rho0'*`sigma0'+`zg'))/normal(-`zg')) if `touse';
	};
	
	if ~missing("`mu0d1'") {;	// mean of status 0 conditional on being in status 1
		gen `typlist' `varlist' = exp(`_xb0')*exp(`sigma0'^2/2)*
			(normal(`rho0'*`sigma0'+`zg')/normal(`zg')) if `touse';
	};
	
	if ~missing("`mu1d0'") {;	// mean of status 1 conditional on being in status 0
		gen `typlist' `varlist' = exp(`_xb1')*exp(`sigma1'^2/2)*
			(normal(-(`rho1'*`sigma1'+`zg'))/normal(-`zg'))	if `touse';
	};
	
	if ~missing("`mu1d1'") {;	// mean of status 0 conditional on being in status 1
		gen `typlist' `varlist' = exp(`_xb1')*exp(`sigma1'^2/2)*
			(normal(`rho1'*`sigma1'+`zg')/normal(`zg'))	if `touse';
	};
	
	
	
	
end; 

