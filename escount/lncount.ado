*! version 1.0.0 March, 2017 Takuya Hasebe
*** wrapper program to execute log-normal count data regression model ***;

#delimit ;

program define lncount, sortpreserve ;
	if replay() {;
		if ("`e(cmd)'" != "lncount") error 301 ;
		Replay `0';
	};
	else Estimate `0';

end;

program Estimate, eclass sortpreserve;
	syntax anything(id="equations" equalok) [aweight pweight iweight fweight] 
		[if] [in], [ model(string) INTPoints(integer 24) 
			NOCOEf NOCONStant EXPosure(varlist) OFFset(varlist) *];
		
	mlopts mlopts, `options';	

	*** noconstant opition ***;
	if "`noconstant'"!="" {; local outnc noconstant; };
	
	/** dealing with margins**/;
	local margout `model';
	if "`margout'" == "" {; local margout "poisson"; };
	if  "`margout'" != "poisson" & "`margout'" != "nb"   {; 
		dis as error "`margout' is not an available marginal distribution."; exit 198;
	};	
	
	tokenize `anything', parse("=");
	if "`2'" != "=" {;
		tokenize `anything';
		local y `1'; macro shift; local x `*';
	};
	else {;
		local y `1'; local x `3';
	};
	capture unab x : `x';
	
	marksample touse;
	markout `touse' `y' `x';
	

	*** checking collinearity ***;
	qui _rmdcoll `y' `x' if `touse', `outnc';
	local result "`r(varlist)'";
	local coll_x: list x - result;
	if ~missing("`coll_x'") {;
		noisily display as text "note: `coll_x' omitted because of collinearity";
		local x `result';
	};
	
	
	*** obtaining initial values ***;
	local eq = 3 ;
	tempname gamma b_init beta;
	
		
	*** obtaining initial values of outcome parts ***;
	if "`margout'"=="poisson" {;
		qui poisson `y' `x' if `touse' [`weight'`exp'], `outnc' exposure(`exposure') offset(`offset');
		matrix `b_init' = e(b);
	};
	else if "`margout'"=="nb" {;
		qui nbreg `y' `x' if  `touse' [`weight'`exp'], `outnc' exposure(`exposure') offset(`offset');
		matrix `beta' = e(b);
		matrix `b_init' = `beta'[1,"`y':"];
		tempname a_init ;	matrix `a_init' = `beta'[1,"/lnalpha:"];
		local ml_a /lnalpha;
	};	
	*********************************************************************************************;
	
	* initial values for variance of errors *;
	tempvar err_init ;
	matrix `err_init' = (0);
	
	*** ML estimation *** ;
	mata nnode = `intpoints';
	
	if "`margout'"=="poisson" local lf "ml model lf2 _lnpoisson_GHQ_lf()";
	else if "`margout'"=="nb" local lf "ml model lf2 _lnnb_GHQ_lf()";
		
	`lf' (`y': `y'=`x', `outnc' exposure(`exposure') offset(`offset')) `ml_a' /lnsigma 
		[`weight'`exp']
		if `touse',
		title("log-normal `margout' regression")
		missing collinear
		init(`b_init' `a_init' `err_init', copy) 
		maximize `vce' `mlopts' search(off)
	;
	
	***********;
	* ereturn *;
	***********;
	* scalars *;
	ereturn scalar n_quad = `intpoints';
		
	local k_aux = 1;
	if "`margout'"=="nb" local k_aux = `k_aux' + 1;
	ereturn scalar k_aux = `k_aux' ;
	
	* locals * ;
	ereturn local cmd lncount ;
	ereturn local model `margout';
	ereturn local offset `offset';
	if ~missing("`exposure'") ereturn local offset "ln(`exposure')";
	
	* prediction *;
	ereturn local predict "lncount_p";	
	ereturn local marginsok default n xb;
	if missing("`nocoef'") Replay, level(`level');
	
end;

program Replay;
	syntax [, Level(cilevel)];
	
	local margin  `e(model)' ;
	
	if "`margin'"=="nb" {;
		local diparm_a diparm(lnalpha, exp label("alpha"));
	};
	
	local diparm_s diparm(lnsigma, exp label("sigma"));
	
	
	ml display, level(`level') 
		`diparm_a' `diparm_s' 
		;
	*dis in smcl "{hline 78}";
end;


