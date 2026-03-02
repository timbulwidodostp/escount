*! version 1.0.0 March, 2017 Takuya Hasebe
*** wrapper program to execute endogenously switching count data regression model ***;

#delimit ;

program define escount, sortpreserve ;
	if replay() {;
		if ("`e(cmd)'" != "escount") error 301 ;
		Replay `0';
	};
	else Estimate `0';

end;

program Estimate, eclass sortpreserve;
	syntax anything(id="equations" equalok) [aweight pweight iweight fweight] 
		[if] [in], SELect(string) 
		[  model(string)  INTPoints(integer 24) 
		CONSTraints(numlist) 
		NOCONStant EXPosure(varlist) OFFset(varlist)
		init_betas(string) init_beta0(string) init_beta1(string)
		init_lnalpha0(string) init_lnalpha1(string) 
		init_lnsigma0(string) init_lnsigma1(string) 
		init_athrho0(string) init_athrho1(string)
		init_lndfs(string)	*];
		
	mlopts mlopts, `options';	

	
	/** dealing with margins**/;
	local margout `model';
	foreach j in 0 1 {;
		if "`margout'" == "" {; local margin`j' "poisson"; };
		else {; local margin`j' `margout';};			
		
		if  "`margin`j''" != "poisson" & "`margin`j''" != "nb"   {; 
			dis as error "`margout' is not an available marginal distribution."; exit 198;
		};	
	};
		
	gettoken reg0 reg1_ : anything, match(parns) bind;
	if "`parns'" != "(" {;
		local reg0 "`anything'";
		local reg1_ "";
	};
	tokenize "`reg0'", parse("()");
	if "`1'" != "`reg0'" {;
		dis as error "If more than one, equations must be enclosed in brackets.";
		exit 198;
	};
	tokenize "`reg0'", parse("=");
	if "`2'" != "=" {;
		tokenize "`reg0'";
		local y0 `1';
		macro shift;
		local x0 `*';
	};
	else {;
		if "`4'" == "=" {;
			dis error "If more than one, equations must be enclosed in brackets.";
			exit 198;
		};
		local y0 `1';
		local x0 `3';
	};
	capture unab x0 : `x0';
	
	if "`reg1_'" != "" {;
		gettoken reg1 : reg1_, match(parns) bind;
		if "`parns'" != "(" {;
			local reg1 `reg1_';
		};
		tokenize "`reg1'", parse("=");
		if "`2'" != "=" {;
			tokenize "`reg1'";
			local y1 `1';
			macro shift;
			local x1 `*';
		};
		else {;
			if "`4'" == "=" {;
				dis as error "If more than one, equations must be enclosed in brackets.";
				exit 198;
			};
			local y1 `1';
			local x1 `3';
		};
	};
	else {;
		local y1 `y0'; local x1 `x0';
	};
	capture unab x1: `x1';
	
	*** selection equation ***;
	Select ys xs selnc seloff : `"`select'"';
	/**
	tokenize `select', parse("=");
	if "`2'" != "=" {;
		tokenize `select'; local ys `1'; macro shift; local xs `*';
	};
	else {;	local ys `1'; local xs `3'; };
	capture unab xs : `xs';
	**/;
	
	marksample touse;
	markout `touse' `y0' `x0' `y1' `x1' `ys' `xs' ;
		
	qui tab `ys' if `touse';
	if r(r) != 2 {; dis as error "`ys' should be binary."; exit 198; };
	qui sum `ys' if `touse';
	if r(max) != 1  & r(min) != 0 {;dis as error "`ys' should be either 0 or 1."; exit 198; };
	
	*** noconstant opition ***;
	if "`noconstant'"!="" {; local outnc noconstant; };
	
	*** checking collinearity ***;
	qui _rmdcoll `ys' `xs' if `touse', `selnc';
	local result "`r(varlist)'";
	local coll_x: list xs - result;
	if ~missing("`coll_x'") {;
		noisily display as text "note: `coll_x' omitted from selection equation because of colliearity";
		local xs `result';
	};
		
	forvalue j = 0/1 {;
		qui _rmdcoll `y`j'' `x`j'' if `touse', `outnc';
		local result "`r(varlist)'";
		local coll_x: list x`j' - result;
		if ~missing("`coll_x'") {;
			noisily display as text "note: `coll_x' omitted from outcome `j' because of colliearity";
			local x`j' `result';
		};
	};

	
	*** obtaining initial values ***;
	local eq = 3 ;
	tempname gamma g_init;
	
	if "`init_betas'"=="" {;
		dis _newline as text "fitting initial values from probit model:";
		probit `ys'  `xs' if `touse' [`weight'`exp'], `mlopts' `selnc'  offset(`seloff') nocoef; //`mlopts' 
		
		local ll_s = `e(ll)';	//log likelihood of selection equation;
			matrix `gamma' = e(b);
			matrix `g_init' = `gamma'[1,"`ys':"];	
	};
	else {;
		matrix `g_init' = `init_betas';
	};
	
	*** obtaining initial values of outcome parts ***;
	foreach j in 0 1 {;
		tempname beta`j' b`j'_init sig`j'_init;
		
		if "`init_beta`j''" == "" {;
			dis _newline as text "fitting initial values from lognormal `margin`j'' regression model:" ;
					
			lncount `y`j'' `x`j'' if `touse' & `ys'==`j' [`weight'`exp'], 
				intp(`intpoints') model(`margin`j'')  `mlopts' `outnc'  
				exposure(`exposure') offset(`offset') nocoef; //
		
			local ll_`j' = `e(ll)'; // log likelihood
			matrix `beta`j'' = e(b);	matrix `b`j'_init' = `beta`j''[1,"`y`j'':"];	
			matrix `sig`j'_init' = _b[/lnsigma];
			if "`margin`j''"=="nb" {;
				tempname a`j'_init ;	matrix `a`j'_init' = `beta`j''[1,"/lnalpha:"];
				local ml_a`j' /lnalpha`j';
			};
		};
		else {;
			matrix `b`j'_init' = `init_beta`j'';
			matrix `sig`j'_init' = `init_lnsigma`j'';
			if "`margin`j''"=="nb" {; 
				tempname a`j'_init ;
				local ml_a`j' /lnalpha`j';
				matrix `a`j'_init' = `init_lnalpha`j'';
			};
		};
	};	// end of foreach j in 0 1 
	*********************************************************************************************;
	* initial values for variance-covariance of errors *;
	tempname err_init;
	matrix `err_init' = (`sig0_init', `sig1_init', 0.0, 0.0);

	
	if (`intpoints'<1 | `intpoints'>512) {;
		dis as error 
			"intpoints() must be an integer greater than 1 and less than 513";
		exit 198;
	};	
	mata nnode = `intpoints';
		
	
	
	*** ML estimation *** ;
	dis _newline as text "Fitting full model:";
	
	if "`margin0'"=="poisson" {;
		local lfmata "ml model lf2 _espoisson_GHQ_lf()";
	};
	else {;
		local lfmata "ml model lf2 _esnb_GHQ_lf()";
	};
	
	`lfmata' (status0: `y0'=`x0', `outnc' exposure(`exposure') offset(`offset')) 
		(status1: `y1'=`x1', `outnc' exposure(`exposure') offset(`offset')) 
		(selection: `ys' = `xs', `selnc' offset(`seloff')) `ml_as' `ml_ds' `ml_a0' `ml_a1' /lnsigma0 /lnsigma1 /athrho0 /athrho1 
		[`weight'`exp']
		if `touse',
		title("endogenously switching `margin0' regression")
		missing collinear
		init(`b0_init' `b1_init' `g_init' `as_init' `ds_init' `a0_init' `a1_init' `err_init', copy) 
		maximize `vce' `mlopts' search(off) constraints(`constraints')
	;
	
	***********;
	* ereturn *;
	***********;
	* scalars *;
	*ereturn scalar ll0 = `ll_s'+`ll_0'+`ll_1';
	ereturn scalar n_quad = `intpoints';
		
	local k_aux = 4;
	if "`margout'"=="nb" local k_aux = `k_aux' + 2;
	
	ereturn scalar k_aux = `k_aux' ;
	
	* macros * ;
	ereturn local cmd escount ;
	ereturn local model `margout';
	ereturn local offset1 `offset';
	ereturn local offset2 `seloff';
	if ~missing("`exposure'") ereturn local offset1 "ln(`exposure')";
	
	* prediction *;
	ereturn local predict "escount_p";	
	ereturn local marginsok default psel xb0 xb1 xbsel mu0 mu1 mu1d1 mu1d0 mu0d1 mu0d0;
	
	Replay, level(`level');
	
end;

program Replay;
	syntax [, Level(cilevel)];
	
	foreach j in 0 1 {;
		local margin`j'  `e(model)' ;
		
		if "`margin`j''"=="nb" {;
			local diparm_a`j' diparm(lnalpha`j', exp label("alpha`j'"));
		};
		
		local diparm_s`j' diparm(lnsigma`j', exp label("sigma`j'"));
		local diparm_r`j' diparm(athrho`j', tanh label("rho`j'")) ;
		
	};
	
	ml display, level(`level') 
		 `diparm_a0' `diparm_a1' `diparm_s0' `diparm_s1' `diparm_r0' `diparm_r1'
		;
end;



program define Select ;
	args seldep selind selnc seloff colon sel_eqn ;

	gettoken dep rest : sel_eqn, parse(" =") ;
	gettoken equal rest : rest, parse(" =") ;

	if "`equal'" == "=" { ;
		tsunab dep : `dep' ;
		c_local `seldep' `dep' ;
	}; 
	else local rest `"`sel_eqn'"';
	
	local 0 `"`rest'"';
	syntax [varlist(numeric default=none)] [, noCONstant OFFset(varname numeric) ] ;

	if "`varlist'" == "" {;
		di in red "no variables specified for selection equation" ;
		exit 198 ;
	};

	c_local `selind' `varlist' ;
	c_local `selnc' `constant' ;
	c_local `seloff' `offset' ;
end ;

