#delimit ;
*** program to calculate the treatment effects after switch count model ;
program define teescount, rclass sortpreserve ;
	syntax  [if] [in] , [ATE ATT ATU LATE(string) MTE(string)
		]	//iv(string) ZLower(string) ZUpper(string) at(string)
		;	//
	
	if "`e(cmd)'" != "escount" error 301;

	local k  = wordcount("`ate' `att' `atu' `mte'");
	if ~missing("`late'") local k = `k'+1;
	if `k'==0 {;
		dis as text "Option ate is assumed";
		local ate "ate"; local k = `k'+1;
	};
	
	if ~missing("`ate'`att'`atu'") {;
		treatment_escount;
		
		if ~missing("`ate'") {;
			tempname bate ;
			matrix `bate' = J(1,4,.);	matrix rownames `bate' = "ATE";
			matrix `bate'[1,1] = r(ate);	matrix `bate'[1,2] =sqrt(r(Vate));
			matrix `bate'[1,3] = r(ate)/sqrt(r(Vate));
			matrix `bate'[1,4] = 2*(1-normal(abs(r(ate)/sqrt(r(Vate)))));
		};	
		if ~missing("`att'") {;
			tempname batt ;
			matrix `batt' = J(1,4,.);	matrix rownames `batt' = "ATT";
			matrix `batt'[1,1] = r(att);	matrix `batt'[1,2] =sqrt(r(Vatt));
			matrix `batt'[1,3] = r(att)/sqrt(r(Vatt));
			matrix `batt'[1,4] = 2*(1-normal(abs(r(att)/sqrt(r(Vatt)))));
		};	
		if ~missing("`atu'") {;
			tempname batu ;
			matrix `batu' = J(1,4,.);	matrix rownames `batu' = "ATU";
			matrix `batu'[1,1] = r(atu);	matrix `batu'[1,2] =sqrt(r(Vatu));
			matrix `batu'[1,3] = r(atu)/sqrt(r(Vatu));
			matrix `batu'[1,4] = 2*(1-normal(abs(r(atu)/sqrt(r(Vatu)))));
		};			
	};
	if ~missing("`late'") {;
		tokenize `late';
		local iv `1';	local zlower = `2'; local zupper = `3';
		qui late_escount, iv(`iv') zlower(`zlower') zupper(`zupper');
		
		tempname blate;
		matrix `blate' = J(1,4,.);	matrix rownames `blate' = "LATE";
		matrix `blate'[1,1] = r(late);		matrix `blate'[1,2] = sqrt(r(Vlate));
		matrix `blate'[1,3] = r(late)/sqrt(r(Vlate));
		matrix `blate'[1,4] = 2*(1-normal(abs(r(late)/sqrt(r(Vlate)))));
	};
	if ~missing("`mte'") {;
		qui mte_escount, at(`mte');
		
		tempname bmte;
		matrix `bmte' = J(1,4,.);	matrix rownames `bmte' = "MTE";
		matrix `bmte'[1,1] = r(mte);		matrix `bmte'[1,2] = sqrt(r(Vmte));
		matrix `bmte'[1,3] = r(mte)/sqrt(r(Vmte));
		matrix `bmte'[1,4] = 2*(1-normal(abs(r(mte)/sqrt(r(Vmte)))));
	};
	
	
	tempname result; 
	
	tokenize "`bate' `batt' `batu' `blate' `bmte'";
	forvalue t = 1/`k' {;
		if `t' == 1 {; matrix `result' = ``t'';};
		else matrix `result' = (`result' \ ``t'');
	};
	
	matrix colnames `result' = Estimate "SE" "z" "P>|z|" ; 
	dis in green "Treatment Effects:";
    matrix list `result', noheader format(%7.5g);
	
	
	return matrix table = `result';
			
end;


	
#delimit ;
*** program to calculate the treatment effects after switch count model ;
*** ATE, ATT, ATU ***;
program define treatment_escount, rclass sortpreserve ;
	syntax  [if] [in] , 
		;	//
	
	if "`e(cmd)'" != "escount" error 301;
	
	tempname beta betas V N Ps Vate Vatt Vatu Vdate Vdatt Vdatu ate att atu S22 G21 S12;
	tempvar zg _mu0 _dmu0db _dmu0da _dmu0dg _dmu0ds _dmu0dr _mu1 _dmu1db _dmu1da _dmu1dg _dmu1ds _dmu1dr ;
	tempvar	touse te tte tue cons;
	
	qui gen `touse' = e(sample) ;
	matrix `beta' = e(b);	matrix `V' = e(V);
		
	local y0: word 1 of `e(depvar)';
	local y1: word 2 of `e(depvar)';
	local s: word 3 of `e(depvar)';
	local outoff `e(offset1)';
	local seloff `e(offset2)';
	local cluster `e(clustvar)';
	
	local _cons _cons;
	qui gen `cons' = _cons ;
		
	matrix `betas' = `beta'[1,"selection:"];	// coef vector of selection equatio
	matrix score double `zg' = `betas' if `touse' ;
	if ~missing("`seloff'") qui replace	`zg' = `zg' + `seloff' if `touse' ;
	
	local xs_temp : colnames `betas' ; 
	local xs: list  xs_temp - _cons;
	if "`xs_temp'"!="`xs'" {; local xs `xs' `cons';};
	mata st_view(Xs=.,.,"`xs'","`touse'");
	
	tempvar Fs nu dnudg dFsdg;
	qui gen double `Fs' = normal(-`zg');
	qui gen double `nu' = invnormal(`Fs');
		
	qui gen double `dFsdg' = -normalden(-`zg');
	qui gen double `dnudg' = `dFsdg'/normalden(`nu');
	
	qui sum `touse' if `touse', meanonly; 	scalar `N' = r(N); 
	qui sum `s' if `touse', meanonly;	scalar `Ps' = r(mean);	
	
	foreach j in 0 1 {;
		local margin`j' `e(model)';
		
		tempvar xb`j' mu`j' dmu`j'db dmu`j'da dmu`j'ds dmu`j'dr;
		tempvar mu`j's0 dmu`j's0db dmu`j's0ds dmu`j's0dr dmu`j's0dg dmu`j's0da dmu`j's0dd ;
		tempvar mu`j's1 dmu`j's1db dmu`j's1ds dmu`j's1dr dmu`j's1dg dmu`j's1da dmu`j's1dd;
		tempname beta`j' alpha`j' lnalpha`j' athrho`j' rho`j' lnsigma`j' sigma`j' dATEdb`j' dATEds`j' ;
		tempname dATTdb`j' dATTds`j' dATTda`j' dATTdr`j' dATTdg`j' dATTda`j' dATTdd`j';
		tempname dATUdb`j' dATUds`j' dATUda`j' dATUdr`j' dATUdg`j' dATUda`j' dATUdd`j';
		
		matrix `beta`j'' = `beta'[1,"status`j':"];
		qui matrix score double `xb`j'' = `beta`j'' if `touse' ;
		if ~missing("`outoff'") qui replace `xb`j'' = `xb`j''+`outoff' if `touse' ;
		
		local x`j'_temp : colnames `beta`j'' ; local x`j': list  x`j'_temp - _cons;
		if "`x`j'_temp'"!="`x`j''" {; local x`j' `x`j'' `cons';};
		mata st_view(X`j'=.,.,"`x`j''","`touse'");
	
		scalar `lnsigma`j'' = _b[/lnsigma`j'];	scalar `sigma`j'' = exp(`lnsigma`j'');
		scalar `athrho`j'' = _b[/athrho`j'];	scalar `rho`j''= tanh(`athrho`j'');
		
		if "`margin`j''"=="nb" {;
			scalar `lnalpha`j'' = _b[/lnalpha`j'];
		};
		else {; scalar `lnalpha`j'' = 0; };
				
	
		qui gen double `mu`j'' = exp(`xb`j'')*exp(`sigma`j''^2/2);
		qui gen double `mu`j's1' = `mu`j''*(normal(`rho`j''*`sigma`j''-`nu')/(1-`Fs'));
		qui gen double `mu`j's0' = `mu`j''*(normal(-(`rho`j''*`sigma`j''-`nu'))/`Fs');
		
		qui gen double `dmu`j'db' = `mu`j'';
		qui gen double `dmu`j'ds' = `mu`j''*(`sigma`j'')*`sigma`j'';
		
		qui gen double `dmu`j's1db' = `mu`j's1';
		qui gen double `dmu`j's1ds' = `mu`j's1'*`sigma`j''
			+ `mu`j''*(normalden(`rho`j''*`sigma`j''-`nu')/(1-`Fs'))*`rho`j'';
		qui replace `dmu`j's1ds' = `dmu`j's1ds'*`sigma`j'';
		qui gen double `dmu`j's1dg' = `mu`j''*(-normalden(`rho`j''*`sigma`j''-`nu')*`dnudg'/(1-`Fs')
			+normal(`rho`j''*`sigma`j''-`nu')*`dFsdg'/((1-`Fs')^2));
		cap gen double `dmu`j's1da' = `mu`j''*(-normalden(`rho`j''*`sigma`j''-`nu')*`dnuda'/(1-`Fs')
			+normal(`rho`j''*`sigma`j''-`nu')*`dFsda'/((1-`Fs')^2));
		cap gen double `dmu`j's1dd' = `mu`j''*(-normalden(`rho`j''*`sigma`j''-`nu')*`dnudd'/(1-`Fs')
			+normal(`rho`j''*`sigma`j''-`nu')*`dFsdd'/((1-`Fs')^2));	
		qui gen double `dmu`j's1dr' = `mu`j''*(normalden(`rho`j''*`sigma`j''-`nu')/(1-`Fs'))*`sigma`j'';
		qui replace `dmu`j's1dr' = `dmu`j's1dr'*(sqrt(1-(`rho`j'')^2))^2 ;
		
		qui replace `dmu`j's1db' = `s'*`dmu`j's1db'/`Ps';
		qui replace `dmu`j's1ds' = `s'*`dmu`j's1ds'/`Ps';
		qui replace `dmu`j's1dg' = `s'*`dmu`j's1dg'/`Ps';
		cap replace `dmu`j's1da' = `s'*`dmu`j's1da'/`Ps';
		cap replace `dmu`j's1dd' = `s'*`dmu`j's1dd'/`Ps';
		qui replace `dmu`j's1dr' = `s'*`dmu`j's1dr'/`Ps';
		
		qui gen double `dmu`j's0db' = `mu`j's0';
		qui gen double `dmu`j's0ds' = `mu`j's0'*`sigma`j''
			+ `mu`j''*(normalden(-(`rho`j''*`sigma`j''-`nu'))/`Fs')*(-`rho`j'');
		qui replace `dmu`j's0ds' = `dmu`j's0ds'*`sigma`j'';
		qui gen double `dmu`j's0dg' = `mu`j''*(normalden(-(`rho`j''*`sigma`j''-`nu'))*`dnudg'/`Fs'
			-normal(-(`rho`j''*`sigma`j''-`nu'))*`dFsdg'/(`Fs'^2));	
		cap gen double `dmu`j's0da' = `mu`j''*(normalden(-(`rho`j''*`sigma`j''-`nu'))*`dnuda'/`Fs'
			-normal(-(`rho`j''*`sigma`j''-`nu'))*`dFsda'/(`Fs'^2));	
		cap gen double `dmu`j's0dd' = `mu`j''*(normalden(-(`rho`j''*`sigma`j''-`nu'))*`dnudd'/`Fs'
			-normal(-(`rho`j''*`sigma`j''-`nu'))*`dFsdd'/(`Fs'^2));		
		qui gen double `dmu`j's0dr' = `mu`j''*(normalden(-(`rho`j''*`sigma`j''-`nu'))/`Fs')*(-`sigma`j'');
		qui replace `dmu`j's0dr' = `dmu`j's0dr'*(sqrt(1-(`rho`j'')^2))^2;
		
		qui replace `dmu`j's0db' = (1-`s')*`dmu`j's0db'/(1-`Ps');
		qui replace `dmu`j's0ds' = (1-`s')*`dmu`j's0ds'/(1-`Ps');
		qui replace `dmu`j's0dg' = (1-`s')*`dmu`j's0dg'/(1-`Ps');
		cap replace `dmu`j's0da' = (1-`s')*`dmu`j's0da'/(1-`Ps');
		cap replace `dmu`j's0dd' = (1-`s')*`dmu`j's0dd'/(1-`Ps');
		qui replace `dmu`j's0dr' = (1-`s')*`dmu`j's0dr'/(1-`Ps');
		
	};	// end of foreach j in 0 1 ;
		
	qui gen double `te' = `mu1' - `mu0'  if  `touse' ;
	qui gen double `tte' = `s'*(`mu1s1'-`mu0s1')/`Ps' if `touse' ;
	qui gen double `tue' = (1-`s')*(`mu1s0' - `mu0s0')/(1-`Ps') if `touse' ;
	
	
	******************************************************************************************************;
	qui sum `te' if `touse', meanonly ;		scalar `ate' = r(mean);		qui replace `te' = `te' - `ate';	
	qui sum `tte' if `touse', meanonly ;	scalar `att' = r(mean);		qui replace `tte' = `tte' - `att';	
	qui sum `tue' if `touse', meanonly ;	scalar `atu' = r(mean);		qui replace `tue' = `tue' - `atu';		
	
	mata st_view(hte=.,.,"`te'","`touse'");
	mata st_view(htt=.,.,"`tte'","`touse'");
	mata st_view(htu=.,.,"`tue'","`touse'");
	
	*** average treatment effect *** ;
	foreach j in 0 1 {;
		mata st_view(dmu`j'db=.,.,tokens("`dmu`j'db'"),"`touse'");
		mata st_matrix("`dATEdb`j''",quadcross(X`j',dmu`j'db)) ;
		
		if "`margin`j''"=="nb" {;
			local _dATEda`j' "\ 0" ; 
		};
		mata st_view(dmu`j'ds=.,.,tokens("`dmu`j'ds'"),"`touse'");
		mata st_matrix("`dATEds`j''",quadcolsum(dmu`j'ds)) ;
		local _dATEds`j' "\ (2*`j'-1)*`dATEds`j''" ;
		
		local _dATEdr`j' "\ 0" ;
	};
	
	if missing("`cluster'") {;
		mata st_matrix("`S22'",quadcross(hte,hte)) ;
	};
	else {;
		sort `cluster';
		mata V = _cluster_var("`cluster'", "`te'", "`touse'");
		mata st_matrix("`S22'",V); 
	};
	matrix `G21' = (-`dATEdb0' \ `dATEdb1' \ J(colsof(`betas'),1,0) `_dATEdas' `_dATEdds' `_dATEda0' `_dATEda1' `_dATEds0' `_dATEds1' `_dATEdr0' `_dATEdr1')/`N'; 
	matrix `Vate' = (`S22'/`N' + `G21''*(`N'*`V')*`G21')/`N' ;
	matrix `Vdate' = (`G21''*(`N'*`V')*`G21')/`N' ;
	
	*** average treatment effect on treated *** ;
	foreach j in 0 1 {;
		mata st_view(dmu`j's1db=.,.,tokens("`dmu`j's1db'"),"`touse'");
		mata st_matrix("`dATTdb`j''",quadcross(X`j',dmu`j's1db)) ;
		
		mata st_view(dmu`j's1dg=.,.,tokens("`dmu`j's1dg'"),"`touse'");
		mata st_matrix("`dATTdg`j''",quadcross(Xs,dmu`j's1dg)) ;
			
		if "`margin`j''"=="nb" {;
			local _dATTda`j' "\ 0" ; 
		};
						
		mata st_view(dmu`j's1ds=.,.,tokens("`dmu`j's1ds'"),"`touse'");
		mata st_matrix("`dATTds`j''",quadcolsum(dmu`j's1ds)) ;	
		local _dATTds`j' "\ (2*`j'-1)*`dATTds`j''" ;

		mata st_view(dmu`j's1dr=.,.,tokens("`dmu`j's1dr'"),"`touse'");
		mata st_matrix("`dATTdr`j''",quadcolsum(dmu`j's1dr)) ;
		local _dATTdr`j' "\ (2*`j'-1)*`dATTdr`j''" ;

	};
	
	if missing("`cluster'") {;
		mata st_matrix("`S22'",quadcross(htt,htt)) ;
	};
	else {;
		sort `cluster';
		mata V = _cluster_var("`cluster'", "`tte'", "`touse'");
		mata st_matrix("`S22'",V); 
	};
	matrix `G21' = (-`dATTdb0' \ `dATTdb1' \ (`dATTdg1'-`dATTdg0') `_dATTdas' `_dATTdds' `_dATTda0' `_dATTda1' `_dATTds0' `_dATTds1'  `_dATTdr0' `_dATTdr1')/`N'; 
	matrix `Vatt' = (`S22'/`N' + `G21''*(`N'*`V')*`G21')/`N'  ;
	matrix `Vdatt' = (`G21''*(`N'*`V')*`G21')/`N' ;
	
	
	*** average treatment effect on untreated *** ;
	foreach j in 0 1 {;
		mata st_view(dmu`j's0db=.,.,tokens("`dmu`j's0db'"),"`touse'");
		mata st_matrix("`dATUdb`j''",quadcross(X`j',dmu`j's0db)) ;
		
		mata st_view(dmu`j's0dg=.,.,tokens("`dmu`j's0dg'"),"`touse'");
		mata st_matrix("`dATUdg`j''",quadcross(Xs,dmu`j's0dg)) ;
		
		if "`margin`j''"=="nb" {;
			local _dATUda`j' "\ 0" ; 
		};
						
		mata st_view(dmu`j's0ds=.,.,tokens("`dmu`j's0ds'"),"`touse'");
		mata st_matrix("`dATUds`j''",quadcolsum(dmu`j's0ds)) ;	
		local _dATUds`j' "\ (2*`j'-1)*`dATUds`j''" ;

		mata st_view(dmu`j's0dr=.,.,tokens("`dmu`j's0dr'"),"`touse'");
		mata st_matrix("`dATUdr`j''",quadcolsum(dmu`j's0dr)) ;
		local _dATUdr`j' "\ (2*`j'-1)*`dATUdr`j''" ;

	};
	
	if missing("`cluster'") {;
		mata st_matrix("`S22'",quadcross(htu,htu)) ;
	};
	else {;
		sort `cluster';
		mata V = _cluster_var("`cluster'", "`tue'", "`touse'");
		mata st_matrix("`S22'",V); 
	};
	matrix `G21' = (-`dATUdb0' \ `dATUdb1' \ (`dATUdg1'-`dATUdg0') `_dATUdas' `_dATUdds' `_dATUda0' `_dATUda1' `_dATUds0' `_dATUds1'  `_dATUdr0' `_dATUdr1')/`N'; 
	matrix `Vatu' = (`S22'/`N' + `G21''*(`N'*`V')*`G21')/`N'  ;
	matrix `Vdatu' = (`G21''*(`N'*`V')*`G21')/`N' ;
	
*************************************************************************************************;
	return scalar Vdatu = `Vdatu'[1,1] ;
	return scalar Vatu = `Vatu'[1,1] ;
	return scalar atu = `atu' ;
		
	return scalar Vdatt = `Vdatt'[1,1] ;		
	return scalar Vatt = `Vatt'[1,1] ;
	return scalar att = `att' ;
	
	return scalar Vdate = `Vdate'[1,1] ;
	return scalar Vate = `Vate'[1,1] ;
	return scalar ate = `ate' ;
	
	mata: mata drop Xs X0 X1 hte htt htu dmu0db dmu1db dmu0s1db dmu1s1db dmu0s0db dmu1s0db dmu0s1dg dmu1s1dg dmu0s0dg dmu1s0dg; 
	cap mata: mata drop dmu0s1dr dmu1s1dr dmu0s0dr dmu1s0dr;
	cap mata: mata drop dmu0ds dmu1ds dmu0s1ds dmu1s1ds dmu0s0ds dmu1s0ds;
	cap mata: mata drop dmu0da dmu1da dmu0s1da dmu1s1da dmu0s0da dmu1s0da;
end;

*** LATE ***;
program define late_escount, rclass sortpreserve ;
	syntax  [if] [in] , iv(varname) ZLower(real) ZUpper(real)
		;	//
	
	if "`e(cmd)'" != "escount" error 301;
	
	tempname beta betas V N LATE Vlate Vdlate S22 G21;
	tempvar zgu zgl z_temp _mul _dmuldb _dmulda _dmulds _dmuldg _dmuldr _muu _dmuudb _dmuuda _dmuuds _dmuudg _dmuudr ;
	tempvar	touse lte cons;
		
	qui gen `touse' = e(sample) ;
	matrix `beta' = e(b);	matrix `V' = e(V);
		
	local y0: word 1 of `e(depvar)';
	local y1: word 2 of `e(depvar)';
	local s: word 3 of `e(depvar)';
	local outoff `e(offset1)';
	local seloff `e(offset2)';
	local cluster `e(clustvar)';
	
	local _cons _cons;
	qui gen `cons' = _cons ;
	
	matrix `betas' = `beta'[1,"selection:"];	// coef vector of selection equation
	qui gen `z_temp' = `iv' ;
	local xs_temp : colnames `betas' ; 
	local xs: list  xs_temp - _cons;
	if "`xs_temp'"!="`xs'" {; local xs `xs' `cons';};
	
	qui replace `iv' = `zupper' if `touse' ;
	matrix score double `zgu' = `betas' if `touse' ;
	if ~missing("`seloff'") qui replace	`zg' = `zgu' + `seloff' if `touse' ;
	mata Xsu = st_data(.,"`xs'","`touse'");
	qui replace `iv' = `zlower' if `touse' ;
	matrix score double `zgl' = `betas' if `touse' ;
	if ~missing("`seloff'") qui replace	`zgl' = `zgl' + `seloff' if `touse' ;
	mata Xsl = st_data(.,"`xs'","`touse'");
	
	
	
	qui sum `touse' if `touse', meanonly; 	scalar `N' = r(N); 
	
	foreach j in 0 1 {;
		local margin`j' `e(model)';
		
		tempvar xb`j' mu`j' dmu`j'db dmu`j'dgu dmu`j'dgl dmu`j'dau dmu`j'dal dmu`j'ddu dmu`j'ddl dmu`j'da dmu`j'ds dmu`j'dr rzgu`j' rzgl`j';
		tempname beta`j' alpha`j' lnalpha`j' lnsigma`j' sigma`j' athrho`j' rho`j' dLATEdb`j' dLATEdgu`j' dLATEdgl`j' dLATEdau`j' dLATEdal`j' dLATEddu`j' dLATEddl`j' dLATEds`j' dLATEda`j' dLATEdr`j' ;
		
		matrix `beta`j'' = `beta'[1,"status`j':"];
		qui matrix score double `xb`j'' = `beta`j'' if `touse' ;
		if ~missing("`outoff'") qui replace `xb`j'' = `xb`j''+`outoff' if `touse' ;
		local x`j'_temp : colnames `beta`j'' ; local x`j': list  x`j'_temp - _cons;
		if "`x`j'_temp'"!="`x`j''" {; local x`j' `x`j'' `cons';};
		mata st_view(X`j'=.,.,"`x`j''","`touse'");
		
		scalar `lnsigma`j'' = _b[/lnsigma`j'];	scalar `sigma`j'' = exp(`lnsigma`j'');
		scalar `athrho`j'' = _b[/athrho`j'];	scalar `rho`j''= tanh(`athrho`j'');
		
		foreach k in u l {;
			tempvar Fs`k' nu`k' dnu`k'dg dFs`k'dg;
			qui gen double `Fs`k'' = normal(-`zg`k'');
			qui gen double `nu`k'' = invnormal(`Fs`k'');
			
			qui gen double `dFs`k'dg' = -normalden(-`zg`k'');
			qui gen double `dnu`k'dg' = `dFs`k'dg'/normalden(`nu`k'');
		};
		
		if "`margin`j''"=="nb" {;
			scalar `lnalpha`j'' = _b[/lnalpha`j'];
		};
		else {; scalar `lnalpha`j''= 0; };
		
		qui gen double `rzgu`j'' = -`rho`j''*`sigma`j''+`nuu';
		qui gen double `rzgl`j'' = -`rho`j''*`sigma`j''+`nul';
		
		qui gen double `mu`j'' = exp(`xb`j'')*exp(`sigma`j''^2/2)*(normal(`rzgl`j'')-normal(`rzgu`j''))/(`Fsl'-`Fsu');		
		
		qui gen double `dmu`j'db' = `mu`j'';
		qui gen double `dmu`j'dgu' = exp(`xb`j'')*exp(`sigma`j''^2/2)*(-normalden(`rzgu`j'')*`dnuudg'/(`Fsl'-`Fsu')
			+(normal(`rzgl`j'')-normal(`rzgu`j''))*`dFsudg'/((`Fsl'-`Fsu')^2));
		qui gen double `dmu`j'dgl' = exp(`xb`j'')*exp(`sigma`j''^2/2)*(normalden(`rzgl`j'')*`dnuldg'/(`Fsl'-`Fsu')
			-(normal(`rzgl`j'')-normal(`rzgu`j''))*`dFsldg'/((`Fsl'-`Fsu')^2));
		cap gen double `dmu`j'dau' = exp(`xb`j'')*exp(`sigma`j''^2/2)*(-normalden(`rzgu`j'')*`dnuuda'/(`Fsl'-`Fsu')
			+(normal(`rzgl`j'')-normal(`rzgu`j''))*`dFsuda'/((`Fsl'-`Fsu')^2));
		cap gen double `dmu`j'dal' = exp(`xb`j'')*exp(`sigma`j''^2/2)*(normalden(`rzgl`j'')*`dnulda'/(`Fsl'-`Fsu')
			-(normal(`rzgl`j'')-normal(`rzgu`j''))*`dFslda'/((`Fsl'-`Fsu')^2));
		cap gen double `dmu`j'ddu' = exp(`xb`j'')*exp(`sigma`j''^2/2)*(-normalden(`rzgu`j'')*`dnuudd'/(`Fsl'-`Fsu')
			+(normal(`rzgl`j'')-normal(`rzgu`j''))*`dFsudd'/((`Fsl'-`Fsu')^2));
		cap gen double `dmu`j'ddl' = exp(`xb`j'')*exp(`sigma`j''^2/2)*(normalden(`rzgl`j'')*`dnuldd'/(`Fsl'-`Fsu')
			-(normal(`rzgl`j'')-normal(`rzgu`j''))*`dFsldd'/((`Fsl'-`Fsu')^2));
			
		qui gen double `dmu`j'ds' = `mu`j''*`sigma`j''
			-`rho`j''*exp(`xb`j'')*exp(`sigma`j''^2/2)*(normalden(`rzgl`j'')-normalden(`rzgu`j''))/(`Fsl'-`Fsu');
		qui replace `dmu`j'ds' = `sigma`j''*`dmu`j'ds' ;
		qui gen double `dmu`j'dr' = -`sigma`j''*exp(`xb`j'')*exp(`sigma`j''^2/2)*(normalden(`rzgl`j'')-normalden(`rzgu`j''))/(`Fsl'-`Fsu');
		qui replace `dmu`j'dr' = `dmu`j'dr'*(sqrt(1-(`rho`j'')^2))^2;
		
	};	// end of foreach j in 0 1 ;
	
	qui gen double `lte' = `mu1' - `mu0'  if  `touse' ;
	******************************************************************************************************;
	qui sum `lte' if `touse', meanonly ;	scalar `LATE' = r(mean);	qui replace `lte' = `lte' - `LATE';		
	
	mata st_view(lte=.,.,"`lte'","`touse'");
	
	*** local average treatment effect *** ;
	foreach j in 0 1 {;
		mata st_view(dmu`j'db=.,.,tokens("`dmu`j'db'"),"`touse'");
		mata st_matrix("`dLATEdb`j''",quadcross(X`j',dmu`j'db)) ;
		
		mata st_view(dmu`j'dgu=.,.,tokens("`dmu`j'dgu'"),"`touse'");
		mata st_matrix("`dLATEdgu`j''",quadcross(Xsu,dmu`j'dgu)) ;	
		
		mata st_view(dmu`j'dgl=.,.,tokens("`dmu`j'dgl'"),"`touse'");
		mata st_matrix("`dLATEdgl`j''",quadcross(Xsl,dmu`j'dgl)) ;
		
		if "`margsel'"=="GTL" {;
			mata st_view(dmu`j'dau=.,.,tokens("`dmu`j'dau'"),"`touse'");
			mata st_matrix("`dLATEdau`j''",quadcolsum(dmu`j'dau)) ;	
		
			mata st_view(dmu`j'dal=.,.,tokens("`dmu`j'dal'"),"`touse'");
			mata st_matrix("`dLATEdal`j''",quadcolsum(dmu`j'dal)) ;
			
			mata st_view(dmu`j'ddu=.,.,tokens("`dmu`j'ddu'"),"`touse'");
			mata st_matrix("`dLATEddu`j''",quadcolsum(dmu`j'ddu)) ;	
		
			mata st_view(dmu`j'ddl=.,.,tokens("`dmu`j'ddl'"),"`touse'");
			mata st_matrix("`dLATEddl`j''",quadcolsum(dmu`j'ddl)) ;
		};
		
		if "`margin`j''"=="nb" {;
			local _dLATEda`j' "\ 0" ; 
		};
		
		mata st_view(dmu`j'ds=.,.,tokens("`dmu`j'ds'"),"`touse'");
		mata st_matrix("`dLATEds`j''",quadcolsum(dmu`j'ds)) ;	
		local _dLATEds`j' "\ (2*`j'-1)*`dLATEds`j''" ;

		mata st_view(dmu`j'dr=.,.,tokens("`dmu`j'dr'"),"`touse'");
		mata st_matrix("`dLATEdr`j''",quadcolsum(dmu`j'dr)) ;	
		local _dLATEdr`j' "\ (2*`j'-1)*`dLATEdr`j''" ;
	};
	if missing("`cluster'") {;
		mata st_matrix("`S22'",quadcross(lte,lte)) ;
	};
	else {;
		sort `cluster';
		mata V = _cluster_var("`cluster'", "`lte'", "`touse'");
		mata st_matrix("`S22'",V); 
	};
	matrix `G21' = (-`dLATEdb0' \ `dLATEdb1' \ ((`dLATEdgu1'+`dLATEdgl1') - (`dLATEdgu0'+`dLATEdgl0')) `_dLATEdas' `_dLATEdds' `_dLATEda0' `_dLATEda1' `_dLATEds0' `_dLATEds1' `_dLATEdr0' `_dLATEdr1')/`N'; 
	matrix `Vlate' = (`S22'/`N' + `G21''*(`N'*`V')*`G21')/`N' ;
	matrix `Vdlate' = (`G21''*(`N'*`V')*`G21')/`N' ;
	
	return scalar Vdlate = `Vdlate'[1,1] ;
	return scalar Vlate = `Vlate'[1,1] ;
	return scalar late = `LATE' ;
	
	qui replace `iv' =  `z_temp';
	
	mata: mata drop Xsl Xsu X0 X1 lte dmu0db dmu1db dmu0dgl dmu0dgu dmu1dgl dmu1dgu ; 
	cap mata: mata drop dmu0ds dmu1ds ;
	cap mata: mata drop dmu0dr dmu1dr ;
	cap mata: mata drop dmu0da dmu1da ;
end;

program define mte_escount, rclass sortpreserve ;
	syntax  [if] [in] , at(string)
		;	//
	
	if "`e(cmd)'" != "escount" error 301;
		
	tempname beta betas V N MTE Vmte Vdmte S22 G21;
	tempvar zg ;
	tempvar	touse mte cons;
		
	qui gen `touse' = e(sample) ;
	matrix `beta' = e(b);	matrix `V' = e(V);
		
	local y0: word 1 of `e(depvar)';
	local y1: word 2 of `e(depvar)';
	local s: word 3 of `e(depvar)';
	local outoff `e(offset1)';
	local seloff `e(offset2)';
	local cluster `e(clustvar)';
	
	local _cons _cons;
	qui gen `cons' = _cons ;
	
	matrix `betas' = `beta'[1,"selection:"];	// coef vector of selection equatio
	matrix score double `zg' = `betas' if `touse' ;
	if ~missing("`seloff'") qui replace	`zg' = `zg' + `seloff' if `touse' ;
	
	local xs_temp : colnames `betas' ; 
	local xs: list  xs_temp - _cons;
	if "`xs_temp'"!="`xs'" {; local xs `xs' `cons';};
	mata st_view(Xs=.,.,"`xs'","`touse'");
	

	qui sum `touse' if `touse', meanonly; 	scalar `N' = r(N); 
	
	tempvar nu ;
	qui gen double `nu' = `at';
	
	foreach j in 0 1 {;
		local margin`j' `e(model)';
		
		tempvar xb`j' mu`j' dmu`j'db dmu`j'da dmu`j'ds dmu`j'dr dmu`j'da dmu`j'dd;
		tempname beta`j' alpha`j' lnalpha`j' athrho`j' rho`j' lnsigma`j' sigma`j' dMTEdb`j' dMTEds`j' dMTEda`j' dMTEdr`j' dMTEda`j' dMTEdd`j' ;
		
		
		matrix `beta`j'' = `beta'[1,"status`j':"];
		qui matrix score double `xb`j'' = `beta`j'' if `touse' ;
		if ~missing("`outoff'") qui replace `xb`j'' = `xb`j''+`outoff' if `touse' ;
		
		local x`j'_temp : colnames `beta`j'' ; local x`j': list  x`j'_temp - _cons;
		if "`x`j'_temp'"!="`x`j''" {; local x`j' `x`j'' `cons';};
		mata st_view(X`j'=.,.,"`x`j''","`touse'");
		
		scalar `lnsigma`j'' = _b[/lnsigma`j'];	scalar `sigma`j'' = exp(`lnsigma`j'');
		scalar `athrho`j'' = _b[/athrho`j'];	scalar `rho`j''= tanh(`athrho`j'');
		
		if "`margin`j''"=="nb" {;
			scalar `lnalpha`j'' = _b[/lnalpha`j'];
		};
		else {; scalar `lnalpha`j''= 0; };
		
		qui gen double `mu`j'' = exp(`xb`j'')*exp(`rho`j''*`sigma`j''*`nu'+0.5*`sigma`j''^2*(1-`rho`j''^2)) if `touse';
		qui gen double `dmu`j'db' = `mu`j'' if `touse';
		cap gen double `dmu`j'da' = `mu`j''*`rho`j''*`sigma`j''*`dnuda' if `touse';
		cap gen double `dmu`j'dd' = `mu`j''*`rho`j''*`sigma`j''*`dnudd' if `touse';
		cap gen double `dmu`j'ds' = `mu`j''*(`rho`j''*`nu'+`sigma`j''*(1-`rho`j''^2))*`sigma`j'' if `touse';
		qui gen double `dmu`j'dr' = `mu`j''*(`sigma`j''*`nu'-`sigma`j''^2*`rho`j'')*(sqrt(1-(`rho`j'')^2))^2 if `touse';
				
	};	// end of foreach j in 0 1 ;
	
	qui gen double `mte' = `mu1' - `mu0'  if  `touse' ;
	******************************************************************************************************;
	qui sum `mte' if `touse', meanonly ;	scalar `MTE' = r(mean);		qui replace `mte' = `mte' - `MTE';		
	
	mata st_view(mte=.,.,"`mte'","`touse'");
	
	*** average treatment effect *** ;
	foreach j in 0 1 {;
		mata st_view(dmu`j'db=.,.,tokens("`dmu`j'db'"),"`touse'");
		mata st_matrix("`dMTEdb`j''",quadcross(X`j',dmu`j'db)) ;
		
		if "`margin`j''"=="nb" {;
			local _dMTEda`j' "\ 0" ; 
		};
	
		mata st_view(dmu`j'ds=.,.,tokens("`dmu`j'ds'"),"`touse'");
		mata st_matrix("`dMTEds`j''",quadcolsum(dmu`j'ds)) ;
		local _dMTEds`j' "\ (2*`j'-1)*`dMTEds`j''" ;
		
		mata st_view(dmu`j'dr=.,.,tokens("`dmu`j'dr'"),"`touse'");
		mata st_matrix("`dMTEdr`j''",quadcolsum(dmu`j'dr)) ;
		local _dMTEdr`j' "\ (2*`j'-1)*`dMTEdr`j''" ;
	};
	
	if missing("`cluster'") {;
		mata st_matrix("`S22'",quadcross(mte,mte)) ;
	};
	else {;
		sort `cluster';
		mata V = _cluster_var("`cluster'", "`mte'", "`touse'");
		mata st_matrix("`S22'",V); 
	};
	matrix `G21' = (-`dMTEdb0' \ `dMTEdb1' \ J(colsof(`betas'),1,0)  `_dMTEdas' `_dMTEdds' `_dMTEda0' `_dMTEda1' `_dMTEds0' `_dMTEds1' `_dMTEdr0' `_dMTEdr1')/`N'; 
	matrix `Vmte' = (`S22'/`N' + `G21''*(`N'*`V')*`G21')/`N' ;
	matrix `Vdmte' = (`G21''*(`N'*`V')*`G21')/`N' ;
	
	return scalar Vdmte = `Vdmte'[1,1] ;
	return scalar Vmte = `Vmte'[1,1] ;
	return scalar mte = `MTE' ;
	
	mata: mata drop Xs X0 X1 mte dmu0db dmu1db ;
	cap mata: mata drop dmu0dr dmu1dr ;
	cap mata: mata drop dmu0da dmu1da ;
	cap mata: mata drop dmu0ds dmu1ds ;
end;

mata;
real matrix _cluster_var(string scalar _clust, string scalar _x, string scalar touse)
{
	st_view(x=.,.,tokens(_x),touse)
	st_view(clustvar = .,.,tokens(_clust),touse)
	clustinfo = panelsetup(clustvar,1)
	Nc = panelstats(clustinfo)[1]	// number of cluster
	K = cols(x)
	S = J(K,K,0)
	for (c=1;c<=Nc;++c) {
		xc = panelsubmatrix(x,c,clustinfo)
		S = S + quadcross(colsum(xc),colsum(xc))
		
	}
	S = (Nc/(Nc-1))*S
	return(S)
}
end;
