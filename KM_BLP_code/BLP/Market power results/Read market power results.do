**************************************************************************
*To accompany Knittel and Metaxoglou (2008)
* Estimation of Random Coefficient Demand Models: 
* Challenges, Difficulties and Warnings
*Knittel      : crknittel@ucdavis.edu
*Metaxoglou   : konstantinos.metaxoglou@bateswhite.com
***************************************************************************

clear
set type double
set more off
set memo 8000m

global inpath  "\\sd-fs\admin\For Kostas\BLP\~BLP delivery\Market power results\"
global outpath "\\sd-fs\admin\For Kostas\BLP\~BLP delivery\Market power results\"


use "$inpath\row_id.dta", clear
sort row_id
save "$inpath\row_id.dta", replace

use "$inpath\model_and_firm_id.dta", clear
capture gen obs=_n
sort obs
save "$inpath\model_and_firm_id.dta", replace

foreach i of numlist 1 2 3 4 5 7 8 9 10 {
	if `i'<10 {
        	capture infile optmethod stvalue fval market model_id price share_est share_obs elast_own mc margin margin_pct elast_cross1-elast_cross150 rcoeff_mean rcoeff_std using "${inpath}\blp_mkt_power_results_0`i'.txt", clear					
	}
	if `i'>=10 {
        	capture infile optmethod stvalue fval market model_id price share_est share_obs elast_own mc margin margin_pct elast_cross1-elast_cross150 rcoeff_mean rcoeff_std using "${inpath}\blp_mkt_power_results_`i'.txt", clear			
	}	
	drop model_id
	gen  row_id=_n
	sort row_id
	merge row_id using "$inpath\row_id.dta"
	tab _merge
	drop _merge

	sort  obs
	merge obs using "$inpath\model_and_firm_id.dta"
	tab _merge
	drop _merge

	qui {
	gen optmethod_str=""
	replace optmethod_str="Qnewton 1"        if `i'==1
	replace optmethod_str="Simplex"          if `i'==2
	replace optmethod_str="Solvopt"          if `i'==3
	replace optmethod_str="Conjgrad"         if `i'==4
	replace optmethod_str="Qnewton 2"        if `i'==5
	replace optmethod_str="GA JBES"          if `i'==6
	replace optmethod_str="SA"               if `i'==7	
	replace optmethod_str="MADS"             if `i'==8		
	replace optmethod_str="GPS"              if `i'==9		
	replace optmethod_str="GA Matlab"        if `i'==10		
	}

	format price share* elast_own mc margin margin_pct elast_cross* rcoeff* %12.5fc
	
	qui compress
	order optmethod optmethod_str stvalue fval market product price share_est share_obs elast_own mc margin margin_pct rcoeff_mean rcoeff_std elast_cross1-elast_cross150 
	
	if `i'<=9 {
		save "$inpath\blp_mkt_power_results0`i'.dta", replace	
	}
	if `i'>=9 {
		save "$inpath\blp_mkt_power_results`i'.dta", replace	
	}
}

use                  "$inpath\blp_mkt_power_results01.dta"
capture append using "$inpath\blp_mkt_power_results02.dta"
capture append using "$inpath\blp_mkt_power_results03.dta"
capture append using "$inpath\blp_mkt_power_results04.dta"
capture append using "$inpath\blp_mkt_power_results05.dta"
capture append using "$inpath\blp_mkt_power_results07.dta"
capture append using "$inpath\blp_mkt_power_results08.dta"
capture append using "$inpath\blp_mkt_power_results09.dta"
capture append using "$inpath\blp_mkt_power_results10.dta"
compress
save                "$inpath\blp_mkt_power_results_all.dta", replace

*eof





