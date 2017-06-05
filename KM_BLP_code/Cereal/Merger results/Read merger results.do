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
set memo 3000m

global inpath  "\\sd-fs\admin\For Kostas\Nevo\~Nevo delivery\Merger results"
global outpath "\\sd-fs\admin\For Kostas\Nevo\~Nevo delivery\Merger results"

foreach i of numlist 1 2 3 4 5 7 8 9 10{
	disp "nevo_merger_results0`i'"
	
	if `i'<=9 {
		qui infile optmethod startval market brand price_pre price_post share_pre share_post mc profit_pre profit_post mean_CV CV_pre1-CV_pre20 CV_post1-CV_post20 using "$inpath\nevo_merger_results_0`i'.txt", clear
	}
	if `i'>9 {
		qui infile optmethod startval market brand price_pre price_post share_pre share_post mc profit_pre profit_post mean_CV CV_pre1-CV_pre20 CV_post1-CV_post20 using "$inpath\nevo_merger_results_`i'.txt", clear
	}
	
	drop CV_pre* CV_post*
	sort optmethod market brand

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
	
	gen firm_pre=1
	qui replace firm_pre=2 if brand>=10 & brand<=18
	qui replace firm_pre=3 if brand>=19 & brand<=20
	qui replace firm_pre=4 if brand>=21 & brand<=23
	qui replace firm_pre=5 if brand>=24
	
	gen firm_post=1
	qui replace firm_post=3 if brand>=19 & brand<=20
	qui replace firm_post=4 if brand>=21 & brand<=23
	qui replace firm_post=5 if brand>=24
			
	qui compress
	
	format price_pre price_post share_pre share_post mc profit_pre profit_post mean_CV %15.6fc
	
	if `i'<=9 {
	save "$inpath\nevo_merger_results0`i'.dta", replace	
	}
	
	if `i'>9 {
	save "$inpath\nevo_merger_results`i'.dta", replace	
	}	
	
}	

use                  "${outpath}\nevo_merger_results01.dta", clear
capture append using "${outpath}\nevo_merger_results02.dta"
capture append using "${outpath}\nevo_merger_results03.dta"
capture append using "${outpath}\nevo_merger_results04.dta"
capture append using "${outpath}\nevo_merger_results05.dta"
capture append using "${outpath}\nevo_merger_results07.dta"
capture append using "${outpath}\nevo_merger_results08.dta"
capture append using "${outpath}\nevo_merger_results09.dta"
capture append using "${outpath}\nevo_merger_results10.dta"
compress
save "${outpath}\nevo_merger_results_all.dta", replace

*eof