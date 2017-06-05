**************************************************************************
*To accompany Knittel and Metaxoglou (2008)
* Estimation of Random Coefficient Demand Models: 
* Challenges, Difficulties and Warnings
* Knittel      : crknittel@ucdavis.edu
* Metaxoglou   : konstantinos.metaxoglou@bateswhite.com
***************************************************************************
clear 
set memo 2000m
set more off
capture log close

global csv_file "\\sd-fs\admin\For Kostas\Nevo\~Nevo delivery\optimization results\Optimization results.csv"
global out_file "\\sd-fs\admin\For Kostas\Nevo\~Nevo delivery\optimization results\Optimization results.dta"
#delimit ;

insheet 

optmethod
stvalue       
fcnevals       
exitinfo       
toc
fval           

price_mean   
brand1_mean    
brand2_mean    
brand3_mean    
brand4_mean    
brand5_mean    
brand6_mean    
brand7_mean    
brand8_mean    
brand9_mean    
brand10_mean   
brand11_mean   
brand12_mean   
brand13_mean   
brand14_mean   
brand15_mean   
brand16_mean   
brand17_mean   
brand18_mean   
brand19_mean   
brand20_mean   
brand21_mean   
brand22_mean   
brand23_mean   
brand24_mean   

const_sigma    
price_sigma    
sugar_sigma    
mushy_sigma    
const_inc_sigma                
price_inc_sigma                
sugar_inc_sigma                
mushy_inc_sigma                
price_inc2_sigma
const_age_sigma
sugar_age_sigma
mushy_age_sigma                
price_child_sigma

price_mean_se   
brand1_mean_se    
brand2_mean_se    
brand3_mean_se    
brand4_mean_se    
brand5_mean_se    
brand6_mean_se    
brand7_mean_se    
brand8_mean_se    
brand9_mean_se    
brand10_mean_se   
brand11_mean_se   
brand12_mean_se   
brand13_mean_se   
brand14_mean_se   
brand15_mean_se   
brand16_mean_se   
brand17_mean_se   
brand18_mean_se   
brand19_mean_se   
brand20_mean_se   
brand21_mean_se   
brand22_mean_se   
brand23_mean_se   
brand24_mean_se   

const_sigma_se    
price_sigma_se    
sugar_sigma_se    
mushy_sigma_se    
const_inc_sigma_se                
price_inc_sigma_se                
sugar_inc_sigma_se                
mushy_inc_sigma_se                
price_inc2_sigma_se
const_age_sigma_se
sugar_age_sigma_se
mushy_age_sigma_se                
price_child_sigma_se

gradients1_const_sigma                
gradients1_price_sigma                
gradients1_sugar_sigma                
gradients1_mushy_sigma                
gradients1_const_inc                
gradients1_price_inc                
gradients1_sugar_inc                
gradients1_mushy_inc                
gradients1_price_inc2                
gradients1_const_age               
gradients1_sugar_age
gradients1_mushy_age
gradients1_price_child
gradients1_norm_inf
                
gradients2_const_sigma                
gradients2_price_sigma                
gradients2_sugar_sigma                
gradients2_mushy_sigma                
gradients2_const_inc                
gradients2_price_inc                
gradients2_sugar_inc                
gradients2_mushy_inc                
gradients2_price_inc2                
gradients2_const_age               
gradients2_sugar_age
gradients2_mushy_age
gradients2_price_child
gradients2_norm_inf

gradients3_const_sigma                
gradients3_price_sigma                
gradients3_sugar_sigma                
gradients3_mushy_sigma                
gradients3_const_inc                
gradients3_price_inc                
gradients3_sugar_inc                
gradients3_mushy_inc                
gradients3_price_inc2                
gradients3_const_age               
gradients3_sugar_age
gradients3_mushy_age
gradients3_price_child
gradients3_norm_inf

hessians_eig1  
hessians_eig2  
hessians_eig3  
hessians_eig4  
hessians_eig5  
hessians_eig6  
hessians_eig7  
hessians_eig8  
hessians_eig9  
hessians_eig10 
hessians_eig11 
hessians_eig12 
hessians_eig13 

hessians2_eig1 
hessians2_eig2 
hessians2_eig3 
hessians2_eig4 
hessians2_eig5 
hessians2_eig6 
hessians2_eig7 
hessians2_eig8 
hessians2_eig9 
hessians2_eig10
hessians2_eig11
hessians2_eig12
hessians2_eig13

using "$csv_file", clear;

local myvarlist

price_mean_se   
brand1_mean_se    
brand2_mean_se    
brand3_mean_se    
brand4_mean_se    
brand5_mean_se    
brand6_mean_se    
brand7_mean_se    
brand8_mean_se    
brand9_mean_se    
brand10_mean_se   
brand11_mean_se   
brand12_mean_se   
brand13_mean_se   
brand14_mean_se   
brand15_mean_se   
brand16_mean_se   
brand17_mean_se   
brand18_mean_se   
brand19_mean_se   
brand20_mean_se   
brand21_mean_se   
brand22_mean_se   
brand23_mean_se   
brand24_mean_se   

const_sigma_se    
price_sigma_se    
sugar_sigma_se    
mushy_sigma_se    
const_inc_sigma_se                
price_inc_sigma_se                
sugar_inc_sigma_se                
mushy_inc_sigma_se                
price_inc2_sigma_se
const_age_sigma_se
sugar_age_sigma_se
mushy_age_sigma_se                
price_child_sigma_se;


local myvarlist2
price_mean   
brand1_mean    
brand2_mean    
brand3_mean    
brand4_mean    
brand5_mean    
brand6_mean    
brand7_mean    
brand8_mean    
brand9_mean    
brand10_mean   
brand11_mean   
brand12_mean   
brand13_mean   
brand14_mean   
brand15_mean   
brand16_mean   
brand17_mean   
brand18_mean   
brand19_mean   
brand20_mean   
brand21_mean   
brand22_mean   
brand23_mean   
brand24_mean   

const_sigma    
price_sigma    
sugar_sigma    
mushy_sigma    
const_inc_sigma                
price_inc_sigma                
sugar_inc_sigma                
mushy_inc_sigma                
price_inc2_sigma
const_age_sigma
sugar_age_sigma
mushy_age_sigma                
price_child_sigma

price_mean_se   
brand1_mean_se    
brand2_mean_se    
brand3_mean_se    
brand4_mean_se    
brand5_mean_se    
brand6_mean_se    
brand7_mean_se    
brand8_mean_se    
brand9_mean_se    
brand10_mean_se   
brand11_mean_se   
brand12_mean_se   
brand13_mean_se   
brand14_mean_se   
brand15_mean_se   
brand16_mean_se   
brand17_mean_se   
brand18_mean_se   
brand19_mean_se   
brand20_mean_se   
brand21_mean_se   
brand22_mean_se   
brand23_mean_se   
brand24_mean_se   

const_sigma_se    
price_sigma_se    
sugar_sigma_se    
mushy_sigma_se    
const_inc_sigma_se                
price_inc_sigma_se                
sugar_inc_sigma_se                
mushy_inc_sigma_se                
price_inc2_sigma_se
const_age_sigma_se
sugar_age_sigma_se
mushy_age_sigma_se                
price_child_sigma_se

gradients1_const_sigma                
gradients1_price_sigma                
gradients1_sugar_sigma                
gradients1_mushy_sigma                
gradients1_const_inc                
gradients1_price_inc                
gradients1_sugar_inc                
gradients1_mushy_inc                
gradients1_price_inc2                
gradients1_const_age               
gradients1_sugar_age
gradients1_mushy_age
gradients1_price_child
gradients1_norm_inf
                
gradients2_const_sigma                
gradients2_price_sigma                
gradients2_sugar_sigma                
gradients2_mushy_sigma                
gradients2_const_inc                
gradients2_price_inc                
gradients2_sugar_inc                
gradients2_mushy_inc                
gradients2_price_inc2                
gradients2_const_age               
gradients2_sugar_age
gradients2_mushy_age
gradients2_price_child
gradients2_norm_inf

gradients3_const_sigma                
gradients3_price_sigma                
gradients3_sugar_sigma                
gradients3_mushy_sigma                
gradients3_const_inc                
gradients3_price_inc                
gradients3_sugar_inc                
gradients3_mushy_inc                
gradients3_price_inc2                
gradients3_const_age               
gradients3_sugar_age
gradients3_mushy_age
gradients3_price_child
gradients3_norm_inf

hessians_eig1  
hessians_eig2  
hessians_eig3  
hessians_eig4  
hessians_eig5  
hessians_eig6  
hessians_eig7  
hessians_eig8  
hessians_eig9  
hessians_eig10 
hessians_eig11 
hessians_eig12 
hessians_eig13 

hessians2_eig1 
hessians2_eig2 
hessians2_eig3 
hessians2_eig4 
hessians2_eig5 
hessians2_eig6 
hessians2_eig7 
hessians2_eig8 
hessians2_eig9 
hessians2_eig10
hessians2_eig11
hessians2_eig12
hessians2_eig13;


#delimit cr

*deal with irregular values
foreach v of local myvarlist {
	qui replace `v'="" if `v'=="NaN"	
	*imaginary numbers
	qui replace `v'="" if index(`v',"i")
	qui destring `v', replace
}

foreach v of local myvarlist2 {
	qui replace `v'=. if fval==10000000000 
}

sort optmethod stvalue
qui gen optmethod_str=""
qui replace optmethod_str="QNewton 1"        if optmethod==1
qui replace optmethod_str="Simplex"          if optmethod==2
qui replace optmethod_str="Solvopt"          if optmethod==3
qui replace optmethod_str="Conjgrad"         if optmethod==4
qui replace optmethod_str="QNewton 2"        if optmethod==5
qui replace optmethod_str="GA JBES"          if optmethod==6
qui replace optmethod_str="SA"  	     if optmethod==7
qui replace optmethod_str="MADS"             if optmethod==8
qui replace optmethod_str="GPS"              if optmethod==9
qui replace optmethod_str="GA Matlab"        if optmethod==10

format *mean* *se* *grad* *hess* %15.4fc

order optmethod_str optmethod stvalue fcnevals exitinfo toc fval

sort optmethod stvalue

compress

save "$out_file", replace

*eof