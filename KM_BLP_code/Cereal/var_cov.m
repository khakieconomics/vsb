function f = var_cov(theta2)

global invA IV x1 mvalold gmmresid 

N = size(x1,1);
Z = size(IV,2);
temp = jacob(mvalold,theta2);
a = [x1 temp]'*IV;
IVres = IV.*(gmmresid*ones(1,Z));
b = IVres'*IVres;
f = inv(a*invA*a')*a*invA*b*invA*a'*inv(a*invA*a');

