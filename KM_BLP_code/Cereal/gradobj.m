function df = gradobj(theta2)

global invA IV gmmresid mvalold

temp = jacob(mvalold,theta2)';

df = 2*temp*IV*invA*IV'*gmmresid;

