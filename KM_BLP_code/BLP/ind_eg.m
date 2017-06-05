function f = ind_eg(expmval,expmu)

global ns 

eg = expmu.*kron(ones(1,ns),expmval);

f = eg;