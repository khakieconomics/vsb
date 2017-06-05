function [f]=merger_focs(p)
global cdid beta ns mc deltajt0 market_no nmkt...
    nbrn x2 own_dummy rcoeffs rcoeffs2

%new price
x22=x2(cdid==market_no,:);
x22(:,2)=p;

sijt_numer=NaN*zeros(1*nbrn,ns);
sijt_denom=NaN*zeros(1*nbrn,ns);

siot_numer=NaN*zeros(1*nbrn,ns);
siot_denom=NaN*zeros(1*nbrn,ns);

sijt=NaN*zeros(1*nbrn,ns);
siot=NaN*zeros(1*nbrn,ns);
deltajt=NaN*zeros(1*nbrn,ns);

for i=1:size(sijt_numer,1)
    k=1;
    for j=1:size(sijt_numer,2)
        deltajt(i)=deltajt0(i)+x22(i,2)*beta(2);        
        sijt_numer(i,j)=exp(deltajt(i)+x22(i,:)*rcoeffs2(i,k:k+3)');
        siot_numer(i,j)=exp(0);
        k=k+4;
    end
end

for i=1:size(sijt_denom,1)
    for j=1:size(sijt_denom,2)
        sijt_denom(i,j)=1+sum(sijt_numer(cdid==cdid(i),j));            
        siot_denom(i,j)=1+sum(sijt_numer(cdid==cdid(i),j));            
    end
end
sijt=sijt_numer./sijt_denom;   
siot=siot_numer./siot_denom;   

sjt=(1/ns)*sum(sijt')';
sot=(1/ns)*sum(siot')';

rcoeffs_price=rcoeffs(cdid==1,2:4:80);
alphai=rcoeffs_price;

l=1;
pjt=p;
deriv=zeros(size(pjt,1),size(pjt,1));
for j=1:size(pjt,1)
    deriv(j,j)=(1/ns)*sum(alphai(j,:).*sijt(j,:).*(ones(1,ns)-sijt(j,:)));
    for k=1:size(pjt,1)
        if k~=j
            deriv(j,k)=-(1/ns)*sum(alphai(j,:).*sijt(j,:).*(sijt(k,:)));
            l=l+1;
        end
    end
end
omega=(deriv).*(own_dummy*own_dummy');
f=sjt+omega'*(p-mc(cdid==market_no));
% fprintf('%15.8f\n',sum(sjt)+max(sot))
return
