%**************************************************************************
%To accompany Knittel and Metaxoglou (2008)
% Estimation of Random Coefficient Demand Models: 
% Challenges, Difficulties and Warnings
%Knittel      : crknittel@ucdavis.edu
%Metaxoglou   : konstantinos.metaxoglou@bateswhite.com
%**************************************************************************

clear all
close all
close hidden
warning off all
clc

%**************************************************************************
%Define globals
%**************************************************************************
global invA ns x1 x2 s_jt IV vfull dfull...
    theta1 theti thetj cdid cdindex optrout...
    omega mvalold oldt2 gmmresid v demogr

%**************************************************************************
%Define paths for input and output folders
%**************************************************************************
code_path                 =pwd;
optim_results_path        =[code_path,'\Optimization results\'];
mkt_power_results_path    =[code_path,'\Market power results\'];
add_path                  =[code_path,'\optimization routines\'];
addpath(add_path);

%**************************************************************************
%Load data and define various parameters, including ownership dummies
%**************************************************************************
load ps2
load iv
IV = [iv(:,2:21) x1(:,2:25)];
clear iv

merger=importdata('own_matrix.xls');
own= merger.data.pre_merger;

ns   = 20;    
nmkt = 94;  
nbrn = 24; 

cdid = kron([1:nmkt]',ones(nbrn,1));    
cdindex = [nbrn:nbrn:nbrn*nmkt]';       

%**************************************************************************
%Logit IV regression
%**************************************************************************
theta2w= [0.377,  3.089,       0,    1.186,         0;
          1.848, 16.598,  -0.659,        0,    11.625;
          0.004, -0.193,       0,    0.029,         0;
          0.081, 1.468,        0,   -1.514,         0];

[theti, thetj, theta2]=find(theta2w);
invA = inv([IV'*IV]);

temp = cumsum(s_jt);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
outshr = 1.0 - sum1(cdid,:);

y = log(s_jt) - log(outshr);
mid = x1'*IV*invA*IV';
t = inv(mid*x1)*mid*y;
mvalold = x1*t;
oldt2 = zeros(size(theta2));
mvalold_logit=mvalold;
mvalold = exp(mvalold);

clear mid y outshr t oldt2 mvalold temp sum1

dfull=demogr(cdid,:);
vfull=v(cdid,:);

%**************************************************************************
%Loop over the various optimization routines
%**************************************************************************
for optrout=1:10

    %Optimization routine 6 (GA-JBES) did not produce reasonable results
    %in the optimization stage
    if optrout~=6
        
        cd(optim_results_path)
        
        if optrout<=9
            matfile=['nevo_0',num2str(optrout),'_data_optim.mat'];
            mkt_power_file=[mkt_power_results_path,'nevo_mkt_power_results_0',num2str(optrout),'.txt'];
        else
            matfile=['nevo_',num2str(optrout),'_data_optim.mat'];
            mkt_power_file=[mkt_power_results_path,'nevo_mkt_power_results_',num2str(optrout),'.txt'];
        end
        
        load (matfile, 'perturbs2','fvals', 'theta1s', 'theta2s','exit_infos',...
                       'hessians','hessians2','gradients', 'gradients2',...
                       'gradients3','deltas' ,'gmmresids' ,'mvalolds2',...
                       'std_errors','counts2','fvals_track','tocs'); 

%        remove comments to perform the analysis only for the "best"
%        set of results
%        [min_fval,min_fval_ind]=min(fvals);
%        theta1s   = theta1s(min_fval_ind,:);
%        theta2s   = theta2s(min_fval_ind,:);
%        deltas    = deltas(:,min_fval_ind);
%        fvals     = fvals(min_fval_ind,:);
%        perturbs2 = perturbs2(min_fval_ind,:);
                                   
        cd(code_path);
        mkt_power_results=[];

        %******************************************************************
        %Loop over the various starting values
        %******************************************************************
        for jj_optrout=1:size(fvals,1);

            price=full(x1(:,1));

            theta1=theta1s(jj_optrout,:)';
            theta2=theta2s(jj_optrout,:)';
            fval=fvals(jj_optrout,:);
            perturb=perturbs2(jj_optrout,:);

            mvalold=exp(mvalold_logit);
            oldt2 = zeros(size(theta2));

            fval=gmmobj(theta2);
            deltajt=log(mvalold);
            gmmresid = deltajt - x1*theta1;

            vcov = full(var_cov(theta2));
            se   = sqrt(diag(vcov));

            theta2w = full(sparse(theti,thetj,theta2));
            t = size(se,1) - size(theta2,1);
            se2w = full(sparse(theti,thetj,se(t+1:size(se,1))));

            omega = inv(vcov(2:25,2:25));
            xmd = [x2(1:24,1) x2(1:24,3:4)];
            ymd = theta1(2:25);

            beta = inv(xmd'*omega*xmd)*xmd'*omega*ymd;
            resmd = ymd - xmd*beta;
            semd = sqrt(diag(inv(xmd'*omega*xmd)));
            mcoef = [beta(1); theta1(1); beta(2:3)];
            semcoef = [semd(1); se(1); semd];

            %the coefficients are as follows: constant price sugar mushy
            coeffs=[mcoef theta2w];
            beta=coeffs(:,1);
            sigma=coeffs(:,2);
            pai=coeffs(:,3:6);

            %construct the various random coefficients
            rcoeffs=NaN*zeros(size(cdid,1),ns);
            rcoeffs2=NaN*zeros(size(cdid,1),ns);

            for i=1:size(rcoeffs,1)
                k=1;
                for j=1:ns
                    tmpd=dfull(i,j:ns:80);
                    tmpv=vfull(i,j:ns:80);
                    rcoeff_tmp=pai*tmpd'+diag(sigma)*tmpv';
                    rcoeffs(i,k:k+3)=(beta+rcoeff_tmp)';
                    rcoeffs2(i,k:k+3)=(rcoeff_tmp)';
                    k=k+4;
                end
            end
            rcoeffs_const=rcoeffs(:,1:4:80);
            rcoeffs_price=rcoeffs(:,2:4:80);
            rcoeffs_sugar=rcoeffs(:,3:4:80);
            rcoeffs_mushy=rcoeffs(:,4:4:80);

            %individual market shares' components
            sijt_numer=NaN*zeros(nmkt*nbrn,ns);
            sijt_denom=NaN*zeros(nmkt*nbrn,ns);
            siot_numer=NaN*zeros(nmkt*nbrn,ns);
            siot_denom=NaN*zeros(nmkt*nbrn,ns);

            sijt=NaN*zeros(nmkt*nbrn,ns);
            siot=NaN*zeros(nmkt*nbrn,ns);

            for i=1:size(sijt_numer,1)
                k=1;
                for j=1:size(sijt_numer,2)
                    sijt_numer(i,j)=exp(deltajt(i)+x2(i,:)*rcoeffs2(i,k:k+3)');
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
            siot=sijt_numer./siot_denom;

            alphai=rcoeffs_price;
            alphaib=alphai;
            sjt=(1/ns)*sum(sijt')';
            sot=(1/ns)*sum(siot')';
            
            sijt_pre=sijt;

            %product market shares and mean utility levels            
            sot_pre=sot;
            sjt_pre=sjt;
            sot_pte=sot;
            deltajt_pre=deltajt;
            
            %derive matrices of price derivatives and elasticities
            deriv_all=zeros(max(nbrn),max(nbrn),nmkt);
            elast_all=zeros(max(nbrn),max(nbrn),nmkt);
            for i=1:max(cdid)
                l=1;
                ind=cdid==i;
                pjt=x1(ind==1,1);
                sjt=s_jt(ind==1,1);
                alpha_i=alphai(ind==1,:);
                s_ijt=sijt(ind==1,:);
                elast=zeros(size(pjt,1),size(pjt,1));
                deriv=zeros(size(pjt,1),size(pjt,1));
                for j=1:size(pjt,1)
                    deriv(j,j)=(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
                    elast(j,j)=(pjt(j)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
                    deriv_own(i,j)=deriv(j,j);
                    elast_own(i,j)=elast(j,j);
                    for k=1:size(pjt,1)
                        if k~=j
                            elast(j,k)=-(pjt(k)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
                            deriv(j,k)=-(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
                            elast_cross(i,l)=elast(j,k);
                            deriv_cross(i,l)=deriv(j,k);
                            l=l+1;
                        end
                    end
                end
                elast_all(1:size(elast,1),1:size(elast,2),i)=elast;
                deriv_all(1:size(deriv,1),1:size(deriv,2),i)=deriv;
            end
            
            temp=[];
            temp2=[];
            for j=1:nmkt
                temp=[temp; (elast_all(:,:,j))];
                temp2=[temp2; diag(elast_all(:,:,j))];
            end
            elast_all=temp;
            elast_own=temp2;

            %Market power calculations
            own_dummy_pre=own;
            price_pre=price;

            mm=[];
            for i=1:max(cdid)
                p=price_pre(cdid==i,:);
                s=sjt_pre(cdid==i,:);
                om=deriv_all(:,:,i).*(own_dummy_pre*own_dummy_pre');
                m=-inv(om')*s;
                mm=[mm;m];
            end
            
            margin_pre=mm;
            mc=price_pre-margin_pre;
            margin_pct_pre=(margin_pre)./price_pre;

            optrout_aux=repmat(optrout,size(price,1),1);
            startval_aux=repmat(jj_optrout,size(price,1),1);
            market= kron((1:1:94)',ones(nbrn,1));
            brand=repmat((1:1:24)',94,1);
            fval_aux=repmat(fval,size(price,1),1);

            mkt_power_results=[mkt_power_results;[optrout_aux,startval_aux,...
                               fval_aux,market,brand,price_pre,sjt_pre,s_jt,elast_own,...
                               mc margin_pre,margin_pct_pre,elast_all,mean(alphai')',std(alphai')']];
            mm=[];
            fprintf('mkt power routine:%3i\t',optrout);
            fprintf('starting value:%3i\t',perturb);
            fprintf('median elast_own:%12.4f\n',median(elast_own));
        end

    end

    cd(mkt_power_results_path);     
    save(mkt_power_file,'mkt_power_results','-ASCII');
    cd(code_path);
     
end