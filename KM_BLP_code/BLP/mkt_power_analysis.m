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
global ns x1 x2 vfull dfull theta1 cdid cdindex IV1 own nbrn alphai mvalold oldt2 s_jt

%**************************************************************************
%Define paths for input and output folders
%**************************************************************************
code_path                 =pwd;
optim_results_path        =[code_path,'\Optimization results\'];
mkt_power_results_path    =[code_path,'\Market power results\'];
merger_results_path       =[code_path,'\Merger results\'];
add_path                  =[code_path,'\optimization routines\'];
addpath(add_path);

%**************************************************************************
%Load data and define various parameters, including ownership dummies
%**************************************************************************
load BLP_data cdid cdindex share outshr price firmid id const hpwt air mpd space mpg trend product model_id own_dummies

load BLP_data_str model_name

s_jt=share;

% demand instruments
sum_other=[];
sum_rival=[];
X=[const,(hpwt),air,(mpg),space];
for i=1:size(id,1)
    other_ind=(firmid==firmid(i)  & cdid==cdid(i) & id~=id(i));
    rival_ind=(firmid~=firmid(i)  & cdid==cdid(i));
    sum_other(i,:)=sum(X(other_ind==1,:));
    sum_rival(i,:)=sum(X(rival_ind==1,:));
end
IV1=[X,sum_other,sum_rival];

%Load N(0,I) drwas
load v

%ownership structure matrix
own=own_dummies; clear owndummies

%variables in demand without random coefficients
x1=[price,X]; clear X

%variables in demand with random coefficients
x2=x1(:,1:size(x1,2)-1);

%#of indivduals and markets
ns =  size(v,2)/5;
nmkt = 20;

%#of brands per market
nbrn=zeros(nmkt,1);
nbrn(1)=cdindex(1);
for i=2:max(nmkt)
nbrn(i)=sum(cdid==i);
end

%demographics and N(0,I) unobservables
demogr=zeros(size(v));
vfull = v(cdid,:);
dfull = demogr(cdid,:);

%**************************************************************************
%Logit IV regression
%**************************************************************************
theta2w=     [3.612 0 0 0 0;
              4.628 0 0 0 0;
              1.818 0 0 0 0;
              1.050 0 0 0 0;
              2.056 0 0 0 0];
[theti, thetj, theta2]=find(theta2w);

invA1=inv(IV1'*IV1);
y = log(s_jt) - log(outshr);
mid = x1'*IV1*invA1*IV1';
t = inv(mid*x1)*mid*y;
oldt2 = zeros(size(theta2));
mvalold=x1*t;
mvalold_logit = mvalold;

mkt_power_results=[];

%**************************************************************************
%Loop over the various optimization routines
%**************************************************************************
for optrout=1:10

    %Optimization routine 6 (GA-JBES) did not produce reasonable results
    %in the optimization stage
    if optrout~=6
        
        cd(optim_results_path)

        if optrout<=9
            matfile=['blp_0',num2str(optrout),'_data_optim'];
            mkt_power_file=[mkt_power_results_path,'blp_mkt_power_results_0',num2str(optrout),'.txt'];            
            merger_file   =[merger_results_path,'blp_merger_results_0',num2str(optrout),'.txt'];            
        else
            matfile=['blp_',num2str(optrout),'_data_optim'];
            mkt_power_file=[mkt_power_results_path,'blp_mkt_power_results_',num2str(optrout),'.txt'];            
            merger_file   =[merger_results_path,'blp_merger_results_',num2str(optrout),'.txt'];                        
        end

        load (matfile, 'theta1s', 'theta2s','deltas','fvals','perturbs2');
        cd(code_path);

%        remove comments to perform the analysis only for the "best"
%        set of results
%        [min_fval,min_fval_ind]=min(fvals);
%        theta1s   = theta1s(min_fval_ind,:);
%        theta2s   = theta2s(min_fval_ind,:);
%        deltas    = deltas(:,min_fval_ind);
%        fvals     = fvals(min_fval_ind,:);
%        perturbs2 = perturbs2(min_fval_ind,:);

        mkt_power_results=[]; merger_results=[];

        for jj_optrout=1:1:size(fvals,1)

            theta1=theta1s(jj_optrout,:)';
            theta2=theta2s(jj_optrout,:)';
            delta=deltas(:,jj_optrout)';
            fval=fvals(jj_optrout,:);
            perturb=perturbs2(jj_optrout,:);

            mvalold=exp(mvalold_logit);
            oldt2 = zeros(size(theta2));
            deltajt=meanval(theta2);

            theta2w(:,1)=theta2;

            sijt=ind_sh(exp(deltajt),exp(mufunc(x2,theta2w)));
            sijt_pre=sijt;
            sjt_pre=(1/ns)*sum(sijt')';
            deltajt_pre=deltajt;

            vfull1=vfull(:,1:ns);
            alpha_i=[];
            for i=1:size(vfull1,1)
                alpha_i(i,:)=vfull1(i,:).*(kron(theta2(1),ones(1,ns)))+...
                    (kron(theta1(1),ones(1,ns)));
            end

            alphai=alpha_i;
            deriv_all=zeros(max(nbrn),max(nbrn),nmkt);
            elast_all=zeros(max(nbrn),max(nbrn),nmkt);

            for i=1:nmkt

                ind=cdid==i;
                pjt=price(ind==1,1);
                sjt=s_jt(ind==1,1);
                alpha_i=alphai(ind==1,:);
                s_ijt=sijt(ind==1,:);

                elast=zeros(size(pjt,1),size(pjt,1));
                deriv=zeros(size(pjt,1),size(pjt,1));

                for j=1:size(pjt,1)

                    deriv(j,j)=(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
                    elast(j,j)=(pjt(j)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));

                    for k=1:size(pjt,1)

                        if k~=j
                            elast(j,k)=-(pjt(k)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
                            deriv(j,k)=-(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));

                        end

                    end

                end

                elast_all(1:size(elast,1),1:size(elast,2),i)=elast;
                deriv_all(1:size(deriv,1),1:size(deriv,2),i)=deriv;

            end

            %store own and cross price elasticities
            elast_own=[];
            elast_cross=[];
            for i=1:nmkt
                temp=diag(elast_all(1:nbrn(i),1:nbrn(i),i));
                elast_own=[elast_own;temp];
                elast_cross=[elast_cross;elast_all(1:nbrn(i),:,i)];
            end

            %Consumer surplus calculations pre-merger
            exp_V=ind_eg(exp(deltajt),exp(mufunc(x2,theta2w)));
            tmp=[];
            CV_pre=[];
            for i=1:nmkt
                alphai_tmp=-alphai(cdid==i,:);
                alphai_tmp=alphai_tmp(1,:);
                tmp(i,:)=log(sum(exp_V(cdid==i,:))+1)./alphai_tmp;
                CV_pre(i,:)=tmp(i,:);
            end

            %Market_power calculation pre_merger
            own_dummy_pre=own;
            price_pre=price;

            mm=[];
            for i=1:max(cdid)
                p=price_pre(cdid==i,:);
                s=sjt_pre(cdid==i,:);
                nn=nbrn(i);
                om=deriv_all(1:nn,1:nn,i).*(own_dummy_pre(cdid==i,:)*own_dummy_pre(cdid==i,:)');
                m=-inv(om')*s;
                mm=[mm;m];
            end

            margin_pre=mm;
            mc=price_pre-margin_pre;
            margin_pct_pre=(margin_pre)./price_pre;
            profit_pre=margin_pre.*sjt_pre;

            optrout_aux=repmat(optrout,size(price,1),1);
            startval_aux=repmat(perturb,size(price,1),1);
            fval_aux=repmat(fval,size(price,1),1);
            market=cdid;

            mkt_power_results=[mkt_power_results;[optrout_aux,startval_aux,fval_aux,market,model_id,price_pre,sjt_pre,s_jt,elast_own,mc margin_pre,margin_pct_pre,elast_cross,mean(alphai')',std(alphai')']];

            %merger ownership matrix
            tmp=own(:,16)+own(:,19);
            own_dummy_post=[own(:,1:15),tmp,own(:,17:18),own(:,20:26)];

            mm=[];
            for i=1:max(cdid)
                p=price_pre(cdid==i,:);
                s=sjt_pre(cdid==i,:);
                nn=nbrn(i);
                om=deriv_all(1:nn,1:nn,i).*(own_dummy_post(cdid==i,:)*own_dummy_post(cdid==i,:)');
                m=-inv(om')*s;
                mm=[mm;m];
            end
            price_approx=mc+mm;
            price_post=price_approx;
            margin_post=mm;
            margin_pct_post=(margin_post)./price_post;

            %individual market shares post-merger
            deltajt_post=deltajt_pre-price_pre*theta1(1)+price_post*theta1(1);
            x2_post=x2;
            x2_post(:,1)=price_post;

            %calculate implied market shares
            theta2w = zeros(5,5);
            theta2w(:,1)=theta2;

            %update component of mu that corresponds to price
            [n k] = size(x2_post);
            j = size(theta2w,2)-1;
            mu = zeros(n,ns);
            for i = 1:ns
                v_i = vfull(:,i:ns:k*ns);
                d_i = dfull(:,i:ns:j*ns);
                mu(:,i) = (x2_post.*v_i*theta2w(:,1));
            end
            mu_post=mu;
            expmu=exp(mu);
            expmval=exp(deltajt_post);

            sijt_post=ind_sh(expmval,expmu);
            sjt_post=(1/ns)*sum(sijt_post')';

            %consumer surplus post-merger
            exp_V=ind_eg(expmval,expmu);
            tmp=[];
            CV_post=[];
            for i=1:nmkt
                alphai_tmp=-alphai(cdid==i,:);
                alphai_tmp=alphai_tmp(1,:);
                tmp(i,:)=log(sum(exp_V(cdid==i,:))+1)./alphai_tmp;
                CV_post(i,:)=tmp(i,:);
            end

            %profit post-merger
            profit_post=margin_post.*sjt_post;

            mean_CV=mean((CV_post-CV_pre)')';

            mean_CV_aux=[];
            for i=1:size(nbrn,1)
                tmp=nbrn(i);
                mean_CV_aux    =[mean_CV_aux;repmat(mean_CV(i,1),tmp,1)];
            end

            fprintf('optim routine : %2i\t',optrout);
            fprintf('start value   : %3i\t',perturb);
            fprintf('median elast  : %12.4f\n',median(diag(elast)));

            merger_results=[merger_results;[optrout_aux,startval_aux,market,model_id,price_pre,price_post,sjt_pre,sjt_post,mc,profit_pre,profit_post,mean_CV_aux]];
        end

        cd(mkt_power_results_path);     
        save(mkt_power_file,'mkt_power_results','-ASCII');
        
        cd(merger_results_path);     
        save(merger_file,'merger_results','-ASCII');
        
        cd(code_path);
                
    end    
end


