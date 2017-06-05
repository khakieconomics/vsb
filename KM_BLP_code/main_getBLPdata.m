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
global ns x1 x2 s_jt vfull dfull theta1 theti thetj...
       cdid cdindex IV1 invA1 nbrn...
       mvalold_logit v mvalold oldt2 gmmresid...
       mymaxfunevals mvalolds mvalold0 mvalold00 mvalolds2...
       fvals_track fcnevals ppp mval_track       

load BLP_data cdid cdindex share outshr price firmid id const hpwt air mpd space mpg trend...
     product model_id own_dummies

load BLP_data_str model_name

% *************************************************************************
% Define paths for codes, optimization results and logs
% *************************************************************************
code_path    =pwd;
results_path =[code_path,'\optimization results\'];
logs_path    =[code_path,'\optimization logs\'];
add_path     =[code_path,'\optimization routines\'];
addpath(add_path);

%**************************************************************************
% Loop over optimization routines
%**************************************************************************

for optrout=6:6
    
    perturbs=   (1:1:50)';
    mytolx=         1e-3;
    mytolfun=       1e-3;
    mymaxiters=   5*10^5;
    mymaxfunevals=  4000;
    
    fvals_track=NaN*ones(mymaxfunevals,size(perturbs,1));    

    if optrout<=9
       outfile=[logs_path,['blp_0',num2str(optrout),'_optim_log.txt']];
       matfile=['blp_0',num2str(optrout),'_data_optim'];
    else
       outfile=[logs_path,['blp_',num2str(optrout),'_optim_log.txt']];
       matfile=['blp_',num2str(optrout),'_data_optim'];
    end

    fid = fopen(outfile,'w'); fclose(fid);

    % *************************************************************************
    counts2    =[];                  %store function evaluations
    deltas     =[];                  %store deltas
    exit_infos =[];                  %store exit info
    fvals      =[];                  %store GMM values
    gmmresids  =[];                  %store gmm residuals
    gradients  =[];                  %store analytical gradients
    gradients2 =[];                  %store numerical  gradients I
    gradients3 =[];                  %store numerical  gradients II   
    hessians   =[];                  %store hessians
    hessians2  =[];                  %store hessians II    
    mvalolds2  =[];                  %store mvalolds
    perturbs2  =[];                  %store perturbation set number
    std_errors =[];                  %store std.errors
    theta1s    =[];                  %store theta1s
    theta2s    =[];                  %store theta2s
    fvals_track=[];                  %store GMM values in all evaluations
    tocs       =[];                  %store time of completion

    % *********************************************************************
    % Demand instruments
    % *********************************************************************
    sum_other=[];
    sum_rival=[];

    X=[const,(hpwt),air,(mpg),space];
    for i=1:size(id,1)
        other_ind=(firmid==firmid(i)  & cdid==cdid(i) & id~=id(i));
        rival_ind=(firmid~=firmid(i)  & cdid==cdid(i));
        total_ind=(cdid==cdid(i));
        sum_other(i,:)=sum(X(other_ind==1,:));
        sum_rival(i,:)=sum(X(rival_ind==1,:));
        sum_total(i,:)=sum(X(total_ind==1,:));
    end
    IV1=[X,sum_other,sum_rival];

    %Load N(0,I) drwas
    load v

    % %ownership structure matrix
    % own=own_dummies; clear owndummies

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

    % *********************************************************************
    % Logit regressions
    % *********************************************************************
    theta2=ones(5,1);
    theta2w=zeros(5,5);
    theta2w(:,1)=theta2;
    [theti, thetj, theta2]=find(theta2w);

    invA1 = inv(IV1'*IV1);

    s_jt=share;
    y = log(s_jt) - log(outshr);

    mid = x1'*IV1*invA1*IV1';
    theta1 = inv(mid*x1)*mid*y;
    mvalold_logit = x1*theta1;

    n=size(x1,1);
    k=size(x1,2);

    ESS=y'*y-2*theta1'*x1'*y+theta1'*x1'*x1*theta1;
    s2=ESS/(n-k);
    A=(x1'*(IV1*inv(IV1'*IV1)*IV1')*x1);
    se=sqrt(diag(s2*inv(A)));

    mvalold=exp(mvalold_logit);
    oldt2=zeros(size(theta2));
	

	datadir = '/Users/susanavasserman/Dropbox/WorkSpaces/Discrete_Choice_Models/KM_blp_data';
	% datadir = '/Users/susanavasserman/Dropbox/WorkSpaces/Discrete_Choice_Models/KM_blp_data';

	cd(datadir) 

	load('BLP_data.mat')


	logshare = log(share);
	logoutshare = log(outshr);


	dlmwrite('cdid.txt',cdid,'delimiter','\t')
	dlmwrite('price.txt',price,'delimiter','\t')
	dlmwrite('x1.txt',x1,'delimiter','\t')
	dlmwrite('x2.txt',x2,'delimiter','\t')
	dlmwrite('IV1.txt',IV1,'delimiter','\t')
	dlmwrite('id.txt',id,'delimiter','\t')
	dlmwrite('firmid.txt',firmid,'delimiter','\t')
	dlmwrite('share.txt',share,'delimiter','\t')
	dlmwrite('outshr.txt',outshr,'delimiter','\t')
	dlmwrite('logshare.txt',logshare,'delimiter','\t')
	dlmwrite('logoutshare.txt',logoutshare,'delimiter','\t')
	

    end %perturbations loop

end %optimization routines loop
		
		
