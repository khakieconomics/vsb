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
global invA ns x1 x2 v demogr s_jt IV vfull dfull...
    theta1 theti thetj cdid cdindex optrout...
    mymaxfunevals mvalold...
    oldt2 gmmresid mvalolds mvalold0 mvalold00 mvalolds2...
    fvals_track fcnevals ppp mval_track

% *************************************************************************
% Define paths for codes, optimization results and logs
% *************************************************************************
code_path    =pwd;
results_path =[code_path,'\optimization results\'];
logs_path    =[code_path,'\optimization logs\'];
add_path     =[code_path,'\optimization routines\'];
addpath(add_path);

% *************************************************************************
% Loop over optimization routines
% *************************************************************************
for optrout=6:6

    mytolx=1e-3;
    mytolfun=1e-3;
    mymaxiters=5*10^5;
    mymaxfunevals=400;
    perturbs=(1:1:1)';
    
    fvals_track=NaN*ones(mymaxfunevals,size(perturbs,1));

    if optrout==1,  perturbs=perturbs(perturbs~=6,:);
        perturbs=perturbs(perturbs~=39,:);
        perturbs=perturbs(perturbs~=22,:);
    end
    if optrout==2,  perturbs=perturbs(perturbs~=6,:);
        perturbs=perturbs(perturbs~=28,:);
        perturbs=perturbs(perturbs~=36,:);
    end

    if optrout<=9
        outfile=[logs_path,['nevo_0',num2str(optrout),'_optim_log.txt']];
        matfile=['nevo_0',num2str(optrout),'_data_optim.mat'];
    else
        outfile=[logs_path,['nevo_',num2str(optrout),'_optim_log.txt']];
        matfile=['nevo_',num2str(optrout),'_data_optim.mat'];
    end

    % *********************************************************************
    % Initialize log files and matrices containing various results
    % *********************************************************************
    fid        =fopen(outfile,'w'); fclose(fid);
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
    % Load data
    % *********************************************************************

    load ps2
    load iv
    IV = [iv(:,2:21) x1(:,2:25)];
    clear iv

    ns   = 20;
    nmkt = 94;
    nbrn = 24;

    cdid    = kron([1:nmkt]',ones(nbrn,1));
    cdindex = [nbrn:nbrn:nbrn*nmkt]';

    % *********************************************************************
    % Logit IV regression
    % *********************************************************************

    dfull=demogr(cdid,:);
    vfull=v(cdid,:);

    invA = inv([IV'*IV]);

    temp = cumsum(s_jt);
    sum1 = temp(cdindex,:);
    sum1(2:size(sum1,1),:) = diff(sum1);
    outshr = 1.0 - sum1(cdid,:);

    y = log(s_jt) - log(outshr);
    mid = x1'*IV*invA*IV';
    t = inv(mid*x1)*mid*y;
    mvalold_logit = x1*t;

    n=size(x1,1);
    k=size(x1,2);

    ESS=y'*y-2*t'*x1'*y+t'*x1'*x1*t;
    s2=ESS/(n-k);
    A=(x1'*(IV*inv(IV'*IV)*IV')*x1);
    se=sqrt(diag(s2*inv(A)));

    % *********************************************************************
    % Start optimization routine with perturb_no different starting values
    % for theta2: normrnd(0,1,size(theta2));
    % for delta:  delta_logit+normrnd(0,stddev(delta_logit),2256,1)
    % *********************************************************************
%     diary(outfile)

    %********* SAVE DATA ********* %
    datadir = '/Users/susanavasserman/Dropbox/WorkSpaces/Discrete_Choice_Models/KM_cereal_data';
% datadir = '/Users/susanavasserman/Dropbox/WorkSpaces/Discrete_Choice_Models/KM_blp_data';

cd(datadir) 

% load('BLP_data.mat')

price = x1(:,1);
logshare = log(s_jt);
logoutshare = log(outshr);


dlmwrite('cdid.txt',full(cdid),'delimiter','\t')
dlmwrite('price.txt',full(price),'delimiter','\t','precision', 16)
dlmwrite('x1.txt',full(x1),'delimiter','\t','precision', 16)
dlmwrite('x2.txt',full(x2),'delimiter','\t','precision', 16)
dlmwrite('IV1.txt',full(IV),'delimiter','\t','precision', 16)
dlmwrite('id.txt',full(id),'delimiter','\t')
% dlmwrite('firmid.txt',full(firmid),'delimiter','\t')
dlmwrite('share.txt',full(s_jt),'delimiter','\t','precision', 16)
dlmwrite('outshr.txt',full(outshr),'delimiter','\t','precision', 16)
dlmwrite('logshare.txt',full(logshare),'delimiter','\t','precision', 16)
dlmwrite('logoutshare.txt',full(logoutshare),'delimiter','\t','precision', 16)

    
end