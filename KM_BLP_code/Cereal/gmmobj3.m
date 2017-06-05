function [f,g] = gmmobj3(theta2,par)

global invA theta1 x1 IV fcnevals delta gmmresid ppp mymaxfunevals...
fcnevals fvals_track

if size(theta2,1)==1
    theta2=theta2';
end

delta = meanval(theta2); par=0;

if max(isnan(delta)) == 1
    f=1e+10;
    g=-999999*ones(size(theta2));
    gmmresid=1e+10*ones(size(delta));
    fcnevals=fcnevals+1;
    if fcnevals<=mymaxfunevals
        fvals_track(fcnevals,ppp)=f;
    end        
else
	temp1 = x1'*IV;
	temp2 = delta'*IV;
    theta1 = inv(temp1*invA*temp1')*temp1*invA*temp2';
	gmmresid = delta - x1*theta1;
	temp1 = gmmresid'*IV;
	f = temp1*invA*temp1';
    fcnevals=fcnevals+1;
    if fcnevals<=mymaxfunevals
        fvals_track(fcnevals,ppp)=f;
    end            
    g=gradobj(theta2);
end

