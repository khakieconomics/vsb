function f = jacob(mval,theta2)
global ns theti thetj cdid cdindex x2 v demogr

theta2w = full(sparse(theti,thetj,theta2));

expmu = exp(mufunc(x2,theta2w));
shares = ind_sh(mval,expmu);
clear expmu

[n,K] = size(x2);
J = size(theta2w,2) - 1;
f1 = zeros(size(cdid,1),K*(J + 1));

% computing (partial share)/(partial sigma)
for i = 1:K
	xv = (x2(:,i)*ones(1,ns)).*v(cdid,ns*(i-1)+1:ns*i);    
	temp = cumsum(xv.*shares);
	sum1 = temp(cdindex,:);
	sum1(2:size(sum1,1),:) = diff(sum1);
	f1(:,i) = mean((shares.*(xv-sum1(cdid,:)))')';
	clear xv temp sum1
end

% If no demogr comment out the next para
% computing (partial share)/(partial pi)
for j = 1:J
d = demogr(cdid,ns*(j-1)+1:ns*j);    
	temp1 = zeros(size(cdid,1),K);
	for i = 1:K
		xd=(x2(:,i)*ones(1,ns)).*d;    
		temp = cumsum(xd.*shares);
		sum1 = temp(cdindex,:);
		sum1(2:size(sum1,1),:) = diff(sum1);
		temp1(:,i) = mean((shares.*(xd-sum1(cdid,:)))')';
		clear xd temp sum1
	end
	f1(:,K*j+1:K*(j+1)) = temp1;
	clear temp1
end

rel = theti + (thetj - 1) * max(theti) ;

% computing (partial delta)/(partial theta2)
f = zeros(size(cdid,1),size(rel,1));
n = 1;
for i = 1:size(cdindex,1)
	temp = shares(n:cdindex(i),:);
	H1 = temp*temp';
	H = (diag(sum(temp')) - H1)/ns;
	f(n:cdindex(i),:) = - inv(H)*f1(n:cdindex(i),rel);
	n = cdindex(i) + 1;
end
