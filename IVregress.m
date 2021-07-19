function [stage2CI,stage1F, stage2t] = IVregress(n, b0, b1, r, pi0, pi1, sigma, mvnmu, mvnsigma)
% IV regress function with data generation 
% Return stage1F and stage2t

% Data generation
F = mvnrnd(mvnmu, mvnsigma, n); %Joint distribution of z, w, episilon, v)

z = F(:,1);
w = F(:,2);
episilon = F(:,3);
v = F(:,4);
x = pi0 + z * pi1 + w * sigma + v;
y = b0 + x * b1 + w * r + episilon;
Z = [ones(n, 1), z, w];

%IV regression
%first stage regression
stage1b = inv(Z'* Z)* Z'* x; %stage 1 beta_hat
xhat = Z * stage1b;
stage1residual = x - xhat;
stage1sigma2 = (x - xhat)'*(x - xhat)/(n - 3);
temp = inv(Z'* Z);
stage1SE = sqrt(temp(2, 2) * stage1sigma2);
stage1t = abs(stage1b(2, 1))/stage1SE; %first stage t statistic of pi1
stage1F = ((stage1b(2,1))/stage1SE)^2; %first stage F statistic of pi1

%second stage regression
X = [ones(n, 1), x];
Xhat = [ones(n, 1), xhat];
stage2b = inv(Xhat'* Xhat)*Xhat'*y; %stage2 beta
yhat = X * stage2b;
stage2sigma2 = (y - yhat)'*(y - yhat)/(n - 3);
temp2 = inv(Xhat'* Xhat);
stage2SE = sqrt(temp2(2, 2) * stage2sigma2);
stage2t = abs(stage2b(2, 1))/stage2SE; %second stage t statistic of beta1
stage2CI = [0; 0];
stage2CI(1, 1) = stage2b(2, 1) - 1.96 * stage2SE;
stage2CI(2, 1) = stage2b(2, 1) + 1.96 * stage2SE; %Confidence interval of beta1

end

