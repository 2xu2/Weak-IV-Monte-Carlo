% Jack Porter Summer Research group
% Monte carlo of Weak IV

%Set initial variables (Preset values)
n = 50;

pi0 = 0;
pi1m = [10, 1, 0.5, 0.1, 0.01];
b0 = 0;
b1 = 1;
r = 0;
sigma = 0;

%Set initial variables (joint distribution of z, w, episilon, v)
mvnmu = [0, 0, 0, 0];
mvnsigma = eye(4); 

%Montecarlo test
loops = 1000;
ro = [0, 0.25, 0.5, 0.9, 0.99]';
powermatrix = [];
for count1 = 1:size(ro, 1)
    %set correlation between episilon and v
    mvnsigma(4, 3) = ro(count1, 1);
    mvnsigma(3, 4) = ro(count1, 1); 
    for count2 = 1:size(pi1m, 2)
        pi1 = pi1m(1, count2);
        result = [];
        for i = 1:loops
            [stage2CI,stage1F, stage2t] = IVregress(n, b0, b1, r, pi0, pi1, sigma, mvnmu, mvnsigma); %IV regress
            [AR, tF] = ARtFtest(stage2t, stage1F); %AR and tF test
            b1inCI = (b1 > stage2CI(1, 1)) && (b1 < stage2CI(2, 1));
            result = [result; b1inCI, AR, tF];
        end
        poweroftest = sum(result(:,1))/loops;
        ARtestrate = sum(result(:,2))/loops;
        tFtestrate = sum(result(:,3))/loops;
        powermatrix(count1, count2) = poweroftest;
        ARmatrix(count1, count2) = ARtestrate;
        tFmatrix(count1, count2) = tFtestrate;
    end
end
resultofpower = [0, pi1m; ro, powermatrix]
resultofAR = [0, pi1m; ro, ARmatrix]
resultoftF = [0, pi1m; ro, tFmatrix]



%Still working 
%Monte carlo that fix everything but beta
% for i = 0:100
%     b0 =  0;
%     b1 = i/2 - 25;
%     [stage2CI,stage1F, stage2t] = IVregress(n, b0, b1, r, pi0, pi1, sigma, mvnmu, mvnsigma);
%     if (b1 > stage2CI(1, 1)) && (b1 < stage2CI(2, 1))
%         power = 
%     end
% end
% [AR, tF] = ARtFtest(stage2t, stage1F) % AR and tF Test

