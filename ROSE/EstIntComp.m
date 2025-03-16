function lp = EstIntComp(intlist, ratio)
% function lp = EstIntComp(intlist, ratio)
% 估计亮度校准值，最小化std(MD)
% intlist: n*3 for 3 int list
% ratio: 0~1 for select ratio or n for max n points to calculate
%
if nargin < 2
    ratio = 1;
end

if ratio>1
    ratio = ratio/size(intlist,1);
end

if ratio>0 
    mask = rand(1,size(intlist,1)) < ratio;
    intlist = intlist(mask,:);
else
    warning('select ratio error, use all data.');
end
%% calculate intratio test
%intlist: temp
%minimize: std(md)
% intlist = temp(:,1:3);
disp(sprintf('Estimating Int Ratio, minimize std(MD)\npcnt = %d', size(intlist,1)));

x0 = [1 1];
options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.0001,'Algorithm','levenberg-marquardt');
lp=fminsearch(@(x0)CalPhaseScore(x0, intlist),x0,options);

%% Calculate md
tintlist1 = intlist;
tintlist1(:,2) = tintlist1(:,2) .*lp(1);
tintlist1(:,3) = tintlist1(:,3).*lp(2);
ret = CalPhaseListLite(tintlist1);

md1 = ret(:,4);
md1 = md1(md1>0.5 & md1<1.3);

disp(sprintf('Estimated int ratio: %0.3f, %0.3f\nmean MD: %0.3f\nstd MD: %0.3f', lp(1), lp(2), mean(md1), std(md1)));
