function [result] = fitGaus(img, initstd)
% [result] = fitGaus(img, initstd)
% result:
% result.x
% result.y
% result.wx
% result.wy
% result.bkg
% result.int

    if nargin <2
        initstd = 1.5;
    end
    s = size(img);
    imglen = size(img, 3);
    [xx yy] = meshgrid((1:s(2)), (1:s(1)));
    xdata=[xx,yy];
    
    result_x = zeros(imglen, 1);
    result_y = zeros(imglen, 1);
    result_wx = zeros(imglen, 1);
    result_wy = zeros(imglen, 1);
    result_int = zeros(imglen, 1);
    result_bkg = zeros(imglen, 1);
    

    parfor i=1:imglen
        Fd=double(img(:,:,i));
        backg=min(Fd(:));
        peak=max(Fd(:));
        [cy cx]=find(Fd==peak);
        cy=mean(cy);
        cx=mean(cx);
        xp=[cx,cy,initstd,initstd,backg,peak];
        options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',1000,'TolFun',0.01,'Algorithm','levenberg-marquardt');
        [lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian_PALM1,xp,xdata,Fd,[],[],options);
    
        result_x(i) = lp(1);
        result_y(i) = lp(2);
        result_wx(i) = lp(3);
        result_wy(i) = lp(4);
        result_bkg(i) = lp(5);
        result_int(i) = lp(6); 

    end

    
    result.x = result_x - (s(2)+1)/2;
    result.y = result_y - (s(1)+1)/2;
    result.wx = result_wx;
    result.wy = result_wy;
    result.bkg = result_bkg;
    result.int = result_int;
end

function f=Gaussian_PALM1(xp,xdata)
    [a b]=size(xdata);
    Xd=xdata(:,1:b/2);
    Yd=xdata(:,b/2+1:b);
    f=xp(6)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(5);
end
% function f=Gaussian_PALM1(xp,xdata)
%     [a b]=size(xdata);
%     Xd=xdata(:,1:b/2);
%     Yd=xdata(:,b/2+1:b);
%     f=xp(6)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(5);
% end