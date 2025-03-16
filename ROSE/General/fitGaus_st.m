function [result] = fitGaus_st(img, initstd)
% [result] = fitGaus(img, initstd)
% result: [x y wx wy bkg int exitflag]
%          1 2  3  4  5   6    7

    if nargin <2
        initstd = 1.5;
    end
    s = size(img);
    imglen = size(img, 3);
    [xx, yy] = meshgrid((1:s(2)), (1:s(1)));
    xdata=[xx,yy];
    
    result = zeros(imglen, 6);
    xoffset = (s(2)+1)/2;
    yoffset = (s(1)+1)/2;
    for i=1:imglen
        Fd=double(img(:,:,i));
        backg=min(Fd(:));
        peak=max(Fd(:));
        [cy, cx]=find(Fd==peak);
        cy=mean(cy);
        cx=mean(cx);
        xp=[cx,cy,initstd,initstd,backg,peak];
        options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',1000,'TolFun',0.01,'Algorithm','levenberg-marquardt');
        [lp,~,~,exitflag]=lsqcurvefit(@Gaussian_PALM1,xp,xdata,Fd,[],[],options);
    
        result(i, 1:6) = lp;
        result(i, 7) = exitflag;
    end
    result(:, 1) = result(:, 1)-xoffset;
    result(:, 2) = result(:, 2)-yoffset;
    result(:, 6) = result(:, 6) .*2*pi.*result(:, 3).*result(:, 4);
    
end

function f=Gaussian_PALM1(xp,xdata)
    [a, b]=size(xdata);
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