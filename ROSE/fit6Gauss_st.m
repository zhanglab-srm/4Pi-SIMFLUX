%% fit six images
%data in timgbuf in saze n*n*6*imglen
%parameter: [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
%           [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]
function result = fit6Gauss_st(timgbuf_)
    
    
    
    s = size(timgbuf_);
    [xx, yy] = meshgrid((1:s(2)), (1:s(1)));
    xdata=[xx,yy];
    xoffset = (s(2)+1)/2;
    yoffset = (s(1)+1)/2;
    fitlen = size(timgbuf_, 4);
    result = zeros(fitlen, 17);
    initparList = fitGaus_st(reshape(sum(timgbuf_,3),[s(1),s(2),fitlen]), 1);
    for m=1:fitlen
        timgbuf = timgbuf_(:,:,:,m);
        initpar = initparList(m,:);
        Fd=double(timgbuf);
        bkglist=min(min(Fd,[],1), [],2);
        bkglist = bkglist(:);
        peaklist=max(max(Fd,[],1), [],2);
        peaklist = peaklist(:);
        initstdx = initpar(3);
        initstdy = initpar(4);
        cx=initpar(1) + xoffset;
        cy=initpar(2) + yoffset;
        xp=[cx,cy,initstdx,initstdy,peaklist',bkglist'];
        options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.0001,'Algorithm','levenberg-marquardt');
        [lp,~,~,exitflag]=lsqcurvefit(@Gaussian_6,xp,xdata,Fd,[],[],options);

        result(m,:) = [lp exitflag];
    end
    result(:,1) = result(:,1) - xoffset;
    result(:,2) = result(:,2) - yoffset;
end

function f=Gaussian_6(xp,xdata)
    [a, b]=size(xdata);
    Xd=xdata(:,1:b/2);
    Yd=xdata(:,b/2+1:b);
    f1=xp(5)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(11);
    f2=xp(6)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(12);
    f3=xp(7)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(13);
    f4=xp(8)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(14);
    f5=xp(9)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(15);
    f6=xp(10)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(16);
    
    f = cat(3, f1,f2,f3,f4,f5,f6);
end