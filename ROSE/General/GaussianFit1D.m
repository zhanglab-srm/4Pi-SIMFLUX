function [x0 y0 width ry] = GaussianFit1D(x,y,s0)
    if nargin < 3
        s0 = 3.0;
    end
    xp = paramest(x,y,s0);
    options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.01,'LargeScale','on');
    [lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian1D,xp,x,y,[],[],options);
    ry = Gaussian1D(lp, x);
    x0 = lp(1);
    width = lp(2);
    y0 = lp(3);
end

%param: mean,std,peak,bkg
function param = paramest(x,y,s0)
    param = [0 0 0 0];
    maxy = max(y);
    miny = min(y);
    t=find(y == maxy);
    maxx = x(t(1));
    
    param(1) = maxx;
    param(2) = s0;
    param(3) = maxy - miny;
    param(4) = miny;
end

function result = Gaussian1D(param, x)
    m = param(1);
    s = param(2);
    p = param(3);
    b = param(4);
    
    result = b + p * exp(-((x - m)/s).^2/2);
end