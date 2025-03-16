function [lp, resnormx, resultp] = FitPhasePlane(x,y,p2, lp0, fast)
    if nargin<5
        fast = 0;
    end
    
    if fast==0
        options = optimset('Display','off','MaxFunEvals',1000,'MaxIter',100,'TolFun',1e-5,'LargeScale','on');
    else
        options = optimset('Display','off','MaxFunEvals',100,'MaxIter',5,'TolFun',1e-4,'LargeScale','off');
    end
	
    [lpx,resnormx,~,~]=lsqcurvefit(@(xp, xdata)CalPhaseSin(xp, x, y), ...
    [lp0(1) lp0(2) 0],x,[sin(p2) cos(p2)],[],[],options);

    lp = [lpx(1), lpx(2), checkPhase(lpx(3))];
    if nargout >=3
        resultp = CalPhase(lp, x, y); 
    end
end