function [result] = fitGaussian_Hex_Oval_dpos(img, r, phase_offset)
    s=size(img);
    if nargin < 3
        phase_offset = 0;
    end
    
    if nargin < 2
        r = 4;
    end
    
    dposmax = 0.2;
%     [xx yy] = meshgrid(1:s(2), 1:s(1));
    img = double(img);
    xdata = img;
    %parameter
    %[x0 y0 r std int1 int2 int3 int4 int5 int6 phase0 bkg]
    %[x0 y0 std int1 int2 int3 int4 int5 int6 bkg]
    xp=[0, 0, 1, 1000, 1000, 1000, 1000, 1000, 1000, 100, r, r, phase_offset, ...
        0,0,0,0,0,0, ...
        0,0,0,0,0,0 ];
%     xp = p0;
    options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',1000,'TolFun',1e-5,'LargeScale','on');
    [lp,resnorm,residual,exitflag]=lsqcurvefit(@(xp, xdata)Gaussian_Hex_Oval_dpos(xp,xdata),xp,xdata,img, ...
        [-2 -2 0.5 0 0 0 0 0 0 0 r-1.5 r-1.5 phase_offset-15 -dposmax -dposmax -dposmax  -dposmax -dposmax -dposmax -dposmax -dposmax -dposmax -dposmax -dposmax -dposmax], ...
        [2 2 1.5 500000 500000 500000 500000 500000 500000 100000 r+1.5 r+1.5 phase_offset+15 dposmax dposmax dposmax dposmax dposmax dposmax dposmax dposmax dposmax dposmax dposmax dposmax],options);
    
    result = lp;
end