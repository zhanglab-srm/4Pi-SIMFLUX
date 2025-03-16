function [fitresult, info1, info2] = GaussFit6_CPU(img, initstd, options, maxiter)
    if nargin < 2 || isempty(initstd)
        initstd = 1.0;
    end
    
    if nargin < 3 || isempty(options)
        options = [1e-3, 1e-15, 1e-15, 1e-15, 1e-15];
    end
    
    if nargin < 4 || isempty(maxiter)
        maxiter = 200;
    end
    
    s = size(img);
    
    [lp, info1] = GaussFit2D_mex(squeeze(sum(img,3)), [options, maxiter], initstd);
    [lp2, info2] = GaussFit6_mex(img, [options, maxiter], lp);
    %parameter: [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
    %           [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]
    fitresult = zeros(size(img,4), 17);
    fitresult(:,1) = lp(:,1) - (s(2)-1)/2;
    fitresult(:,2) = lp(:,2) - (s(1)-1)/2;
    fitresult(:,3:4) = lp(:,3:4);
    fitresult(:,5:16) = lp2(:,1:12);
    fitresult(:,17) = info2(:,7);
end