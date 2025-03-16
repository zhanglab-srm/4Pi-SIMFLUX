function [currzang_cor, lp] = phaseCorrection_v3(currx,curry,currzang,currzfresult,currI,currcrlb,currt,num_images,lp0)

mask = currt<(num_images*12) & currI>2000 & currcrlb<0.05;
zd = double(currzfresult(mask));
pd = double(currzang(mask));
xd = double(currx(mask));
yd = double(curry(mask));

if lp0==0
    [resultp, lp, resnormx] = fitPhaseXYZ(xd,yd,zd,pd,[0.02,0,0,0]);
else
    lp=lp0;
end

% pd_corr = checkPhase(PhaseXYCorrection(pd, lp, xd, yd, [0 0]));
% resultp_corr = checkPhase(PhaseXYCorrection(resultp, lp, xd, yd, [0 0]));
pdata_corr1 = checkPhase(PhaseXYCorrection(currzang, lp, currx, curry, [0 0]));
currzang_cor = pdata_corr1;

[xx, yy] = meshgrid(1:256,1:256);
ret = PhaseXYCorrection(zeros(size(xx)), lp, xx, yy, [0 0]);
corrvalue = max(ret(:)) - min(ret(:))

figure_flag = 1;
doXYcorr = 1;
if figure_flag>0
    figure();
%     set(gcf,'position',[150,0,1200,800]);
%     subplot(1,2,1);
%     [xx, yy] = meshgrid(1:250,1:250);
%     ret = PhaseXYCorrection(zeros(size(xx)), lp, xx, yy, [0 0]);
    imagesc(ret)
    colorbar
    if doXYcorr>0
        title('step1: XY corr');
    else
        title('step1: XY corr (SKIPPED)');
    end
end

end

