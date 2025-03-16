%% display DAFL results
%fitresult:    [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
%              [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]

%smInfo: [x y frame int phase1 phase2 md1 md2  idx]
%        [1 2   3    4    5      6     7   8    9 ]

%befoe grid fitting
% xcorr, ycorr, xcorr2, ycorr2
%after grid fitting
% x_gcorr, y_gcorr, x_gcorr2, y_gcorr2
%

disprange = [80 160 80 160];

hf = figure();
% ROSE before grid
dispx = xcorr;
dispy = ycorr;
dispf = smInfo(:,3);
dmask = dispx>disprange(1) & dispx<disprange(2) & dispy > disprange(3) & dispy<disprange(4);

subplot(2,2,1)
[dlist2, max_x] = doDAFL(dispx(dmask), dispy(dmask),dispf(dmask),0);
title(sprintf('ROSE DAFL, std = %0.2f', max_x/sqrt(2)));



% Gaussian before grid
dispx = xcorr2;
dispy = ycorr2;
dispf = smInfo(:,3);
dmask = dispx>disprange(1) & dispx<disprange(2) & dispy > disprange(3) & dispy<disprange(4);

subplot(2,2,2)
[dlist2, max_x] = doDAFL(dispx(dmask), dispy(dmask),dispf(dmask),0);
title(sprintf('Gaussian DAFL, std = %0.2f', max_x/sqrt(2)));

if doGrid>0
% ROSE with grid
dispx = x_gcorr;
dispy = y_gcorr;
dispf = smInfo(mask,3);
dmask = dispx>disprange(1) & dispx<disprange(2) & dispy > disprange(3) & dispy<disprange(4);

subplot(2,2,3)
[dlist2, max_x] = doDAFL(dispx(dmask), dispy(dmask),dispf(dmask),0);
title(sprintf('ROSE Grid DAFL, std = %0.2f', max_x/sqrt(2)));

% Gaussian with grid
% dispx = x_gcorr2;
% dispy = y_gcorr2;
% dispf = smInfo(mask,3);
% dmask = dispx>disprange(1) & dispx<disprange(2) & dispy > disprange(3) & dispy<disprange(4);
% 
% subplot(2,2,4)
% [dlist2, max_x] = doDAFL(dispx(dmask), dispy(dmask),dispf(dmask),0);
% title(sprintf('Gaussian Grid DAFL, std = %0.2f', max_x/sqrt(2)));
end
%% save as img
saveas(hf, [detFileInfo(1:end-22) '_result2.png']);