%% display result
%fitresult:    [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
%              [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]

%smInfo: [x y frame int phase1 phase2 md1 md2  idx]
%        [1 2   3    4    5      6     7   8    9 ]

smInfoD_int1 = sum(fitresult(:,5:7),2);
smInfoD_int2 = sum(fitresult(:,8:10),2);
smInfoD_int = smInfo(:,4);
smInfoD_intR = (smInfoD_int1 - smInfoD_int2) ./ (smInfoD_int1 + smInfoD_int2);

smInfoD_md1 = smInfo(:,7);
smInfoD_md2 = smInfo(:,8);
smInfoD_wx = fitresult(:,3);
smInfoD_wy = fitresult(:,4);

smInfoD_dx = smInfo_dx;
smInfoD_dy = smInfo_dy;

hf = figure();
smInfoD_int = smInfoD_int(smInfoD_int>0 & smInfoD_int<40000);
subplot(3,3,1)
hist(smInfoD_int,20);
title(sprintf('Photon mean = %0.1f',mean(smInfoD_int)));


smInfoD_intR = smInfoD_intR(smInfoD_intR>-1 & smInfoD_intR<1);
subplot(3,3,2)
hist(smInfoD_intR,20);
title(sprintf('intR mean = %0.2f',mean(smInfoD_intR)));

smInfoD_md1 = smInfoD_md1(smInfoD_md1>0.2 & smInfoD_md1<1.7);
subplot(3,3,3)
hist(smInfoD_md1,20);
title(sprintf('MD1 mean = %0.2f\nstd=%0.2f',mean(smInfoD_md1), std(smInfoD_md1)));

smInfoD_md2 = smInfoD_md2(smInfoD_md2>0.2 & smInfoD_md2<1.7);
subplot(3,3,4)
hist(smInfoD_md2,20);
title(sprintf('MD2 mean = %0.2f\nstd=%0.2f',mean(smInfoD_md2), std(smInfoD_md2)));


smInfoD_wx = smInfoD_wx(smInfoD_wx>0.4 & smInfoD_wx<1.7);
subplot(3,3,5)
hist(smInfoD_wx,20);
title(sprintf('wx mean = %0.2f',mean(smInfoD_wx)));

smInfoD_wy = smInfoD_wy(smInfoD_wy>0.4 & smInfoD_wy<1.7);
subplot(3,3,6)
hist(smInfoD_wy,20);
title(sprintf('wy mean = %0.2f',mean(smInfoD_wy)));


smInfoD_dx = smInfoD_dx(smInfoD_dx>-1 & smInfoD_dx<1);
subplot(3,3,7)
hist(smInfoD_dx,20);
title(sprintf('dx'));

smInfoD_dy = smInfoD_dy(smInfoD_dy>-1 & smInfoD_dy<1);
subplot(3,3,8)
hist(smInfoD_dy,20);
title(sprintf('dy'));

subplot(3,3,9)
doDAFL(xcorr,ycorr,smInfo(:,3),0);

%% save as img
saveas(hf, [detFileInfo(1:end-22) '_result.png']);

hf1=figure();
doDAFL(xcorrg,ycorrg,smInfo(:,3),0);
saveas(hf1, [detFileInfo(1:end-22) '_DAFLg.png']);