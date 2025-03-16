%load('Y:\ROSE\20190124_10grid\sample7\spool_7189_Data.mat')
% tx = smInfo(:,1);
% ty = smInfo(:,2);
% p1 = smInfo(:,5);
% p2 = smInfo(:,6);
% [lp1, lp2, resultp1, resultp2, resultx, resulty, resultdx, resultdy, A] ...
% = CalculatePhaseAndFineXY_faster_6(tx, ty, p1, p2);
[tp1, tidx] = sort(p1);
tdy = resultdy(tidx);
[tp2, tidx2] = sort(p2);
tdx = resultdx(tidx2);

th = figure();
plot(tp1, smooth(tdy, 200),'b')
hold on
plot(tp2, smooth(tdx, 200),'r')
hold off
legend phase1 phase2
xlabel phase
ylabel dxy

saveas(th, [detFileInfo(1:end-22) '_result_dp.png']);
