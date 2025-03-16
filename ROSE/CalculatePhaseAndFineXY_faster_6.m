function [lp1, lp2, resultp1, resultp2, resultx, resulty, resultdx, resultdy, A] ...
    = CalculatePhaseAndFineXY_faster_6(x, y, p1, p2)
%% Calculate phase plane and  calculate the fine coordinates X and Y
%data : x y p1 p2
%results: lp1, lp2, resultp1, resultp2, resultx, resulty, 
%         resultdx, resultdy
%parameters: samplelength
% x = smInfo(:,1);
% y = smInfo(:,2);
% p1 = smInfo(:,13);
% p2 = smInfo(:,14);

% if nargin < 5 || isempty(IsEM)
%     disp('parameter IsEM is empty, use default 1(EM On)')
%     IsEM = 1;
% end

% if nargin < 6 || isempty(samplelength)
    samplelength = 10000;
% end


%% fitting phase plane
disp('Fitting Phase Plane 1 ..');
[lp1, resultp1] = SearchPhasePlane6(x,y,p1,samplelength,[-0.1 0.1]); % -3.2 -3.1 % -0.05 0.03 640laser & 4.8-560pupil 
title('Fitting P1');

disp('Fitting Phase Plane 2 ..'); 
[lp2, resultp2] = SearchPhasePlane6(x,y,p2,samplelength,[-1.7 -1.5]);%[-1.7 -1.4] for EM, [1.55 1.65] for Conventional; -1.65 -1.55 for 640  [-1.7 -1.5]); % Cell09: [-1.7 -1.5]
title('Fitting P2');
drawnow

%% calculate fine X and Y
dp1 = checkPhase(resultp1 - p1);
dp2 = checkPhase(resultp2 - p2);
a1 = lp1(1); % 每个点都计算出 相位和 位置的关系 a 和 b
b1 = lp1(2);
a2 = lp2(1);
b2 = lp2(2);
A = [b1*cos(a1), b1*sin(a1); b2*cos(a2), b2*sin(a2)];
B = [dp1, dp2];
% iA = inv(A);
% x * A = B   =>   x = B * inv(A) = B/A
resX = B/A;
resultdx = resX(:,1);
resultdy = resX(:,2);
resultx = x - resultdx;
resulty = y - resultdy;

