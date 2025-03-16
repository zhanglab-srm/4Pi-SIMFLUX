function [fitresult,smInfo,smInfomask] = ROSE_SMFitting(subimgbuf,subimginfo)

fitresult = fitGauss6GPU(subimgbuf);

subimgcnt = size(subimgbuf,4);
intCompensationData = [1, 1, 1, 1, 1, 1];

phase0 = 120/180*pi;
phasedata = zeros(subimgcnt, 8); %[phase1 phase2 md1 md2 amp1 amp2 offset1 offset2]
maskI = zeros(subimgcnt,2);
for m=1:subimgcnt
    intlist = fitresult(m,5:10).*intCompensationData;
    [phase1, amp1, offset1] = MyCalPhase3_new(intlist(4:6), -phase0); % 4PX -phase0
    [phase2, amp2, offset2] = MyCalPhase3_new(intlist(1:3), phase0);
    phasedata(m,:) = [phase1, phase2, amp1/offset1, amp2/offset2, amp1, amp2, offset1, offset2];
    if intlist(1)>10 && intlist(6)>10 && all(intlist>0)
        maskI(m,1) = 1;
    end
    if abs((sum(intlist(1:3))-sum(intlist(4:6)))/sum(intlist))<0.2
        maskI(m,2) = 1;
    end
end
% 此时根据三幅图强度关系计算出的相位，被包裹

smInfo = zeros(subimgcnt, 9); 
%[x y frame int phase1 phase2 md1 md2  idx]
%[1 2   3    4    5      6     7   8    9 ]
smInfo(:,1) = fitresult(:,1) + subimginfo(:,1);
smInfo(:,2) = fitresult(:,2) + subimginfo(:,2);
smInfo(:,3) = subimginfo(:,3);
smInfo(:,4) = sum(fitresult(:,5:10), 2).*2.*pi.*fitresult(:,3).*fitresult(:,4);
smInfo(:,5:8) = phasedata(:,1:4);
smInfo(:,9) = 1:subimgcnt;
smInfo(:,10:13) = 0;
smInfo(:,14:19) = fitresult(:,5:10);

% smInfomask = fitresult(:,1)<3 & fitresult(:,1)>-3 & fitresult(:,2)<3 & fitresult(:,2)>-3 & ...%fit mask
%             smInfo(:,4)>-100 & smInfo(:,4)<maxPhoton & ...
%             fitresult(:,3) >= minstd & fitresult(:,3) <= maxstd & fitresult(:,4) >= minstd & fitresult(:,4) <= maxstd & ...
%             smInfo(:,7)>0.2 & smInfo(:,7)<1.7 & smInfo(:,8)>0.2 &
%             smInfo(:,8)<1.7;  % md mask ~ROSE

% smInfomask = fitresult(:,1)<1 & fitresult(:,1)>-1 & fitresult(:,2)<1 & fitresult(:,2)>-1 & ...%fit mask
%             smInfo(:,4)>500 & ...
%             fitresult(:,3) >= minstd & fitresult(:,3) <= maxstd & fitresult(:,4) >= minstd & fitresult(:,4) <= maxstd & ...
%             smInfo(:,7)>0.5 & smInfo(:,7)<1.2 & smInfo(:,8)>0.5 & smInfo(:,8)<1.2 & maskI(:,1)>0 & maskI(:,2)>0; % md mask & maskI(:,1)>0 & maskI(:,2)>0

% smInfomask = smInfo(:,7)>0.5 & smInfo(:,7)<1.2 & smInfo(:,8)>0.5 & smInfo(:,8)<1.2 & maskI(:,1)>0 & maskI(:,2)>0 & fitresult(:,1)<3 & fitresult(:,1)>-3 & fitresult(:,2)<3 & fitresult(:,2)>-3;

r = size(subimgbuf,1);
hd = (r-1)/2;
smInfomask = abs(fitresult(:,1)-hd)<3 & abs(fitresult(:,2)-hd)<3;
% 根据干涉对比度强弱去掉一些点   ;  % md mask & maskI(:,1)>0 & maskI(:,2)>0

%%
% smInfo=smInfo; % 不加mask
% smInfo = smInfo(smInfomask, :);
% fitresult = fitresult(smInfomask, :);

end