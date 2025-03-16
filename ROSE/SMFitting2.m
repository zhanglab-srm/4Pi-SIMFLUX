function [fitresult, smInfo, intCompensationData] = SMFitting2(imglist, detfile, isEM, gain, doIntComp)
%fitresult:    [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
%              [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]

%smInfo: [x y frame int phase1 phase2 md1 md2  idx]
%        [1 2   3    4    5      6     7   8    9 ]
if nargin <3
    warning('isEM is NOT set, using default: 1');
    isEM = 1;
end
if nargin <4
    warning('EMgain is NOT set, using default: 100');
    gain = 100;
end
if nargin <5
    warning('doIntComp is NOT set, using default: 1');
    doIntComp = 1;
end
%% Single Molecule Fitting 
% intCompensationData = [1.0000    1.132    1.206    1.0000    1.054    1.002]; %intesity compensation
imgfile1 = imglist{1};
imgfile2 = imglist{2};
imgfile3 = imglist{3};
imgfile4 = imglist{4};
imgfile5 = imglist{5};
imgfile6 = imglist{6};

subimgsize = 7;
maxPhoton = 50000; % ~ROSE 50000
maxstd = 1.5;
minstd = 0.5; % 0.5 G

%% load images
imgbuf1 = LoadTiff16bitRaw(imgfile1);
imgbuf2 = LoadTiff16bitRaw(imgfile2);
imgbuf3 = LoadTiff16bitRaw(imgfile3);
imgbuf4 = LoadTiff16bitRaw(imgfile4);
imgbuf5 = LoadTiff16bitRaw(imgfile5);
imgbuf6 = LoadTiff16bitRaw(imgfile6);

%% load detection results
tempdata = load(detfile);
detectResult = tempdata.detectResult;
detectCnt = tempdata.detectCnt;
detectCntList = tempdata.detectCntList;

%% subimg
s = size(imgbuf1);
subimgbuf = zeros(subimgsize,subimgsize,6,detectCnt);
subimginfo = zeros(detectCnt, 3);%[x y frame]
windowsize = round((subimgsize-1)/2);

disp('Make sub imges ..');

subimgcnt = 0;
for m=1:len(detectResult)
    temp = round(detectResult{m});%[x y ~]
    tlen = size(temp,1);
    for n=1:tlen
        tx = temp(n,1);
        ty = temp(n,2);
        xs = tx - windowsize;
        xe = tx + windowsize;
        ys = ty - windowsize;
        ye = ty + windowsize;
        if xs>0 && ys >0 && xe<=s(2) && ye<=s(1)
            subimgcnt = subimgcnt +1;
            subimgbuf(:,:,1,subimgcnt) = imgbuf1(ys:ye, xs:xe, m);
            subimgbuf(:,:,2,subimgcnt) = imgbuf2(ys:ye, xs:xe, m);
            subimgbuf(:,:,3,subimgcnt) = imgbuf3(ys:ye, xs:xe, m);
            subimgbuf(:,:,4,subimgcnt) = imgbuf4(ys:ye, xs:xe, m);
            subimgbuf(:,:,5,subimgcnt) = imgbuf5(ys:ye, xs:xe, m);
            subimgbuf(:,:,6,subimgcnt) = imgbuf6(ys:ye, xs:xe, m);
            subimginfo(subimgcnt,:) = [tx ty m]; % 第几个子图 它的 x y 坐标 和属于哪一张 z
        end
    end
end
% 把探测到的点 应该是最亮点 + - window，框出一小块来

subimgbuf = subimgbuf(:,:,:,1:subimgcnt);  % [x框，  y框，  1-6幅图，  1-所有子图]
subimginfo = subimginfo(1:subimgcnt,:); % [1-所有子图， 检测x坐标， 检测y坐标，  z坐标] 
 
% 去掉6幅图中，任一张图光强小于100的点
% for i=1:subimgcnt
%     if subimgbuf(4,4,1,i)>=200 && subimgbuf(4,4,2,i)>=200 && subimgbuf(4,4,3,i)>=200 && subimgbuf(4,4,4,i)>=200 && subimgbuf(4,4,5,i)>=200 && subimgbuf(4,4,6,i)>=200
%         subimgbuf(:,:,:,i)=subimgbuf(:,:,:,i);
%         subimginfo(i,:)=subimginfo(i,:);
%     else
%         subimgbuf(:,:,:,i)=[];
%         subimginfo(i,:)=[];
%     end 
%     disp (i)
% end

% tiffwrite(subimgbuf(:,:,:,i),'0.tif');

%% fitting sub images
%result:    [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
%           [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]
fitresult = fitGauss6GPU(subimgbuf); % fitGauss6GPU fit6Gauss4_mp  %输入 [x框，  y框，  1-6幅图，  1-所有子图]  输出 【37346,16】

%% estimate intensity compensation data
if doIntComp>0
    EstIngCompRatio = 5000;
    lp1 = EstIntComp(fitresult(:,5:7), EstIngCompRatio);
    lp2 = EstIntComp(fitresult(:,8:10), EstIngCompRatio);
    intCompensationData = [1, lp1(1), lp1(2), 1, lp2(1), lp2(2)];
    fprintf('Int Compensation data: [%0.2f %0.2f %0.2f %0.2f]\n', lp1(1), lp1(2), lp2(1), lp2(2));
else
%     intCompensationData = [1, 1, 1, 1, 1, 1];
    intCompensationData = [1, 1, 1, 1, 1, 1];
%     intCompensationData = [1.0000  1.0906  1.1364  1.0000  1.1107  0.9916]; % F:\ROSE\20190219_10grid\sample2_reg\ ~ROSE
%     intCompensationData = [1.0000  1.0664  1.1022  1.0000  1.1299  0.9810]; % F:\ROSE\20190215_5Lgrid\sample2_reg\
%     intCompensationData = [1.0000  1.1340  1.1148  1.0000  1.0772  0.9682]; % E:\ROSE\20190131_IBP\sample4_reg\
%     intCompensationData = [1.0000    1.1178    1.1746    1.0000    1.1495    1.0117]; % E:\ROSE\20190131_5Lgrid\sample2_reg
%     intCompensationData = [1.0000    1.0927    1.1514    1.0000    1.1158    0.9907]; % E:\ROSE\20190131_5Lgrid\sample2_reg
%     intCompensationData = [1.0000    1.0607    1.1201    1.0000    1.0817    0.9701]; % E:\ROSE\20190130_5Lgrid\sample6_reg
%     intCompensationData = [1.0000    1.0308    1.0912    1.0000    1.0454    0.9462]; % E:\ROSE\20190130_5Lgrid\sample4_reg
%     intCompensationData = [1.0000    1.0194    1.0864    1.0000    1.0397    0.9473]; % E:\ROSE\20190130_5Lgrid\sample2_reg
    disp('NO int comp applied');
end

%% post process
%smInfo
%[x y frame int p1 p2 md1 md2]
%
% gain = 100;
e2cfactor = 10.5;
e2cfactor_EMoff = 1;%e 2 c factor when EM off

phase0 = 120/180*pi;
phasedata = zeros(subimgcnt, 8); %[phase1 phase2 md1 md2 amp1 amp2 offset1 offset2]
maskI = zeros(subimgcnt,2);
for m=1:subimgcnt
    intlist = fitresult(m,5:10).*intCompensationData;
    [phase1, amp1, offset1] = MyCalPhase3_new(intlist(1:3), -phase0); % 4PX -phase0
    [phase2, amp2, offset2] = MyCalPhase3_new(intlist(4:6), phase0);
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
if isEM >0  %EM on
    smInfo(:,4) = sum(fitresult(:,5:10), 2).*2.*pi.*fitresult(:,3).*fitresult(:,4)./gain.*e2cfactor;
else    %EM off
    smInfo(:,4) = sum(fitresult(:,5:10), 2).*2.*pi.*fitresult(:,3).*fitresult(:,4).*e2cfactor_EMoff;
end
smInfo(:,5:8) = phasedata(:,1:4);
smInfo(:,9) = 1:subimgcnt;

% smInfomask = fitresult(:,1)<3 & fitresult(:,1)>-3 & fitresult(:,2)<3 & fitresult(:,2)>-3 & ...%fit mask
%             smInfo(:,4)>-100 & smInfo(:,4)<maxPhoton & ...
%             fitresult(:,3) >= minstd & fitresult(:,3) <= maxstd & fitresult(:,4) >= minstd & fitresult(:,4) <= maxstd & ...
%             smInfo(:,7)>0.2 & smInfo(:,7)<1.7 & smInfo(:,8)>0.2 &
%             smInfo(:,8)<1.7;  % md mask ~ROSE

smInfomask = fitresult(:,1)<1 & fitresult(:,1)>-1 & fitresult(:,2)<1 & fitresult(:,2)>-1 & ...%fit mask
            smInfo(:,4)>500 & smInfo(:,4)<maxPhoton & ...
            fitresult(:,3) >= minstd & fitresult(:,3) <= maxstd & fitresult(:,4) >= minstd & fitresult(:,4) <= maxstd & ...
            smInfo(:,7)>0.5 & smInfo(:,7)<1.2 & smInfo(:,8)>0.5 & smInfo(:,8)<1.2 & maskI(:,1)>0 & maskI(:,2)>0; % md mask & maskI(:,1)>0 & maskI(:,2)>0

% 根据干涉对比度强弱去掉一些点   ;  % md mask & maskI(:,1)>0 & maskI(:,2)>0

%%
% smInfo=smInfo; % 不加mask
smInfo = smInfo(smInfomask, :);
% fitresult = fitresult(smInfomask, :);
end
