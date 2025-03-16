function [fitresult, smInfo] = SMFitting(imglist, detfile, isEM)
%fitresult:    [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
%              [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]

%smInfo: [x y frame int phase1 phase2 md1 md2  idx]
%        [1 2   3    4    5      6     7   8    9 ]
if nargin <3
    warning('isEM is NOT set, using default: 1');
    isEM = 1;
end

%% Single Molecule Fitting 
intCompensationData = [1.0000    1.132    1.206    1.0000    1.054    1.002]; %intesity compensation
imgfile1 = imglist{1};
imgfile2 = imglist{2};
imgfile3 = imglist{3};
imgfile4 = imglist{4};
imgfile5 = imglist{5};
imgfile6 = imglist{6};


subimgsize = 7;
maxPhoton = 50000;
maxstd = 1.5;
minstd = 0.5;

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
            subimginfo(subimgcnt,:) = [tx ty m];
        end
    end
end

subimgbuf = subimgbuf(:,:,:,1:subimgcnt);
subimginfo = subimginfo(1:subimgcnt,:);

%% fitting sub images
%result:    [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
%           [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]
fitresult = fit6Gauss_mp(subimgbuf);

%% post process
%smInfo
%[x y frame int p1 p2 md1 md2]
%
gain = 100;
e2cfactor = 10;
e2cfactor_EMoff = 1.85;%e 2 c factor when EM off

phase0 = 120/180*pi;
phasedata = zeros(subimgcnt, 8); %[phase1 phase2 md1 md2 amp1 amp2 offset1 offset2]
for m=1:subimgcnt
    intlist = fitresult(m,5:10).*intCompensationData;
    [phase1, amp1, offset1] = MyCalPhase3_new(intlist(1:3), phase0);
    [phase2, amp2, offset2] = MyCalPhase3_new(intlist(4:6), phase0);
    phasedata(m,:) = [phase1, phase2, amp1/offset1, amp2/offset2, amp1, amp2, offset1, offset2];
end

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

smInfomask = fitresult(:,1)<3 & fitresult(:,1)>-3 & fitresult(:,2)<3 & fitresult(:,2)>-3 & ...%fit mask
            smInfo(:,4)>0 & smInfo(:,4)<maxPhoton & ...
            fitresult(:,3) >= minstd & fitresult(:,3) <= maxstd & fitresult(:,4) >= minstd & fitresult(:,4) <= maxstd & ...
            smInfo(:,7)>0.5 & smInfo(:,7)<1.3 & smInfo(:,8)>0.5 & smInfo(:,8)<1.3;  % md mask

smInfo = smInfo(smInfomask, :);
