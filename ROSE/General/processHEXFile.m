function [smInfo, procresult] = processHEXFile(filename, IsEM, xydrift, forcefit)
% [smInfo, procresult] = processHEXFile(filename, IsEM, xydrift, forcefit)
%Process file and save the results
%smInfo:
% X, Y,newX, newY, Frame, std, int1, int2, int3, int4 ,int5 ,int6 ,phase1,    phase2, amp1, anp2, bkg1, bkg2
%1   2   3    4    5       6    7     8     9     10    11    12    13          14       15   16     17    18
%
%xydrift for Galvo drift compensation
%forcefit = 0(default) for load detection and fitting from saved mat files
%                       if possible
%forcefit = 1 for calculating detection and fitting regardless of mat files
%procresult:
% procresult.filename = filename;
% procresult.x = smInfo_x;
% procresult.y = smInfo_y;
% procresult.xraw = smInfo(:,1);
% procresult.yraw = smInfo(:,2);
% procresult.std = smInfo_std;
% procresult.p1 = smInfo_p1;
% procresult.p2 = smInfo_p2;
% procresult.int = smInfo_int;
% procresult.int1 = smInfo_int1;
% procresult.int2 = smInfo_int2;
% procresult.md1 = smInfo_md1;
% procresult.md2 = smInfo_md2;
% procresult.rx = smInfo_rx;
% procresult.ry = smInfo_ry;
% %phase plane
% procresult.lp1 = lp1;
% procresult.lp2 = lp2;
% procresult.dx = smInfo_dx;
% procresult.dy = smInfo_dy;
% procresult.A = A;
% procresult.resultp1 = resultp1;
% procresult.resultp2 = resultp2;
% procresult.phase2pos_func = phase2pos_func;
%
%% parameters
%detection
detection_threshold = 3;
detection_size = 17;
detection_disk_radius = 4.00; %4.75;
%fit
fit_r = 4.0;
%subimage generation
subimghw = 8; %-w~+w
%post processing fit data
rotation = 90;%rotation of XY
% misc
gain = 100;%100
e2count_factor = 10;%10


if nargin < 2 || isempty(IsEM)
    warning('parameter IsEM is empty, use default 1(EM on)');
    IsEM = 1;
end
% xydrift for galvo drift compensation
if nargin < 3 
    xydrift = [];
end

if nargin < 4 || isempty(forcefit)
    forcefit = 0;
end


tpos = strfind(filename, '\');
foldername = filename(1:tpos(end));
matfilename = [filename(1:end-4) '_procdata*.mat'];
matfilelist = dir(matfilename);
if ~isempty(matfilelist)
    tmatfile = matfilelist(1).name;
else
    tmatfile = '';
    forcefit = 1;
end
if len(matfilelist)>1
    disp(['Multi mat file using: ' tmatfile])
end
if forcefit == 0
    tmatfile = [foldername tmatfile];
end
%% load image
if forcefit == 0 && exist(tmatfile,'file')
    imgbuf = tiffread(filename, 1);
    tempbuf = load(tmatfile, 's');
    s = tempbuf.s;
else
%     imgbuf = tiffread(filename);
    imgbuf = LoadTiff16bit(filename);
    s = size(imgbuf);
end

%% detection
% detection_result=[];
if forcefit == 0 && exist(tmatfile,'file')%load data
    disp('Loading detection data');
    tempbuf = load(tmatfile, 'detection_result');
    detection_result = tempbuf.detection_result;
else
    tic
    [detection_result] = HEXDetection3(imgbuf, detection_threshold, detection_size, detection_disk_radius);
    toc
end

%% Post-processing & Duplicate sub-images
moleculecnt = 0;
if forcefit == 0 && exist(tmatfile,'file')%load data
    %moleculecnt
    tempbuf = load(tmatfile, 'moleculecnt');
    moleculecnt = tempbuf.moleculecnt;
    %subimginfo
    tempbuf = load(tmatfile, 'subimginfo');
    subimginfo = tempbuf.subimginfo;
else
    for m=1:length(detection_result)
        moleculecnt = moleculecnt + size(detection_result{m},1);
    end

    subimgsize = subimghw*2+1;
    % s = size(imgbuf);
    subimgbuf = zeros(subimgsize,subimgsize,moleculecnt);
    subimginfo = zeros(moleculecnt, 3);%x,y,frameidx
    subimglen = 0;
    for m=1:length(detection_result)
        tbuf = detection_result{m};
        for n=1:size(tbuf,1)
            tx = tbuf(n,1);
            ty = tbuf(n,2);
            txs = tx - subimghw;
            txe = tx + subimghw;
            if txs <1
                txs = 1;
                txe = txs + subimghw + subimghw;
            end
            if txe >s(2)
                txe = s(2);
                txs = txe - subimghw - subimghw;
            end

            tys = ty - subimghw;
            tye = ty + subimghw;
            if tys <1
                tys = 1;
                tye = tys + subimghw + subimghw;
            end
            if tye >s(1)
                tye = s(1);
                tys = tye - subimghw - subimghw;
            end

            subimg = imgbuf(tys:tye, txs:txe, m);

            subimglen = subimglen+1;
            subimgbuf(:,:,subimglen) = subimg;
            subimginfo(subimglen, 1) = tx;
            subimginfo(subimglen, 2) = ty;
            subimginfo(subimglen, 3) = m;
        end
    end
    moleculecnt = subimglen;
end
fprintf('Detected Counts: %d\n', moleculecnt);

%% fitting
if forcefit == 0 && exist(tmatfile,'file')%load data
    disp('Loading Fitting data');
    tempbuf = load(tmatfile, 'fitresult','fitrawdata','phasedata');
    fitresult = tempbuf.fitresult;
    fitrawdata = tempbuf.fitrawdata;
    phasedata = tempbuf.phasedata;
else
    tic
    disp(['Fitting.. (mcnt=' num2str(moleculecnt) ')']);
    [fitresult, fitrawdata, phasedata] = fitHEX(subimgbuf, fit_r);
    %[result fitdata phasedata] = fitHEX(imgbuf)
    %result: [x y p1 p2 std int md1 md2]
    %         1 2 3   4  5   6   7  8  
    %fitdata: [x0 y0 std int1 int2 int3 int4 int5 int6 bkg rx ry phase_offset]
    %          1  2  3   4    5    6    7     8    9    10 11 12     13
    %phasedata: [phase1, amp1, offset1, phase2, amp2, offset2]
    %              1      2       3       4      5      6 
    fitmask = fitresult(:,7)>0 & fitresult(:,7)<2 & fitresult(:,8)>0 & fitresult(:,8)<2;
    fitresult = fitresult(fitmask,:);
    fitrawdata = fitrawdata(fitmask,:);
    phasedata = phasedata(fitmask,:);
    subimginfo = subimginfo(fitmask,:);
    moleculecnt = sum(fitmask);
    toc
end

%% PostProcessing
theta = rotation/180*pi;
rotMat = [cos(theta) -sin(theta); sin(theta) cos(theta)];
smInfo = zeros(moleculecnt, 18);
% X, Y,newX, newY, Frame, std, int1, int2, int3, int4 ,int5 ,int6 ,phase1,    phase2, amp1, anp2, bkg1, bkg2
%1   2   3    4    5       6    7     8     9     10    11    12    13          14       15   16     17    18
smInfo(:,1) = fitresult(:,1) + subimginfo(:,1);
smInfo(:,2) = fitresult(:,2) + subimginfo(:,2);
smInfo(:,3:4) = smInfo(:,1:2) * rotMat;
smInfo(:,5) = subimginfo(:,3);
smInfo(:,6) = fitresult(:,5);
smInfo(:,7:12) = fitrawdata(:,4:9);
smInfo(:,13) = phasedata(:,1);
smInfo(:,14) = phasedata(:,4);
smInfo(:,15) = phasedata(:,2);
smInfo(:,16) = phasedata(:,5);
smInfo(:,17) = phasedata(:,3);
smInfo(:,18) = phasedata(:,6);

%% Post Processing
smInfo_frame = smInfo(:,5);
smInfo_p1 = smInfo(:,13);
smInfo_p2 = smInfo(:,14);
%Photon Number
smInfo_int1 = sum(smInfo(:,7:9),2)/gain*e2count_factor;%photon number of Phase 1
smInfo_int2 = sum(smInfo(:,10:12),2)/gain*e2count_factor;%photon number of Phase 2
smInfo_int = smInfo_int1 + smInfo_int2;%photon number
%Modulatioin Depth
smInfo_md1 = smInfo(:,15)./smInfo(:,17);
smInfo_md2 = smInfo(:,16)./smInfo(:,18);
%r
smInfo_rx = fitrawdata(:,11);
smInfo_ry = fitrawdata(:,12);
smInfo_std = smInfo(:,6);

%% fit phase plane
disp('---------- Fitting Phase Plane ----------');
tx = smInfo(:,1);
ty = smInfo(:,2);
p1 = smInfo(:,13);
p2 = smInfo(:,14);
if ~isempty(xydrift) && size(xydrift,1) >= max(smInfo(:,5))
    framelist = smInfo(:,5);
    xydriftdata = xydrift(framelist,:);
    tx = tx - xydriftdata(:,1);
    ty = ty - xydriftdata(:,2);
end
[lp1, lp2, resultp1, resultp2, resultx, resulty, resultdx, resultdy, A] ...
    = CalculatePhaseAndFineXY_faster(tx, ty, p1, p2, IsEM);

%function for phase-position translate
phase2pos_func = @(p1, p2)[p1(:), p2(:)]/A;
%--------- smInfo_x and smInfo_y as resultx and resulty ----------%
smInfo_x = resultx;
smInfo_y = resulty;
%dx and dy
smInfo_dx = resultdx;
smInfo_dy = resultdy;

%% clean data
procresult = [];
procresult.filename = filename;
procresult.s = s;
procresult.x = smInfo_x;
procresult.y = smInfo_y;
procresult.frame = smInfo_frame;
procresult.xraw = smInfo(:,1);
procresult.yraw = smInfo(:,2);
procresult.std = smInfo_std;
procresult.p1 = smInfo_p1;
procresult.p2 = smInfo_p2;
procresult.int = smInfo_int;
procresult.int1 = smInfo_int1;
procresult.int2 = smInfo_int2;
procresult.md1 = smInfo_md1;
procresult.md2 = smInfo_md2;
procresult.rx = smInfo_rx;
procresult.ry = smInfo_ry;
%phase plane
procresult.lp1 = lp1;
procresult.lp2 = lp2;
procresult.dx = smInfo_dx;
procresult.dy = smInfo_dy;
procresult.A = A;
procresult.resultp1 = resultp1;
procresult.resultp2 = resultp2;
procresult.phase2pos_func = phase2pos_func;
%% save results
fighandle = figure();
subplot(3,3,1)
hist(smInfo_md1,50)
title('MD of Phase1')
subplot(3,3,4)
hist(smInfo_md2,50)
title('MD of Phase2')

subplot(3,3,2)
hist(smInfo_int,50)
title('Phton Number')

subplot(3,3,5)
plot(smInfo_int1, smInfo_int2, '.')
title('Int1 vs. Int2')

subplot(3,3,3)
hist(smInfo_std,50)
title('std of PSF')

subplot(3,3,6)
hist([smInfo_rx smInfo_ry],20)
title('r of PSF')

subplot(3,3,7)
hist(resultdx,50);
title('DX')

subplot(3,3,8)
hist(resultdy,50);
title('DY')

subplot(3,3,9)
hist((smInfo_int1-smInfo_int2)./(smInfo_int1+smInfo_int2), 50);
title('Int Ratio')

saveFigure(fighandle.Number, [filename(1:end-4) '_result.png'])
%% save data
if forcefit == 1 || ~exist(tmatfile,'file')%load data
    t = clock();
    clear imgbuf
    % clear subimgbuf
    clear fighandle
    save(sprintf('%s_procdata_%d-%d-%d_%d-%d-%.0f.mat',filename(1:end-4), t(1),t(2),t(3),t(4),t(5),t(6)));
end

end