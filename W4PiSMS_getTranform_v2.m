function W4PiSMS_getTranform_v2(handles)

close all
posfile = get(handles.pathMainfolder, 'String');
center1 = str2double(get(handles.center1,'string'));
center2 = str2double(get(handles.center2,'string'));
center3 = str2double(get(handles.center3,'string'));
center4 = str2double(get(handles.center4,'string'));
centers = [center1 center2 center3 center4];

c = posfile(end);       %get last character of path to image file
while c ~= '\'          %get path only
    posfile(end) = [];  %delete filename character
    if isempty(posfile) %break if no path present
        break;
    end
    c = posfile(end);   %get last character of path to image file
end

[imageName,dataFolder] = uigetfile([posfile '*.dcimg'],'Open Calibration Images','MultiSelect','on');

if ~isequal(imageName,0)    
    if ~iscell(imageName)
        imageName={imageName};
    end
    L=numel(imageName);
    
    if contains(imageName{1},'_642_')        
        resultpath = [pwd '\calibration_files\'];
        namestr = ['align' '_642_'];
        mkdir(resultpath);
        channel = '642';
    elseif contains(imageName{1},'_561_')
        resultpath = [pwd '\calibration_files\'];
        namestr = ['align' '_561_'];
        mkdir(resultpath);
        channel = '561';
    else
        x = inputdlg({'Result Path','channel'},'File Directory',2,{[pwd '\calibration_files\'],''},options);
        resultpath = x{1};
        channel= x{2};
        namestr = ['align' '_' channel '_'];
        mkdir(resultpath);
    end
else
    return
end

qd1=[];qd2=[];qd3=[];qd4=[];
for ii=1:L
    [~,qds]=W4PiSMS_readdcimg([dataFolder imageName{ii}],centers);
    meanqds1 = mean(qds(:,:,:,1),3);
    meanqds2 = mean(qds(:,:,:,2),3);
    meanqds3 = mean(qds(:,:,:,3),3);
    meanqds4 = mean(qds(:,:,:,4),3);
    qd1(:,:,ii)=meanqds1-median(meanqds1(:));
    qd2(:,:,ii)=meanqds2-median(meanqds2(:));
    qd3(:,:,ii)=meanqds3-median(meanqds3(:));
    qd4(:,:,ii)=meanqds4-median(meanqds4(:));
end
% qd1=tiffread('G:\4Pi_paper\SUM_Cell17_642_000_080_q1.tif');
% qd2=tiffread('G:\4Pi_paper\SUM_Cell17_642_000_080_q2.tif');
% qd3=tiffread('G:\4Pi_paper\SUM_Cell17_642_000_080_q3.tif');
% qd4=tiffread('G:\4Pi_paper\SUM_Cell17_642_000_080_q4.tif');

%% find Fourier-Mellin transform
qdall=cat(4,qd1,qd2,qd3,qd4);
qdall(qdall<=0)=1e-10;
zm_all=[];
trans_all=[];
ang_all=[];
R=[];
invR=[];
for ii=1:1:4
    im1=mean(qdall(:,:,:,1),3);
    im2=mean(qdall(:,:,:,ii),3);
    
%     [im1,im2]=GetLocalizationImage(im1,im2);
    
    [zm,trans,ang] = fmmatch(im2,im1);
    [out,R(:,:,ii)] = find_affine_trans(im2, im1, [[zm zm],trans,ang]);
    zm_fin=out(1:2);
    trans_fin=out(3:4);
    ang_fin=out(5);
    
%     [imout]=affine_trans(im2,zm_fin,trans_fin,ang_fin);
    zm_all(ii,:)=zm_fin;
    trans_all(ii,:)=trans_fin;
    ang_all(ii,:)=ang_fin;
    
    [zm,trans,ang] = fmmatch(im1,im2);
    [~,invR(:,:,ii)] = find_affine_trans(im1, im2, [[zm zm],trans,ang]);

    tform=affinetform2d(invR(:,:,ii));
    [imout2]=imwarp(im2,tform,'OutputView',imref2d(size(im1)),'interp','linear');
    joinchannels('RBG',im1,imout2)

    set(handles.programStatus, 'String', ['Apply Affine Transform to Channel ' num2str(ii) ' of 4']); 
    drawnow update;
end

para.centers = centers;
para.file = [dataFolder imageName];
datestring = datestr(now,'yyyymmdd');
save([resultpath namestr 'FMTtransform_' datestring],'R','invR','para','zm_all','trans_all','ang_all');
I=find(dataFolder=='\',2,'last');
parentFolder = dataFolder(1:I);
save([parentFolder namestr 'FMTtransform_' datestring],'R','invR','para','zm_all','trans_all','ang_all');
set(handles.programStatus, 'String', {'Channel Alignment', ['Color Channel:' channel],['Input:' dataFolder imageName{1}],['Output:' parentFolder namestr 'FMTtransform' datestring '.mat']}); %show progress
drawnow update; 
