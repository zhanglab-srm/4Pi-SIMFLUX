function [imgbuf1, imgbuf2b, imgbuf3b, imgbuf4b, imgbuf5b, imgbuf6b] = ...
    dupRegImagesLite(imgfile, matfile, outfile_prefix, isEM, isFT)
%% image duplicate and registration
% img_7189 = 'E:\ROSE\20181015_2CCD_Registration\spool_7189.tif';
% img_7002 ='E:\ROSE\20181015_2CCD_Registration\spool_7002.tif';
% matfile = '.\RegInfo.mat';
% 
% outfile_prefix = 'E:\ROSE\20181015_2CCD_Registration\proc_out_';
if nargin<2 || isempty(matfile)
    matfile = '..\PreProcess\RegInfo.mat';
end
if nargin <3 || isempty(outfile_prefix)
    outfile_prefix = [imgfile(1:end-4) '_out'];
end

if nargin <4
    warning('isEM is NOT set, using default: 1');
    isEM = 1;
end

if nargin <5
    warning('isFT is NOT set, using default: 1');
    isFT = 1;
end

is7002 = 0;
is7189 = 0;

if contains(imgfile, '7002')
    is7002 = 1;
end
if contains(imgfile, '7189')
    is7189 = 1;
end

%% check if file exist
if is7189>0 %7189: _out1, _out2, _out3, _outsum
    if exist([outfile_prefix '1.tif'], 'file') && ...
            exist([outfile_prefix '2.tif'], 'file') && ...
            exist([outfile_prefix '3.tif'], 'file') && ...
            exist([outfile_prefix 'sum.tif'], 'file')
        disp('Transformed file exist, skip image duplicate and registration ..');
        return
    end
else %7002: _out4, _out5, _out6, _outsum
    if exist([outfile_prefix '4.tif'], 'file') && ...
            exist([outfile_prefix '5.tif'], 'file') && ...
            exist([outfile_prefix '6.tif'], 'file') && ...
            exist([outfile_prefix 'sum.tif'], 'file')
        disp('Transformed file exist, skip image duplicate and registration ..');
        return
    end
end
%% load parameters
temp = load(matfile);
roibuf = temp.roibuf;
tform2 = temp.tform2;
tform3 = temp.tform3;
tform4 = temp.tform4;
tform5 = temp.tform5;
tform6 = temp.tform6;

%% Load Images
if is7002
    imgbuf_7002 = LoadTiff16bitRaw(imgfile);
    if isEM == 0%if no EM then flip image
        imgbuf_7002 = fliplr(imgbuf_7002);
    end
    s1 = size(imgbuf_7002);
    s=s1;
    imglen = s1(3);
end

if is7189
    imgbuf_7189 = LoadTiff16bitRaw(imgfile);
    if isEM == 0%if no EM then flip image
        imgbuf_7189 = fliplr(imgbuf_7189);
    end
    s2 = size(imgbuf_7189);
    s=s2;
    imglen = s2(3);
end

%% make sub-images
if isFT >0 % frame tranfer OFF, 7189 2:end, 7002 1:end-2
%     warning('Notice£ºThis function duplicate 1:end-1 for 7189data; 2:end for 7002data')
end
if is7189
    if isFT ==0 % frame tranfer OFF, 7189 2:end, 7002 1:end-2
        imgbuf1 = single(imgbuf_7189(roibuf(1,2):roibuf(1,4), roibuf(1,1):roibuf(1,3),1:end));
        imgbuf2 = single(imgbuf_7189(roibuf(2,2):roibuf(2,4), roibuf(2,1):roibuf(2,3),1:end));
        imgbuf3 = single(imgbuf_7189(roibuf(3,2):roibuf(3,4), roibuf(3,1):roibuf(3,3),1:end));
    else % frame tranfer ON, 7189 1:end-1, 7002 2:end
        imgbuf1 = single(imgbuf_7189(roibuf(1,2):roibuf(1,4), roibuf(1,1):roibuf(1,3),1:end-1));
        imgbuf2 = single(imgbuf_7189(roibuf(2,2):roibuf(2,4), roibuf(2,1):roibuf(2,3),1:end-1));
        imgbuf3 = single(imgbuf_7189(roibuf(3,2):roibuf(3,4), roibuf(3,1):roibuf(3,3),1:end-1));
    end
%     imgbuf1 = single(imgbuf_7189(roibuf(1,2):roibuf(1,4), roibuf(1,1):roibuf(1,3),:));
%     imgbuf2 = single(imgbuf_7189(roibuf(2,2):roibuf(2,4), roibuf(2,1):roibuf(2,3),:));
%     imgbuf3 = single(imgbuf_7189(roibuf(3,2):roibuf(3,4), roibuf(3,1):roibuf(3,3),:));
    clear imgbuf_7189
    s = size(imgbuf1);
    imglen = s(3);
end
if is7002
    if isFT ==0 % frame tranfer OFF, 7189 2:end, 7002 1:end-2
        imgbuf4 = single(imgbuf_7002(roibuf(4,2):roibuf(4,4), roibuf(4,1):roibuf(4,3),1:end));
        imgbuf5 = single(imgbuf_7002(roibuf(5,2):roibuf(5,4), roibuf(5,1):roibuf(5,3),1:end));
        imgbuf6 = single(imgbuf_7002(roibuf(6,2):roibuf(6,4), roibuf(6,1):roibuf(6,3),1:end));
    else % frame tranfer ON, 7189 1:end-1, 7002 2:end
        imgbuf4 = single(imgbuf_7002(roibuf(4,2):roibuf(4,4), roibuf(4,1):roibuf(4,3),2:end));
        imgbuf5 = single(imgbuf_7002(roibuf(5,2):roibuf(5,4), roibuf(5,1):roibuf(5,3),2:end));
        imgbuf6 = single(imgbuf_7002(roibuf(6,2):roibuf(6,4), roibuf(6,1):roibuf(6,3),2:end));
    end
%     imgbuf4 = single(imgbuf_7002(roibuf(4,2):roibuf(4,4), roibuf(4,1):roibuf(4,3),:));
%     imgbuf5 = single(imgbuf_7002(roibuf(5,2):roibuf(5,4), roibuf(5,1):roibuf(5,3),:));
%     imgbuf6 = single(imgbuf_7002(roibuf(6,2):roibuf(6,4), roibuf(6,1):roibuf(6,3),:));
    clear imgbuf_7002
    s = size(imgbuf4);
    imglen = s(3);
end

% if is7189
%     tiffwriteStack(imgbuf1, [outfile_prefix '1.tif']);
%     tiffwriteStack(imgbuf2, [outfile_prefix '2.tif']);
%     tiffwriteStack(imgbuf3, [outfile_prefix '3.tif']);
%     tiffwriteStack(imgbuf1+imgbuf2+imgbuf3, [outfile_prefix 'sum.tif']);
% end
% 
% if is7002
%     tiffwriteStack(imgbuf4, [outfile_prefix '4.tif']);
%     tiffwriteStack(imgbuf5, [outfile_prefix '5.tif']);
%     tiffwriteStack(imgbuf6, [outfile_prefix '6.tif']);
%     tiffwriteStack(imgbuf4+imgbuf5+imgbuf6, [outfile_prefix 'sum.tif']);
% end

%% Tranform images
xdata = [1 s(2)];
ydata = [1 s(1)];
imref = imref2d(s(1:2));
if is7189
    imgbuf2b = zeros(size(imgbuf2));
    for m=1:imglen
%         imgbuf2b(:,:,m) = imtransform(imgbuf2(:,:,m), tform2,'XData',xdata,'YData',ydata);
        imgbuf2b(:,:,m) = imwarp(imgbuf2(:,:,m), tform2,'OutputView',imref);
    end

    imgbuf3b = zeros(size(imgbuf3));
    for m=1:imglen
%         imgbuf3b(:,:,m) = imtransform(imgbuf3(:,:,m), tform3,'XData',xdata,'YData',ydata);
        imgbuf3b(:,:,m) = imwarp(imgbuf3(:,:,m), tform3,'OutputView',imref);
    end
end

if is7002
    imgbuf4b = zeros(size(imgbuf4));
    for m=1:imglen
%         imgbuf4b(:,:,m) = imtransform(imgbuf4(:,:,m), tform4,'XData',xdata,'YData',ydata);
        imgbuf4b(:,:,m) = imwarp(imgbuf4(:,:,m), tform4,'OutputView',imref);
    end

    imgbuf5b = zeros(size(imgbuf5));
    for m=1:imglen
%         imgbuf5b(:,:,m) = imtransform(imgbuf5(:,:,m), tform5,'XData',xdata,'YData',ydata);
        imgbuf5b(:,:,m) = imwarp(imgbuf5(:,:,m), tform5,'OutputView',imref);
    end

    imgbuf6b = zeros(size(imgbuf6));
    for m=1:imglen
%         imgbuf6b(:,:,m) = imtransform(imgbuf6(:,:,m), tform6,'XData',xdata,'YData',ydata);
        imgbuf6b(:,:,m) = imwarp(imgbuf6(:,:,m), tform6,'OutputView',imref);
    end
end
%% save images
if is7189
    tiffwriteStack(imgbuf1, [outfile_prefix '1.tif']);
    tiffwriteStack(imgbuf2b, [outfile_prefix '2.tif']);
    tiffwriteStack(imgbuf3b, [outfile_prefix '3.tif']);
    tiffwriteStack(imgbuf1+imgbuf2b+imgbuf3b, [outfile_prefix 'sum.tif']);
end

if is7002
    tiffwriteStack(imgbuf4b, [outfile_prefix '4.tif']);
    tiffwriteStack(imgbuf5b, [outfile_prefix '5.tif']);
    tiffwriteStack(imgbuf6b, [outfile_prefix '6.tif']);
    tiffwriteStack(imgbuf4b+imgbuf5b+imgbuf6b, [outfile_prefix 'sum.tif']);
end


