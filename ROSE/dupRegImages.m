function [imgbuf1, imgbuf2b, imgbuf3b, imgbuf4b, imgbuf5b, imgbuf6b] = ...
    dupRegImages(img_7189, img_7002, matfile, outfile_prefix)
%% image duplicate and registration
% img_7189 = 'E:\ROSE\20181015_2CCD_Registration\spool_7189.tif';
% img_7002 ='E:\ROSE\20181015_2CCD_Registration\spool_7002.tif';
% matfile = '.\RegInfo.mat';
% 
% outfile_prefix = 'E:\ROSE\20181015_2CCD_Registration\proc_out_';

%% load parameters
temp = load(matfile);
roibuf = temp.roibuf;
tform2 = temp.tform2;
tform3 = temp.tform3;
tform4 = temp.tform4;
tform5 = temp.tform5;
tform6 = temp.tform6;

%% Load Images
imgbuf_7002 = LoadTiff16bit(img_7002);
imgbuf_7189 = LoadTiff16bit(img_7189);
s1 = size(imgbuf_7002);
s2 = size(imgbuf_7189);
imglen = s1(3);

%% make sub-images
imgbuf1 = imgbuf_7189(roibuf(1,2):roibuf(1,4), roibuf(1,1):roibuf(1,3),:);
imgbuf2 = imgbuf_7189(roibuf(2,2):roibuf(2,4), roibuf(2,1):roibuf(2,3),:);
imgbuf3 = imgbuf_7189(roibuf(3,2):roibuf(3,4), roibuf(3,1):roibuf(3,3),:);
imgbuf4 = imgbuf_7002(roibuf(4,2):roibuf(4,4), roibuf(4,1):roibuf(4,3),:);
imgbuf5 = imgbuf_7002(roibuf(5,2):roibuf(5,4), roibuf(5,1):roibuf(5,3),:);
imgbuf6 = imgbuf_7002(roibuf(6,2):roibuf(6,4), roibuf(6,1):roibuf(6,3),:);

clear imgbuf_7189
clear imgbuf_7002
%% Tranform images
xdata = [1 size(imgbuf1,2)];
ydata = [1 size(imgbuf1,1)];

imgbuf2b = zeros(size(imgbuf2));
parfor m=1:imglen
    imgbuf2b(:,:,m) = imtransform(imgbuf2(:,:,m), tform2,'XData',xdata,'YData',ydata);
end

imgbuf3b = zeros(size(imgbuf3));
parfor m=1:imglen
    imgbuf3b(:,:,m) = imtransform(imgbuf3(:,:,m), tform3,'XData',xdata,'YData',ydata);
end

imgbuf4b = zeros(size(imgbuf4));
parfor m=1:imglen
    imgbuf4b(:,:,m) = imtransform(imgbuf4(:,:,m), tform4,'XData',xdata,'YData',ydata);
end

imgbuf5b = zeros(size(imgbuf5));
parfor m=1:imglen
    imgbuf5b(:,:,m) = imtransform(imgbuf5(:,:,m), tform5,'XData',xdata,'YData',ydata);
end

imgbuf6b = zeros(size(imgbuf6));
parfor m=1:imglen
    imgbuf6b(:,:,m) = imtransform(imgbuf6(:,:,m), tform6,'XData',xdata,'YData',ydata);
end

%% save images
tiffwriteStack(imgbuf1, [outfile_prefix '1.tif']);
tiffwriteStack(imgbuf2b, [outfile_prefix '2.tif']);
tiffwriteStack(imgbuf3b, [outfile_prefix '3.tif']);
tiffwriteStack(imgbuf4b, [outfile_prefix '4.tif']);
tiffwriteStack(imgbuf5b, [outfile_prefix '5.tif']);
tiffwriteStack(imgbuf6b, [outfile_prefix '6.tif']);
tiffwriteStack(imgbuf1+imgbuf2b+imgbuf3b+imgbuf4b+imgbuf5b+imgbuf6b, [outfile_prefix 'sum.tif']);
