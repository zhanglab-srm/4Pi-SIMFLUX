function result = ProcessFilePair(imgpath_7189, varargin)
% result = ProcessFilePair(imgpath_7189, varargin)
% Process 2 image pair in ROSE, image name format: spool_7189_X2.tif,
%                                                   spool_7002_X2.tif
% sub-image duplicate, registration, detection, fitting
% results save into matfile
%
%options:
%
%
%

%% init default parameters
result = -1;
tpos = strfind(imgpath_7189, '7189');
if isempty(tpos) || len(tpos)>1
    warning(sprintf('Filename err: 7189 find: %d', len(tpos)));
    imgpath_7002 = '';
else
    imgpath_7002 = [imgpath_7189(1:tpos-1), '7002', imgpath_7189(tpos+4:end)];
end

isEM = 0;
tform_matfile = [];
gain = 100;
isFT = 1;
doIntComp = 0; 
%detection
threshold = 3;
windowWidth = 7;

%% apply input parameters
parnum = len(varargin);
if mod(parnum,2) >0
    error('number of parameter error! (%d)', parnum)
end

for m=1:2:parnum
    parname = varargin{m};
    pardata = varargin{m+1};
    switch parname
        case 'isEM'
            isEM = str2double(pardata);
        case 'isFT'
            isFT = str2double(pardata);
        case 'doIntComp'
            doIntComp = str2double(pardata);
        case 'gain'
            gain = str2double(pardata);
        case 'tform_matfile'
            tform_matfile = pardata;
        case 'imgpath_7002'
            imgpath_7002 = pardata;
        case 'threshold'
            threshold = str2double(pardata);
        case 'windowWidth'
            windowWidth = str2double(pardata);
    end
end


%% start process
fprintf('Processing Files: \n%s\n%s\n', imgpath_7189, imgpath_7002);

imgpath_prefix = imgpath_7189(1:end-4);
imgpath_prefix2 = imgpath_7002(1:end-4);
%% duplicate and registration
dupRegImagesLite(imgpath_7189, tform_matfile, [], isEM, isFT);
dupRegImagesLite(imgpath_7002, tform_matfile, [], isEM, isFT);

%% detection
SMDetection([imgpath_prefix '_outsum.tif'], threshold, windowWidth);
%% fitting
detfile = [imgpath_prefix '_outsum_Detection.mat'];
imglist = cell(1, 6);
imglist{1} = [imgpath_prefix '_out1.tif'];
imglist{2} = [imgpath_prefix '_out2.tif'];
imglist{3} = [imgpath_prefix '_out3.tif'];
imglist{4} = [imgpath_prefix2 '_out4.tif'];
imglist{5} = [imgpath_prefix2 '_out5.tif'];
imglist{6} = [imgpath_prefix2 '_out6.tif'];
[fitresult, smInfo, intCompensationData] = SMFitting2(imglist, detfile, isEM, gain, doIntComp);

%% save to file
fitfile = [imgpath_prefix '_outsum_Fitting.mat'];
save(fitfile, 'fitresult', 'smInfo', 'detfile', 'imgpath_7189', 'imgpath_7002', 'intCompensationData');

result = 0;
end