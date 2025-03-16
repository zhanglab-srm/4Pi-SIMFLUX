function imgbuf = LoadTiff16bitRaw(filename, framerange)
if nargin<2
    framerange = [1 -1];
end
if ~exist(filename,'file')
    error(['No file found! filename: ' filename]);
else
    imgbuf = TIFFloadFrame_16bit_CPU(filename,framerange);
end