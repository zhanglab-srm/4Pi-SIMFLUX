function save2DHF5(procresult)
filename = procresult.filename;
s = procresult.s;
tdata = [];
tdata.x = procresult.x-1;
tdata.y = procresult.y-1;
tdata.frame = procresult.frame-1;
tdata.photon = procresult.int;
save_hdf5([filename(1:end-4) '.hdf5'],tdata);
fp = fopen([filename(1:end-4) '.yaml'],'w');
fprintf(fp,'Byte Order: ''>''\nData Type: uint16\nFile: %s\nFrames: %d\nHeight: %d\nWidth: %d\n---\nBox Size: 7\nGenerated by: Picasso Localize\nMin. Net Gradient: 5000\nROI: null\n', ...
filename, s(3),s(1),s(2));
fclose(fp);