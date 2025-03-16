function [qd1,qd2,qd3,qd4]=W4PiSMS_readdcimg_v2(filename,center)

files=dir(filename);
if numel(files)>1
    error('Multiple files with the same assigned name are detected');
elseif numel(files)==0
    error('No file detected');
end

obj=dcimgReaderMatlab(filename);
% obj.metadata.num_frames=3000;
qd1=single(zeros(obj.metadata.num_rows,obj.metadata.num_rows,obj.metadata.num_frames));
qd2=single(zeros(obj.metadata.num_rows,obj.metadata.num_rows,obj.metadata.num_frames));
qd3=single(zeros(obj.metadata.num_rows,obj.metadata.num_rows,obj.metadata.num_frames));
qd4=single(zeros(obj.metadata.num_rows,obj.metadata.num_rows,obj.metadata.num_frames));
totalframes=obj.metadata.num_frames;

flipsigns=[0 0 1 1];
for ii=1:1:totalframes
    if mod(ii,1000)==0
        disp([num2str(ii) ' out of ' num2str(totalframes) ' is done...']);
    end
    ims=getSpecificFrames(obj, ii)';
    ims=(ims-100)/4.35;
    qd=single(W4PiSMS_scmos_makeqds(ims,center,flipsigns));  
    qd1(:,:,ii)=qd(:,:,:,1);
    qd2(:,:,ii)=qd(:,:,:,2);
    qd3(:,:,ii)=qd(:,:,:,3);
    qd4(:,:,ii)=qd(:,:,:,4);
end

clear mex