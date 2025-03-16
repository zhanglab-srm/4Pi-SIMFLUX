function [sumim,qd1,qd2,qd3,qd4]=W4PiSMS_RotAlign_FMT_v3(qd1,qd2,qd3,qd4,fmtname)

tmp=load(fmtname);
invR=tmp.invR;
tform2=affinetform2d(invR(:,:,2));
tform3=affinetform2d(invR(:,:,3));
tform4=affinetform2d(invR(:,:,4));

% q1=single(qd1);
q2=single(zeros(size(qd2)));
q3=single(zeros(size(qd3)));
q4=single(zeros(size(qd4)));
frames=size(qd1,3);
sz=size(qd1(:,:,1));

interp='linear';
parfor ii=1:frames
    q2(:,:,ii)=single(imwarp(qd2(:,:,ii),tform2,'OutputView',imref2d(sz),'interp',interp));
    q3(:,:,ii)=single(imwarp(qd3(:,:,ii),tform3,'OutputView',imref2d(sz),'interp',interp));
    q4(:,:,ii)=single(imwarp(qd4(:,:,ii),tform4,'OutputView',imref2d(sz),'interp',interp));
    if mod(ii,1000)==0
        disp([num2str(ii) ' out of ' num2str(size(qd1,3)) ' is done...']);
    end
end

sumim=qd1+q2+q3+q4;
q2=[];q3=[];q4=[];
sumim(sumim<=1e-6)=1e-6;

frames1=frames/6;
qtemp=qd1;
qd1=single(zeros(sz(1),sz(2),frames1));
for ii=1:frames1
    id=(1:6)+(ii-1)*6;
    qd1(:,:,ii)=sum(qtemp(:,:,id),3);
end

qtemp=qd2;
qd2=single(zeros(sz(1),sz(2),frames1));
for ii=1:frames1
    id=(1:6)+(ii-1)*6;
    qd2(:,:,ii)=sum(qtemp(:,:,id),3);
end

qtemp=qd3;
qd3=single(zeros(sz(1),sz(2),frames1));
for ii=1:frames1
    id=(1:6)+(ii-1)*6;
    qd3(:,:,ii)=sum(qtemp(:,:,id),3);
end

qtemp=qd4;
qd4=single(zeros(sz(1),sz(2),frames1));
for ii=1:frames1
    id=(1:6)+(ii-1)*6;
    qd4(:,:,ii)=sum(qtemp(:,:,id),3);
end
