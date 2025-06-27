function [P,CRLB,LL]=GPUmleFit_LM_4Pi_old(imstack,sigma)

q1=imstack(:,:,:,1);
q2=imstack(:,:,:,2);
q3=imstack(:,:,:,3);
q4=imstack(:,:,:,4);
sub_regions=q1+q2+q3+q4;

[P,CRLB,LL]=mleFit_LM(sub_regions,4,50,[1.3,1.3],0,0,0); 

for ii=1:4
    ycenter(:,ii)=P(:,2);
    xcenter(:,ii)=P(:,1);
end
bg_f=P(:,4);

subims=cat(4,q1,q2,q3,q4);

% [rm1,rm2]=iPALMast_findmom_givenC(subims,xcenter(:,1),ycenter(:,1));

sz=size(subims,1);
Model=single(zeros(sz,sz,size(subims,3),4));
for ii=1:4
    yf=ycenter(:,ii);
    xf=xcenter(:,ii);
    ROI=finitegausspsf(sz,sigma,1,0,[xf,yf]);
    Model(:,:,:,ii)=single(ROI);
end

int_est=zeros(size(subims,3),size(subims,4));
for ii=1:1:size(subims,3)
    for jj=1:1:size(subims,4)
        model=Model(:,:,ii,jj);
        subim=subims(:,:,ii,jj);
        bg=bg_f(ii);
        nom=(subim-bg./4).*model;
        denom=model.^2;
        int_est(ii,jj)=sum(nom(:))./sum(denom(:));                        
    end
end

rm1=(int_est(:,1)-int_est(:,3))./(int_est(:,1)+int_est(:,3));
rm2=(int_est(:,2)-int_est(:,4))./(int_est(:,4)+int_est(:,2));

ang2=MyCalPahse_new(rm1,rm2,0,pi/2);
P(:,6)=wrapTo2Pi(ang2(:,2));