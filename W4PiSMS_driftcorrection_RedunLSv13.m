function [finaldrift]=W4PiSMS_driftcorrection_RedunLSv13(xout,yout,zout,tout,frmnum)

x=xout;
y=yout;
z=zout;
frame=tout;
p.framestart=0;
p.framestop=ceil((max(frame)/frmnum))*frmnum-1;
points=ceil((max(frame)/frmnum));
p.correctxy=true;
p.drift_timepoints=points;
p.correctz=true;
p.drift_timepointsz=points;
p.zrange=[min(z) max(z)];

[drift,driftinfo,fieldc]=driftcorrection3D_so(x,y,z,frame,p);
shiftz=driftinfo.z.dz;
shiftz=diff(shiftz);
shiftx=driftinfo.xy.dx;
shifty=driftinfo.xy.dy;
shiftx=diff(shiftx);
shifty=diff(shifty);

% cutmeth='nocut';
% pixelsz=16; % nm
% thresh=8; % nm
% if sum(shiftx.^2)<sum(shifty.^2)
%     [shiftx2,shiftz]=correct_drift_LS(single(xout),single(zout),tout,frmnum,pixelsz,thresh,cutmeth);
% else
%     [shifty2,shiftz]=correct_drift_LS(single(yout),single(zout),tout,frmnum,pixelsz,thresh,cutmeth);
% end

% [shiftx,shifty,shiftz,~,~,~,~]=W4PiSMS_driftcorrection_RedunLS3D(single(xout),single(yout),single(zout),pixelsz,thresh,cutmeth,frmnum,tout);

ntotalframe=max(tout)+1;
nbinframe=length(shiftx)+1;
indexinterp=zeros(nbinframe+2,1);
indexinterp(1)=1;
indexinterp(nbinframe+2)=ntotalframe;
indexinterp(2:nbinframe+1)=round(frmnum/2):frmnum:frmnum*nbinframe-1;
drift=cumsum(shiftx);
finaldrift(:,1)=interp1(indexinterp,[0 0 drift' drift(end,1)],1:ntotalframe,'linear')';
drift=cumsum(shifty);
finaldrift(:,2)=interp1(indexinterp,[0 0 drift' drift(end,1)],1:ntotalframe,'linear')';
drift=cumsum(shiftz);
finaldrift(:,3)=interp1(indexinterp,[0 0 drift' drift(end,1)],1:ntotalframe,'linear')';

