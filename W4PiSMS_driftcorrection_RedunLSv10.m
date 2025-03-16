function [xout2,yout2,zout2,shifts]=W4PiSMS_driftcorrection_RedunLSv10(xout,yout,zout,tout,frmnum,reverseflag,interpflag,str)

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

if exist(str,'file')
    load(str);
else
    [drift,driftinfo,fieldc]=driftcorrection3D_so(x,y,z,frame,p);
    shiftz=driftinfo.z.dz;
    shiftz=diff(shiftz);
    shiftx=driftinfo.xy.dx;
    shifty=driftinfo.xy.dy;
    shiftx=diff(shiftx);
    shifty=diff(shifty);
    save(str,'shiftx','shifty','shiftz');
end

if ~interpflag
    [xout2]=shiftcoords_LS(xout,shiftx,tout,frmnum,reverseflag);
    [yout2]=shiftcoords_LS(yout,shifty,tout,frmnum,reverseflag);  
    [zout2]=shiftcoords_LS(zout,shiftz,tout,frmnum,reverseflag);
    shifts=[cumsum(shiftx(:)),cumsum(shifty(:)),cumsum(shiftz(:))];
else
%     id=abs(shiftz)>40;
%     if sum(id)>0
%         sum(id)
%         shiftz(id)=0;
%     end
%     id=abs(shiftx)>30;
%     sum(id)
%     shiftx(id)=0;    
%     id=abs(shifty)>50;
%     sum(id)
%     shifty(id)=0;
    
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
    if reverseflag
        finaldrift=finaldrift-ones(length(finaldrift),1)*finaldrift(end,:);
    end
    shift=finaldrift(tout+1,1);
    xout2=xout-shift;
    shift=finaldrift(tout+1,2);
    yout2=yout-shift;
    shift=finaldrift(tout+1,3);
    zout2=zout-shift;
    shifts=finaldrift;
    
    pixelsz=128;
    drift_xyz=finaldrift;
    drift_xyz(:,1:2) = drift_xyz(:,1:2)/pixelsz;    % pixel
    drift_xyz(:,3) = drift_xyz(:,3)/1000;           % um
    drift_xyz(:,1) = drift_xyz(:,1)-mean(drift_xyz(:,1));
    drift_xyz(:,2) = drift_xyz(:,2)-mean(drift_xyz(:,2));
    drift_xyz(:,3) = drift_xyz(:,3)-mean(drift_xyz(:,3));       
    rccstr=insertBefore(str, ".mat", "_RCC");
    save(rccstr,'drift_xyz'); 
    
%     x=xout2;
%     y=yout2;
%     z=zout2;
%     p.correctxy=true;
%     p.zrange=[min(z) max(z)];
%     [drift,driftinfo,fieldc]=driftcorrection3D_so(x,y,z,frame,p);
%     shiftx=driftinfo.xy.dx;
%     shifty=driftinfo.xy.dy;
%     shiftx=diff(shiftx);
%     shifty=diff(shifty);
%     [xout2]=shiftcoords_LS(xout2,shiftx,tout,frmnum,reverseflag);
%     [yout2]=shiftcoords_LS(yout2,shifty,tout,frmnum,reverseflag);
%     shiftz=driftinfo.z.dz;
%     shiftz=diff(shiftz);
%     [zout2]=shiftcoords_LS(zout2,shiftz,tout,frmnum,reverseflag);
    
    x=xout2;
    y=yout2;
    z=zout2;
    p.correctxy=false;
    p.zrange=[min(z) max(z)];
    [drift,driftinfo,fieldc]=driftcorrection3D_so(x,y,z,frame,p);
    shiftz=driftinfo.z.dz;
    shiftz=diff(shiftz);
    [zout2]=shiftcoords_LS(zout2,shiftz,tout,frmnum,reverseflag);

%     z=zout2;
%     p.zrange=[min(z) max(z)];
%     [drift,driftinfo,fieldc]=driftcorrection3D_so(x,y,z,frame,p);
%     shiftz=driftinfo.z.dz;
%     shiftz=diff(shiftz);
%     [zout2]=shiftcoords_LS(zout2,shiftz,tout,frmnum,reverseflag);
end
