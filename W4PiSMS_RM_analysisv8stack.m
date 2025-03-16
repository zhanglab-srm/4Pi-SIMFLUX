function [xresult, yresult, tresult, zfresult, zangresult, llresult, CRLBresult, Iresult, zast_err_result, stacktot, num_images_com, imagesz, zangctrresult, bgresult, subimstot, sxtot, sytot, smInfotot]=W4PiSMS_RM_analysisv8stack(foldername,nameroot,fmfile,astfile,anglefile,llrthresh,scmos_cali_file,det_thresh,CRLB_thresh,sigma_thresh,centers,handles)
%% parameters  
fit_flag=str2double(get(handles.fit_flag,'string')); % 1: 4pi, 2: 2D, 3: 3D
if fit_flag==2
    subsz=7;   % sub region size
else
    subsz=13;
end

%% camera calibration file
% load(handles.scmos_cali_file);
% offsetim=ones(200,2304)*100;
% varim=ones(200,2304);
% gainim=ones(200,2304)*4.35;
% caliims=cat(3,offsetim,varim,gainim);

%%
xresult=[];
yresult=[];
zfresult=[];
zangresult=[];
tresult=[];
zast_err_result=[];
llresult=[];
CRLBresult=[];
Iresult=[];
stacktot=[];
zangctrresult=[];
bgresult=[];
sxtot=[];
sytot=[];
smInfotot=[];

%%
files=dir([foldername nameroot]);
% [~,idx] = sort([files.datenum]);
lastf=0;
subimstot=[];
for ff=1:numel(files)
    disp(['Reading dcimg file: ',num2str(ff)]);
%     filename=files(idx(ff)).name;
    filename=files(ff).name;
    tic
    [qd1,qd2,qd3,qd4]=W4PiSMS_readdcimg_v2([foldername filename],centers);
    toc
%     tmpim=calicrops(:,:,3,:);
%     tmpim(tmpim<1.3|tmpim>3.5)=mean(mean(mean(calicrops(:,:,3,:))));
%     calicrops(:,:,3,:)=tmpim;
%     offsetim=squeeze(calicrops(:,:,1,:));
%     varim=squeeze(calicrops(:,:,2,:));
%     gainim=squeeze(calicrops(:,:,3,:));
    if handles.PAINT
        stacknum=getstepnum_v2(filename);   % PAINT
    else
        stacknum=getstepnum(filename);   
    end
%     qd1=(qds(:,:,:,1)-repmat(offsetim(:,:,1),[1 1 size(qds,3) 1]))./repmat(gainim(:,:,1),[1 1 size(qds,3) 1]);
%     qd2=(qds(:,:,:,2)-repmat(offsetim(:,:,2),[1 1 size(qds,3) 1]))./repmat(gainim(:,:,2),[1 1 size(qds,3) 1]);
%     qd3=(qds(:,:,:,3)-repmat(offsetim(:,:,3),[1 1 size(qds,3) 1]))./repmat(gainim(:,:,3),[1 1 size(qds,3) 1]);
%     qd4=(qds(:,:,:,4)-repmat(offsetim(:,:,4),[1 1 size(qds,3) 1]))./repmat(gainim(:,:,4),[1 1 size(qds,3) 1]);
%     num_images=size(qd1,3);
%     num_images_com=num_images/6;
 %
%  qds=[];
    
    %% median filtering
%     if handles.MedFilter.Value
%         disp('Median filtering...');
%         tic
%         [qd1,qd2,qd3,qd4]=medFilter(qd1,qd2,qd3,qd4,num_images);
%         toc
%     end
    
    %% rotate and align
    disp('Channel alignment...');
    tic
    [sumim,qd1,qd2,qd3,qd4]=W4PiSMS_RotAlign_FMT_v3(qd1,qd2,qd3,qd4,fmfile);
    toc
    num_images=size(sumim,3);
    num_images_com=size(qd1,3);
    imagesz=size(qd1,1);

    %% check shift
%     [shifts]=getShifts(q1,q2,q3,q4);
%     [q1,q2,q3,q4]=checkShifts(q1,q2,q3,q4,shifts);
    
    %%
%     sumim1=q1+q2+q3+q4;
%     sumim1(sumim1<=1e-6)=1e-6;
%     q1=[];q2=[];q3=[];q4=[];

    %% combine 6 frames
    id=1:6:num_images;
    img1=sumim(:,:,id);
    id=2:6:num_images;
    img2=sumim(:,:,id);
    id=3:6:num_images;
    img3=sumim(:,:,id);
    id=4:6:num_images;
    img4=sumim(:,:,id);
    id=5:6:num_images;
    img5=sumim(:,:,id);
    id=6:6:num_images;
    img6=sumim(:,:,id);
    sumim=[];
    sumim1=img1+img2+img3+img4+img5+img6;
   
    %% save to tif images
    if ff==1
        filestr=[foldername filename(1:end-5) 'tif'];
        if ~exist(filestr,'file')
            tiffwrite(sumim1,filestr,[1,500]);
        end
    end
    
    %%
    pick_flag=1; % use range
    tic
    [sub_regions,tlz,locmaxc]=W4PiSMS_sumim_seg(sumim1,det_thresh,subsz,pick_flag,fit_flag);
    toc
    disp(['A total of ' num2str(size(sub_regions,3)) ' subregions were detected. Start sCMOS_sigmaxy fitting']);

    %% show detection result
    f=100;
    id=locmaxc(:,3)==f-1;
    V=locmaxc(id,:)+1;
    close all
    figure;imshow(sumim1(:,:,f),[0 1000],'InitialMagnification',250);hold on;plot(V(:,1),V(:,2),'bo');pause(1)
    sumim1=[];
        
    %%
    tic
%     Nt=size(sub_regions,3);
%     Nf=ceil(Nt/1000000);
%     P=[];
%     CRLB=[];
%     LL=[];
%     for k=1:Nf
%         st=(k-1)*1000000+1;
%         et=min(k*1000000,Nt);
%         P1=[];
%         CRLB1=[];
%         LL1=[];
%         if fit_flag==2
%             [P1,CRLB1,LL1]=mleFit_LM(single(sub_regions(:,:,st:et)),2,50,1.4,0,0,0);
%             P1(:,6)=P1(:,5);
%         else
%             [P1,CRLB1,LL1]=mleFit_LM(single(sub_regions(:,:,st:et)),4,50,1.4,0,0,0);
%         end         
%         P=cat(1,P,P1);
%         CRLB=cat(1,CRLB,CRLB1);
%         LL=cat(1,LL,LL1);
%     end
    [P,CRLB,LL]=mleFit_LM(single(sub_regions),4,50,1.4,0,0,0);
    toc
    disp('Fitting finished...Start Z_ast estimation');
    
    xco=P(:,2);
    yco=P(:,1);
    I=P(:,3);
    llr=-2*LL;
    sigmax=P(:,5);
    sigmay=P(:,6);
    bg=P(:,4);

    %%
    sum_ims=single(zeros(subsz,subsz,6,size(sub_regions,3)));
    sub_regions=[]; 
    sum_ims(:,:,1,:)=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(permute(img1,[1 2 3])));
    sum_ims(:,:,2,:)=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(permute(img2,[1 2 3])));
    sum_ims(:,:,3,:)=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(permute(img3,[1 2 3])));
    sum_ims(:,:,4,:)=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(permute(img4,[1 2 3])));
    sum_ims(:,:,5,:)=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(permute(img5,[1 2 3])));
    sum_ims(:,:,6,:)=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(permute(img6,[1 2 3])));
    img1=[];img2=[];img3=[];img4=[];img5=[];img6=[];

    %%
    [~,smInfo,smInfomask]=ROSE_SMFitting(sum_ims,tlz);
    sum_ims=[];

    %% filter
    rr=(subsz-1)/2;
    maskxy=abs(xco-rr)<3&abs(yco-rr)<3;
    maskothers=CRLB(:,1)>0&CRLB(:,2)>0&CRLB(:,1)<1&CRLB(:,2)<1&I>0&bg>0&llr>0&sigmax>0&sigmax<5&sigmay>0&sigmay<5;   
    mask=maskxy&maskothers&smInfomask;
    
    xf=xco(mask);
    yf=yco(mask);
    llr_f=llr(mask);
    CRLB(:,1:2)=sqrt(CRLB(:,1:2));
    CRLB_f=CRLB(mask,:);
    I_f=I(mask);
    tlz_f=tlz(mask,:);
    sigmaxf=sigmax(mask);
    sigmayf=sigmay(mask);
    locmaxc_f=locmaxc(mask,:);
    stacknum_f=repmat(stacknum(1),size(xf));
    bg_f=bg(mask);

   %% xy-phase
    smInfo=smInfo(mask,:);
    pixelsize=128;

    id1=(80<=smInfo(:,1)&smInfo(:,1)<=120&smInfo(:,3)<=1000);
    xx1=smInfo(id1,1)*pixelsize;
    yy1=smInfo(id1,5);
    f1=figure(2);
    scatter(xx1,yy1,'filled');

    id2=(80<=smInfo(:,2)&smInfo(:,2)<=120&smInfo(:,3)<=1000);
    xx2=smInfo(id2,2)*pixelsize;
    yy2=smInfo(id2,6);
    f2=figure(3);
    scatter(xx2,yy2,'filled');

%     saveas(f1, [detFileInfo(1:end-22) '_phase1.png']);
%     saveas(f2, [detFileInfo(1:end-22) '_phase2.png']);

    %% fit phase plane
    disp('---------- Fitting Phase Plane ----------');
    tx1 = smInfo(:,2);
    ty1 = smInfo(:,1);
%     tx1 = double(xf+tlz_f(:,2));
%     ty1 = double(yf+tlz_f(:,1));
    p1 = smInfo(:,6);
    p2 = smInfo(:,5);
    
    [lp1, lp2, resultp1, resultp2, resultx, resulty, resultdx, resultdy, A] ...
        = CalculatePhaseAndFineXY_faster_6(tx1, ty1, p1, p2);
    
    %function for phase-position translate
%     phase2pos_func = @(p1, p2)[p1(:), p2(:)]/A;

    smInfo(:,11)=resultx;
    smInfo(:,10)=resulty;
    smInfo(:,13)=resultdy;
    smInfo(:,12)=resultdx;
    
    %% z ast initial guess
    ast_config=astfile;
    zf=[];
    zerr=[];
    
    %% z estimation by astigmatism by look up table
    if ff==1
        tmp=load(ast_config);
        zdata=(0:1200)';
        Yf=[];
        params=tmp.estx;
        w=params(1);
        c=params(2);
        d=params(3);
        A=params(4);
        B=params(5);
        Yf(:,1) = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
        params=tmp.esty;
        w=params(1);
        c=params(2);
        d=params(3);
        A=params(4);
        B=params(5);
        Yf(:,2) = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
        dif=abs(Yf(:,1)-Yf(:,2));
        id=find(dif==min(dif(:)));
        
        % center
        zdata=id-800:id+800;
        Yf=[];
        params=tmp.estx;
        w=params(1);
        c=params(2);
        d=params(3);
        A=params(4);
        B=params(5);
        Yf(:,1) = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
        params=tmp.esty;
        w=params(1);
        c=params(2);
        d=params(3);
        A=params(4);
        B=params(5);
        Yf(:,2) = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
        M=length(zdata);
    end
    
    if fit_flag==2
        zerr=zeros(length(xf),1);
        zf=zeros(length(xf),1);
    else
        N=length(sigmaxf);
        zf=zeros(N,1);
        zerr=zeros(N,1);
        tic
        parfor ii=1:N
            w0=[sigmaxf(ii),sigmayf(ii)];
            dist=Yf.^0.5-(ones(M,1)*w0).^0.5;
            dist=sum(dist.^2,2);
            id=dist==min(dist);
            zerr(ii,1)=min(dist);
            zf(ii,1)=mean(zdata(id));
        end
        toc
    end    
    disp('Z_ast fitting finished...Start phase estimation');
    
    %% determine the phase
    if isempty(locmaxc_f)
        continue
    end
    
    if fit_flag==1
        tic
        [z_ang,ang_ctr]=W4PiSMS_RM_zangv4(anglefile,qd1,qd2,qd3,qd4,locmaxc_f,subsz,xf,yf,fmfile,tlz_f,bg_f,handles.channel);
        toc
    else
        z_ang=zeros(length(locmaxc_f),1);
        ang_ctr=ones(length(locmaxc_f),1);
    end

    qd1=[];qd2=[];qd3=[];qd4=[];
    disp('Phase estimation finished...Start next file...');
    
    %% collect data
    xest=xf+tlz_f(:,2);
    yest=yf+tlz_f(:,1);
    sxtot=cat(1,sxtot,sigmaxf(:));
    sytot=cat(1,sytot,sigmayf(:));
    xresult=cat(1,xresult,xest);
    yresult=cat(1,yresult,yest);
    bgresult=cat(1,bgresult,bg_f);
    zfresult=cat(1,zfresult,zf(:));
    zangresult=cat(1,zangresult,z_ang(:));
    zangctrresult=cat(1,zangctrresult,ang_ctr(:));
    tresult=cat(1,tresult,tlz_f(:,3)+(ff-1)*num_images_com);
    zast_err_result=cat(1,zast_err_result,zerr(:));
    llresult=cat(1,llresult,llr_f);
    CRLBresult=cat(1,CRLBresult,CRLB_f);
    Iresult=cat(1,Iresult,I_f);
    lastf=lastf+tlz(end,3);
    stacktot=cat(1,stacktot,stacknum_f); 
    smInfotot=cat(1,smInfotot,smInfo);
end

%% filtering 
mask=zangctrresult>0&zangctrresult<sqrt(2)&Iresult>0&llresult>0;
sxtot=sxtot(mask);
sytot=sytot(mask);
xresult=xresult(mask);
yresult=yresult(mask);
bgresult=bgresult(mask);
llresult=llresult(mask);
CRLBresult=CRLBresult(mask,:);
Iresult=Iresult(mask);
tresult=tresult(mask);
zfresult=single(zfresult(mask));
zangresult=single(zangresult(mask));
zangctrresult=single(zangctrresult(mask));
zast_err_result=single(zast_err_result(mask));
stacktot=single(stacktot(mask));
smInfotot=smInfotot(mask,:);