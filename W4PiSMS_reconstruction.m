function W4PiSMS_reconstruction(handles)

mainFolder = get(handles.pathMainfolder,'string');
index_selected = get(handles.pathSubfolder,'value');
filelist = get(handles.pathSubfolder,'string');
file_selected = filelist(index_selected);

photons=str2double(get(handles.photons,'string'));
llrthresh=str2double(get(handles.llrthresh2,'string'));
mephiBin = str2double(get(handles.mephiBin,'string'));
num_images = str2double(get(handles.numImages,'string'));
crlb_thresh = str2double(get(handles.crlb,'string'));

angctr=[0.4 1.0];
st_frame=3000*0;      %startframe
frmthresh=3000*300;   %stopframe
stopval=str2double(get(handles.stopValue,'string'));     %stop value for Dmap
zthresh=[0 1200];
zasterr=[0 0.05];

%%
freqld=load('frequency_simulate4Pi_594_690.mat');
for ff=1:numel(file_selected)

    currentfile=file_selected{ff};
    I=find(currentfile=='\',1,'last');
    currentfolder=currentfile(1:I);
    fileName=currentfile(I+1:end);
    if contains(fileName,'_642_')
        channel = '642';
        reverseflag=0;
        freqstr='f_690';
        eval(['freq=freqld.' freqstr ';']);
    elseif contains(fileName,'_561_')
        channel = '561';
        reverseflag=0;
        freqstr='f_594';
         eval(['freq=freqld.' freqstr ';']);
    else
        channel = 'not 561 or 642'
    end
    savename=fileName(1:end-27);%name format foldername_642_tempresult_yymmdd.mat
    datestring = datestr(now,'yyyymmdd');
    set(handles.programStatus,'string',{'Current Analized File:',currentfile, ['Current Channel:' channel]})
    drawnow update
    
    tmpld=load(file_selected{ff});
    imagesz=tmpld.para.imagesz;
    pixelsz=tmpld.para.pixelsz;
%     imagesz=256;pixelsz=128;num_images=3000;
    if num_images==0
        num_images=tmpld.para.num_images;
    end
    
    tmpld.para.CRLB_thresh = crlb_thresh;
    tmpld.para.photon_thresh = photons;
    tmpld.para.llrthresh = llrthresh;
    tmpld.para.unwrap_thresh = stopval;

    %% Drift correction by beads
    if handles.beadsDriftCorr.Value
        Vin=[];
        Vin(:,1)=tmpld.xresult;
        Vin(:,2)=tmpld.yresult;
        Vin(:,3)=tmpld.tresult+1;
        Vin(:,5)=tmpld.Iresult;
        Vin(:,6)=tmpld.zangresult;
        Vin(:,7)=tmpld.zfresult;
        [Vout, drift, flag]=Beads_correction(Vin,reverseflag);
        if flag
            tmpld.xresult=Vout(:,1);
            tmpld.yresult=Vout(:,2);
            tmpld.zangresult=wrapToPi(Vout(:,3));
        end
    end

    %% reject beads
     Vin=[];
     Vin(:,1)=tmpld.xresult;
     Vin(:,2)=tmpld.yresult;
     Vin(:,3)=tmpld.tresult+1;
     maskbead=reject_beads(single(Vin),num_images,imagesz);

    %% rejection
    mtc=(tmpld.sxtot.^3./tmpld.sytot-tmpld.sytot.^3./tmpld.sxtot)./40*2*pi;
    id=abs(mtc)<(3*pi);
    mtc=mtc/std(mtc(id,:))*pi/2;
%     mtc=(tmpld.sytot.^3./tmpld.sxtot-tmpld.sxtot.^3./tmpld.sytot)./40*2*pi;
    % mtc mask
    maskmtc=mtc<pi&mtc>-pi &(tmpld.sxtot+tmpld.sytot)<6;% & tmpld.sxtot<2.5;
%     maskmtc=mtc<median(mtc)+1*std(mtc) & mtc>median(mtc)-1*std(mtc) & (tmpld.sxtot+tmpld.sytot)<6;
    % ll threshold
    maskll=tmpld.llresult<llrthresh;
    % crlb
    maskcrlb=tmpld.CRLBresult(:,1)<crlb_thresh & tmpld.CRLBresult(:,2)<crlb_thresh;% & tmpld.CRLBresult(:,1)>0.045 & tmpld.CRLBresult(:,2)>0.045;
    % frame number threshold
    maskt=tmpld.tresult>=st_frame & tmpld.tresult<frmthresh;
%     maskt=ones(length(mtc),1);
%     id=tmpld.tresult>=204000&tmpld.tresult<206000;
%     maskt(id)=0;
    % contrast threshold
    maskzctr=tmpld.zangctrresult>angctr(1)&tmpld.zangctrresult<angctr(2);
    % background threshold
    maskbg=tmpld.bgresult<(median(tmpld.bgresult)+2*std(tmpld.bgresult))&tmpld.bgresult>(median(tmpld.bgresult)-2*std(tmpld.bgresult));
    % astigmatism fitting threshold
    maskzast_err=tmpld.zast_err_result<zasterr(2)&tmpld.zast_err_result>zasterr(1);
    % z range threshold
    maskzf=tmpld.zfresult<zthresh(2)&tmpld.zfresult>zthresh(1);
    % # photons range threshold
    maskI=tmpld.Iresult>photons&tmpld.Iresult<100000;
    % ROI mask
    maskxy=tmpld.xresult>2000/128&tmpld.xresult<22000/128&tmpld.yresult>2000/128&tmpld.yresult<18000/128;
    % mask ROI
%     str='E:\Two-color\20231227\Cell02_1644-1674.roi';
%     roi=ReadImageJROI(str);
%     Rc=roi.mnCoordinates;
%     Rc(end+1,:)=Rc(1,:);
%     I=zeros(3276,3276);
%     BW=roipoly(I,Rc(:,1),Rc(:,2));
%     V=[tmpld.xresult*12.8,tmpld.yresult*12.8];
%     V=round(V);
%     L=length(V);
%     maskroi=zeros(L,1);
%     for jj=1:L
%         maskroi(jj)=BW(V(jj,2),V(jj,1));
%     end
    % apply
    mask=maskI&maskll&maskzctr&maskmtc&maskcrlb&maskt&maskbead;%&maskxy;%&maskroi;%&maskbead
    %mask=maskI&maskll&maskzctr&maskmtc&maskcrlb&maskt; % retain beads
    
    smInfo=tmpld.smInfotot;
    smInfomask1=smInfo(:,7)>0.4 & smInfo(:,7)<1.2 & smInfo(:,8)>0.4 & smInfo(:,8)<1.2;
    smInfomask2=smInfo(:,1)>0 & smInfo(:,2)>0 & smInfo(:,1)<imagesz & smInfo(:,2)<imagesz & abs(smInfo(:,12))<0.5 & abs(smInfo(:,13))<0.5;
    N=length(smInfomask1);
    smInfomask3=zeros(length(smInfomask1),1);
    difI=(sum(smInfo(:,14:16),2)-sum(smInfo(:,17:19),2))./sum(smInfo(:,14:19),2);
    difI=difI(difI>-0.5&difI<0.5);
    figure;histogram(difI);
    medianI=median(difI)
%     medianI=-0.05;
    for i=1:N
        intlist=smInfo(i,14:19);
        difi=(sum(intlist(1:3))-sum(intlist(4:6)))/sum(intlist);
        if abs(difi-medianI)<=0.3 && intlist(1)>0 && intlist(6)>0 && all(intlist>0)
            smInfomask3(i)=1;
        end
    end
    smInfomask=smInfomask1&smInfomask2&smInfomask3;
    mask=mask&smInfomask;

    PR=(1:length(mask))';
    PR=PR(mask);
 
    rej=[];
    nl=length(mask);
    rej(1,1)=sum(maskmtc)/nl;
    rej(2,1)=sum(maskll)/nl;
    rej(3,1)=sum(maskcrlb)/nl;
    rej(4,1)=sum(maskzctr)/nl;
    rej(5,1)=sum(maskbg)/nl;
    rej(6,1)=sum(maskI)/nl;
    rej(7,1)=sum(smInfomask)/nl;
    rej(8,1)=sum(mask)/nl;

    astfile=tmpld.astfile;
    anglefile=tmpld.anglefile;
    Iresult=tmpld.Iresult(mask);
    llresult=tmpld.llresult(mask);
    stacktot=tmpld.stacktot(mask);
    tresult=tmpld.tresult(mask)-st_frame;
%     xresult=tmpld.xresult(mask);    
%     yresult=tmpld.yresult(mask);    
    zangresult=tmpld.zangresult(mask);
    zast_err_result=tmpld.zast_err_result(mask);
    zfresult=tmpld.zfresult(mask);
    mtcresult=mtc(mask);
    sxresult=tmpld.sxtot(mask);
    syresult=tmpld.sytot(mask);
    zangctrresult=tmpld.zangctrresult(mask);
    crlbresult=tmpld.CRLBresult(mask,:);
    bgresult=tmpld.bgresult(mask,:);

    xresult=smInfo(mask,11);        %% SIMFLUX
    yresult=smInfo(mask,10);        %% SIMFLUX
    smInfonew=smInfo(mask,:);
    simid=(1:size(smInfonew(:,1)))';

    %%
    coords=[];
    coords(:,1)=xresult;
    coords(:,2)=yresult;
    im=cHistRecon(imagesz*3,imagesz*3,single(coords(:,2))*3,single(coords(:,1))*3,0);
    gaussim=gaussf(im,[1 1]);
    dipshow(gaussim); pause(eps)

    %% Phase tip/tilt correction
    currcrlb=sqrt((crlbresult(:,1).^2+crlbresult(:,2).^2)/2);
    str=[currentfolder,savename,'_pmap.mat'];
    if exist(str,'file')
        load(str);
        lp0=lp;
    else
        lp0=0;
    end
    [currzang_cor,lp] = phaseCorrection_v3(smInfo(mask,11),smInfo(mask,10),zangresult,zfresult,Iresult,currcrlb,tresult,num_images,lp0);
    save(str,'lp');
    zangresult=currzang_cor;
%     mtcresult=mtc_cor;

    %% Phase drift correction
if handles.phaseDriftCorr.Value == 1
    str=[currentfolder,savename,'_phase.mat'];
    if exist(str,'file')
        load(str);
    else
        [zangresult,mtcresult]=PhaseDriftCorrection(zangresult,mtcresult,tresult,stacktot,num_images,reverseflag);
        save(str,'zangresult','mtcresult');
    end
end

    %% dmap test
    segnum=ceil((max(tresult))/num_images);
    zest=[];
    zerr=[];
    zmask=[];
    Dmap0=[];
    Dmap=[];
    for ii=1:1:segnum
        close all
        st=(ii-1)*num_images;
        if ii==segnum
            ed=max(tresult);
        else
            ed=(ii)*num_images-1;
        end
        maskt=tresult>=st&tresult<=ed;
        mtc_seg=mtcresult(maskt);
        zang_seg=zangresult(maskt);
        [dmap]=build_dmap(mtc_seg,zang_seg,256,5);
        Dmap(:,:,ii)=dmap/max(dmap(:));
    end

    if ~handles.PAINT
        unistack=unique(stacktot);
        N=numel(unistack);
        if N>1
            OD=[];
            for j=1:N
                od=j:N:segnum;
                OD=cat(1,OD,od');
            end
            for ii=1:1:segnum
                Dmap0(:,:,ii)=Dmap(:,:,OD(ii));
            end
            Dmap=Dmap0;
        end
    end

    if handles.phaseDriftCorr.Value ==1
        str=([currentfolder savename '_dmap_DC.tif']);
    else
        str=([currentfolder savename '_dmap.tif']);
    end
    tiffwrite(Dmap*10000,str);
    programStatus_string = get(handles.programStatus,'string');
    programStatus_string = [programStatus_string;['Save Dmap Image to:' str]];
    set(handles.programStatus,'string',programStatus_string)
%     continue

    %% phase estimate for every optical section separately
    unistack=unique(stacktot);
    xout=[];
    yout=[];
    zout=[];
    shifts=[];
    z_err_out=[];
    zf_out=[];
    t_out=[];
    ll_out=[];
    I_out=[];
    bg_out=[];
    pr_out=[];
    Z_ast_err_out=[];
    crlb_out=[];
    zangctrl_out=[];
    close all
    for ss=1:numel(unistack)
        close all
        currst=unistack(ss);
        maskst=(stacktot==currst);
        currt=tresult(maskst);
        currI=Iresult(maskst);
        currll=llresult(maskst);
        currx=xresult(maskst);
        curry=yresult(maskst);
        currzang=zangresult(maskst);
        currzangctrl=zangctrresult(maskst);
        currcrlb=sqrt((crlbresult(maskst,1).^2+crlbresult(maskst,2).^2)/2);
        %tuotuo
        currrealcrlb=crlbresult(maskst,:);
        %tuotuo
        currzasterr=zast_err_result(maskst);
        currzfresult=zfresult(maskst);
        currmtc=mtcresult(maskst);
        currbg=bgresult(maskst);
        currpr=PR(maskst);
        currid=simid(maskst);

        if handles.PAINT
            currt=currt-min(currt(:));  % PAINT
        else
            if numel(unistack)>1
                currt=currt-(ss-1)*num_images;
                cycle=floor(currt/(num_images*numel(unistack)));
                remn=rem(currt,num_images*numel(unistack));
                currt=cycle*num_images+remn;
            end
        end      

        % phase correction
%         [currzang_cor]=phaseCorrection(currx,curry,currzang,currzfresult,currI,currcrlb,currt,num_images);
%         currzang=currzang_cor;

        centermtc=[];
        stopval1=stopval;
        [currzresult,z_err,mephimask]=Mephi_z_4PiSMS(currzang,currmtc,currt,num_images,mephiBin,stopval1,centermtc,freq);
        maskzerr=z_err<60;
        maskall=maskzerr&(mephimask>0)&abs(currzresult)<700;
        currzresult=-currzresult(maskall>0);
        currx=currx(maskall>0);
        curry=curry(maskall>0);
        currt=currt(maskall>0);

        z_err=z_err(maskall>0);
        currzfresult=currzfresult(maskall>0);
        currll=currll(maskall>0);
        currI=currI(maskall>0);
        currzangctrl=currzangctrl(maskall>0);
        currcrlb=currcrlb(maskall>0);
        %tuotuo
        currrealcrlb=currrealcrlb(maskall>0,:);
        crlbxyz=currrealcrlb(:,[1,2,5]);
        %tuotuo
        currzasterr=currzasterr(maskall>0);
        currbg=currbg(maskall>0);
        currpr=currpr(maskall>0);
        currid=currid(maskall>0);
        
        % astigmatic 3D 
        currzfresult=-currzfresult;
        difz=currzfresult-currzresult;
        currzfresult=currzfresult-mean(difz);

        driftstr=[currentfolder,savename,'_',num2str(ss),'_drift.mat'];
        interpflag=handles.Interpolation.Value;
%         if length(unistack)>1
%             interpflag=0;
%         end
        %tuotuo
        dc_algo=handles.DriftCorrectionAlgo.Value;     % 0:RCC 1:DME
        if dc_algo==1
            [xout{ss},yout{ss},zout{ss},shifts{ss}]=W4PiSMS_driftcorrection_DME(currx,curry,currzresult,currt,num_images,pixelsz,crlbxyz,driftstr);
        else
            [xout{ss},yout{ss},zout{ss},shifts{ss}]=W4PiSMS_driftcorrection_RedunLSv10(currx.*pixelsz,curry.*pixelsz,currzresult,currt,num_images,reverseflag,interpflag,driftstr);
        end
        %tuotuo
        %[xout{ss},yout{ss},zout{ss},shifts{ss}]=tuotuotest(currx.*pixelsz,curry.*pixelsz,currzresult,currt,pixelsz,[currentfolder savename '_tuotuo' num2str(ss) '.mat']);

        %tuotuo pass intermediate results to python
        %svcurrz=currzresult;
        %save([currentfolder savename '_dctmp' num2str(ss)],'currx','curry','svcurrz','currt','num_images','currrealcrlb','imagesz','pixelsz');
        %tuotuo

        % No drift correction
        % xout{ss}=currx.*pixelsz; yout{ss}=curry.*pixelsz;zout{ss}=currzresult;shifts{ss}=0;
       
        z_err_out{ss}=z_err;
        zf_out{ss}=currzfresult;
        ll_out{ss}=currll;
        t_out{ss}=currt;
        I_out{ss}=currI;
        zangctrl_out{ss}=currzangctrl;
        crlb_out{ss}=currcrlb;
        Z_ast_err_out{ss}=currzasterr;
        bg_out{ss}=currbg;
        pr_out{ss}=currpr;
        simid_out{ss}=currid;
    end

    save([currentfolder savename '_DCresult'],'shifts','xout','yout','zout','z_err_out','zf_out','ll_out','t_out','pr_out');
    programStatus_string =[programStatus_string;['Save DCresult to:' currentfolder savename '_DCresult']];
    set(handles.programStatus,'string',programStatus_string)
    drawnow update

    %% align stack
    if numel(xout)==1
        shiftx_st=0;
        shifty_st=0;
        shiftz_st=0;
    else
        disp('Aligning Z stacks...');
        pixelsz2=25;
        errorthresh=15;
        cutmeth='nocut';
        if handles.PAINT
            iniguess=[0 0 500];
        else
            iniguess=[0 0 500];
        end        
        maskflag=0;
        stksort=[];
        [shiftx_st,shifty_st,shiftz_st,rankf]=W4PiSMS_stack_RedunLSv9z(stksort,xout,yout,zout,pixelsz2,errorthresh,cutmeth,iniguess,maskflag,[]);
        disp('shiftz_st=' + string(shiftz_st));
        save([currentfolder savename '_alignStacks.mat'],'shiftx_st','shifty_st','shiftz_st','-mat');
    end
    [xcof]=shiftcoords_stack(xout,shiftx_st);
    [ycof]=shiftcoords_stack(yout,shifty_st);
    [zcof]=shiftcoords_stack(zout,shiftz_st);

    %% 2nd fileter result again
    zflow=-700;
    zfhigh=700;
    zerrcut=100;
    xoutf=[];
    youtf=[];
    zoutf=[];
    toutf=[];
    lloutf=[];
    Ioutf=[];
    zconf=[];
    crlbf=[];
    zasterf=[];
    bgf=[];
    prf=[];
    zerrf=[];
    simidf=[];

    for jj=1:1:numel(xcof)
        maskzerr=z_err_out{jj}<zerrcut;
        maskzf=zout{jj}<zfhigh&zout{jj}>zflow;
        mask2=maskzerr&maskzf;
        zerrf{jj}=z_err_out{jj}(mask2);
        xoutf{jj}=xcof{jj}(mask2);
        youtf{jj}=ycof{jj}(mask2);
        zoutf{jj}=zcof{jj}(mask2);
        toutf{jj}=t_out{jj}(mask2);
        Ioutf{jj}=I_out{jj}(mask2);
        lloutf{jj}=ll_out{jj}(mask2);
        zconf{jj}=zangctrl_out{jj}(mask2);
        crlbf{jj}=crlb_out{jj}(mask2);
        zasterf{jj}=Z_ast_err_out{jj}(mask2);
        bgf{jj}=bg_out{jj}(mask2);
        prf{jj}=pr_out{jj}(mask2);
        simidf{jj}=simid_out{jj}(mask2);
    end

    %% output to Vutara
    vutarax=cat(1,xoutf{:});
    vutaray=cat(1,youtf{:});
    vutaraz=cat(1,zoutf{:});
    vutarat=cat(1,toutf{:});
    vutaraI=cat(1,Ioutf{:});
    vutarall=cat(1,lloutf{:});
    vutarabg=cat(1,bgf{:});
    vutarazcon=cat(1,zconf{:});
    vutaracrlb=cat(1,crlbf{:});
    vutarazaster=cat(1,zasterf{:});
    vutarazerr=cat(1,zerrf{:});
    vutarapr=cat(1,prf{:});

    vutarasim=cat(1,simidf{1});
    smInfofinal=smInfonew(vutarasim,:);
    vutaraconx=smInfofinal(:,7);
    vutaracony=smInfofinal(:,8);
    
    % SIMFLUX
    [flag]=W4PiSMS2vutarav3(currentfolder,[savename '_ll'],1,{vutarax},{vutaray},{vutaraz},{floor(vutarat/100)},{vutaraI},{vutaracrlb},{vutarall},{vutarabg},{vutarazcon},{vutarazerr},{vutaraconx},{vutaracony});

    % save result
    para=tmpld.para;
    save([currentfolder savename '_' channel 'v20_60'],'vutarax','vutaray','vutaraz','vutarat','vutarall','vutaraI','vutarabg','vutarazcon','vutaracrlb','vutarazaster','vutarazerr','vutarapr','smInfofinal','vutaraconx','vutaracony','para');
    programStatus_string =[programStatus_string;['Save Vutarax File to:' currentfolder savename '_ll\' 'particles.csv']];
    set(handles.programStatus,'string',programStatus_string)
    drawnow update

    %% output image
    coords=[];
    pixel_SR=10;
    coords(:,1)=vutarax/pixel_SR;
    coords(:,2)=vutaray/pixel_SR;
    sz=imagesz*pixelsz/pixel_SR;
    im=cHistRecon(sz,sz,single(coords(:,2)),single(coords(:,1)),0);
    gaussim=gaussf(im,[1 1]);
    str2=([currentfolder savename '_gauss_1.tif']);
    writeim(gaussim,str2,'tiff',1);
    
    % color coded image
    bin=100;
    segnum=ceil((max(vutaraz)-min(vutaraz))/bin);
    [rch,gch,bch]=srhist_color(imagesz,pixelsz/pixel_SR,vutaray/pixelsz,vutarax/pixelsz,min(vutaraz)-vutaraz,segnum);
    rchsm=gaussf(rch,1);
    gchsm=gaussf(gch,1);
    bchsm=gaussf(bch,1);
    max_chsm=max([ceil(max(rchsm)),ceil(max(gchsm)),ceil(max(bchsm))]);
    max_chsm=min(max_chsm,255);
    rchsmst=imstretch_linear(rchsm,0,max_chsm,0,1000);
    gchsmst=imstretch_linear(gchsm,0,max_chsm,0,1000);
    bchsmst=imstretch_linear(bchsm,0,max_chsm,0,1000);
    colorim=joinchannels('RGB',rchsmst,gchsmst,bchsmst);
    str3=([currentfolder savename '_gauss_colorim_1.tif']);
    writeim(colorim,str3,'TIFF');
    
    disp('Reconstruction done');
end
