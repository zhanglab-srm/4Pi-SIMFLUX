function W4PiSMS_getPosition(handles)

mainFolder = get(handles.pathMainfolder,'string');
index_selected = get(handles.pathSubfolder,'value');
folderlist = get(handles.pathSubfolder,'string');
folder_selected = folderlist(index_selected);

%parameters
det_thresh=str2double(get(handles.det_thresh,'string'));
center1 = str2double(get(handles.center1,'string'));
center2 = str2double(get(handles.center2,'string'));
center3 = str2double(get(handles.center3,'string'));
center4 = str2double(get(handles.center4,'string'));

centers = [center1 center2 center3 center4];
para.det_thresh = det_thresh;
para.centers = centers;
para.pixelsz = handles.pixelsz;   % pixel size, nm

for ii=1:numel(folder_selected)
    close all;
    
    currentfolder=folder_selected{ii};
    I=find(currentfolder=='\',1,'last');
    childfoldername=currentfolder(I+1:end);
    
    I=find(currentfolder=='\',1,'last');
    parentFolder = currentfolder(1:I);

    files=dir([currentfolder '\*.dcimg']);
    imageName = files(1).name;
    if contains(imageName,'_642_')
        channel = '642';
    elseif contains(imageName,'_561_')
        channel = '561';
    else
        x = inputdlg({'Channel'},'Please specify the channel',1,{'488'});
        channel = x{1};
    end
    para.channel = channel;
    handles.channel = channel;
   
    %% 4PiSMS sanalysis
    foldername=[currentfolder '\'];
    nameroot=('*.dcimg');
    savename=childfoldername;    
    
    if ii==1
        tmpf=dir([mainFolder '*_' channel '_*Ast*.mat']);
        if numel(tmpf)>1
            set(handles.programStatus,'string','More than one *Ast*.mat calibration file detected!');
            return;
        elseif numel(tmpf)==0
            set(handles.programStatus,'string','No *Ast*.mat calibration file detected!');
            return;
        end    
        astfile=tmpf.name;
        astfile=[mainFolder astfile];

        tmpf=dir([mainFolder '*_' channel '_*dphi*.mat']);    
        if numel(tmpf)>1
            set(handles.programStatus,'string','More than one *dphi*.mat calibration file detected!');
            return;
        elseif numel(tmpf)==0
            set(handles.programStatus,'string','No *dphi*.mat calibration file detected!');
            return;
        end    
        anglefile=tmpf.name;
        anglefile=[mainFolder anglefile];

        tmpf=dir([mainFolder '*_' channel '_*FMTtransform*.mat']);    
        if numel(tmpf)>1
            set(handles.programStatus,'string','More than one *FMTtransform*.mat calibration file detected!');
            return;
        elseif numel(tmpf)==0
            set(handles.programStatus,'string','No *FMTtransform*.mat calibration file detected!');
            return;
        end    
        fmfile=tmpf.name;
        fmfile=[mainFolder fmfile];
        tmp=load(fmfile);
        centers=tmp.para.centers;
        para.centers=centers;
    end
    
    datestring = datestr(now,'yyyymmdd');
    set(handles.programStatus,'string',{'Current Analized Folder:',currentfolder, ['Current Channel:' channel],'Calibration Files Used:',astfile,anglefile,fmfile,'Output:',[mainFolder savename '_' channel '_tmpresult_' datestring '.mat']})
    drawnow update
    
    [xresult yresult tresult zfresult zangresult llresult CRLBresult Iresult ...
        zast_err_result stacktot num_images imagesz zangctrresult bgresult subimstot sxtot sytot smInfotot...
        ]=W4PiSMS_RM_analysisv8stack(foldername,nameroot,fmfile,astfile,anglefile,0,handles.scmos_cali_file,det_thresh,0,0,centers,handles);
    
    para.num_images=num_images;
    para.imagesz=imagesz;
    save([parentFolder savename '_' channel '_tmpresult_' datestring],'foldername','sxtot','sytot','tresult','anglefile','astfile','xresult','yresult','zfresult',...
        'zangctrresult','bgresult','zangresult','llresult','CRLBresult','Iresult','zast_err_result','stacktot','smInfotot','para','-v7.3');
    
    %% output image
    coords=[];
    pixelsz=para.pixelsz;  
    pixel_SR=10;          % nm      
    coords(:,1)=xresult*pixelsz/pixel_SR;
    coords(:,2)=yresult*pixelsz/pixel_SR;
    sz=para.imagesz*pixelsz/pixel_SR;
    im=cHistRecon(sz,sz,single(coords(:,2)),single(coords(:,1)),0);
    gaussim=gaussf(im,[1 1]);
    str2=([parentFolder savename '_gauss_1.tif']);
    writeim(gaussim,str2,'tiff',1);   
end