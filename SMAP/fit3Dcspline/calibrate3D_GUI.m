%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  This file is part of GPUmleFit_LM Fitter.
%  
%  GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  GPUmleFit_LM Fitter is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with GPUmleFit_LM Fitter.  If not, see <http://www.gnu.org/licenses/>.
%  
%  
%  Additional permission under GNU GPL version 3 section 7
%  
%%
classdef calibrate3D_GUI<handle
    properties
        guihandles
        smappos
    end
    methods
        function obj=calibrate3D_GUI(varargin)  
            %constructur: make GUI   
            if ~isdeployed
            addpath('shared')
            addpath('bfmatlab')
            javaaddpath('ImageJ/plugins/bioformats_package.jar')
            end
            
            if nargin>0 %called from our propriety fitting software SMAP: extended funtionality. Hidden if called directly
                extended=true;
                figureheight=730;
            else
                extended=false;
                figureheight=630;
            end
            
            h=figure('Name','3D calibration','MenuBar','none','ToolBar','none');
            initPosition = h.Position;
            h.Position=[initPosition(1), initPosition(2)- figureheight+initPosition(4),450, figureheight];
            top=h.Position(4)-10;
            vsep=24;
            if ispc
                fontsize=12;
            else 
                fontsize=14;
            end
            xpos1=10;
            xw=100;
            hatitle='left';
            obj.guihandles.handle=h;
            obj.guihandles.title=uicontrol('style','text','String','Calibrate PSF model for MLE fit from bead stacks. (c) 2017 Ries lab','Position',[xpos1,top-vsep+10,xw*4.5,vsep],'FontSize',10,'HorizontalAlignment',hatitle,'FontWeight','bold');
            
            obj.guihandles.selectfiles=uicontrol('style','pushbutton','String','Select camera files','Position',[xpos1,top-2*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.selectfiles_callback);
            obj.guihandles.selectfiles.TooltipString='Select image files with bead stacks. You can select several files from different locations with the file select dialog box opend';
            obj.guihandles.filelist=uicontrol('style','listbox','String','','Position',[xpos1+1.5*xw,top-4*vsep,xw*2.5,vsep*3],'FontSize',fontsize);
            obj.guihandles.filelist.TooltipString='List of image files used for calibration. To change this list, use select camera files';
            obj.guihandles.selectoutputfile=uicontrol('style','pushbutton','String','Select output file','Position',[xpos1,top-5*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.selectoutputfile_callback);
            obj.guihandles.selectoutputfile.TooltipString='Select file name for output calibration file. E.g. bead_astig_3dcal.mat or bead2d_3dcal.mat';
            obj.guihandles.outputfile=uicontrol('style','edit','String','bead_3dcal.mat','Position',[xpos1+1.5*xw,top-5*vsep,xw*2.5,vsep],'FontSize',fontsize);
            obj.guihandles.outputfile.TooltipString='Name of the output file';
            
          
            ha='right';
            obj.guihandles.csplinet=uicontrol('style','text','String','General parameters: ','Position',[xpos1,top-7*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            obj.guihandles.dzt=uicontrol('style','text','String','Distance between frames (nm)','Position',[xpos1,top-8*vsep,xw*2.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.dz=uicontrol('style','edit','String','10','Position',[xpos1+2.5*xw,top-8*vsep,xw*1,vsep],'FontSize',fontsize);
            obj.guihandles.dz.TooltipString=sprintf('Distance in nm between frames. By convention, these are objective positions (not corrected for refractive index mismatch). \n A spacing between 10 nm and 50 nm works well ');
            obj.guihandles.dzt.TooltipString=obj.guihandles.dz.TooltipString;
            
            obj.guihandles.modalityt=uicontrol('style','text','String','3D modality ','Position',[xpos1,top-9*vsep,xw*2.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.modality=uicontrol('style','popupmenu','String',{'astigmatism','arbitrary','2D PSF'},'Value',1,'Position',[xpos1+2.5*xw,top-9*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.modality_callback);
            obj.guihandles.modality.TooltipString='Select the kind of PSF. Astigmatic, arbitrary (e.g. saddle-point, double-helix), or unmodified 2D';
            obj.guihandles.modalityt.TooltipString=obj.guihandles.modality.TooltipString;
            
            
            obj.guihandles.corrzt=uicontrol('style','text','String','Correct bead z-positions using ','Position',[xpos1,top-10*vsep,xw*2.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.corrzselect=uicontrol('style','popupmenu','String',{'none','cross-correlation','shape (astig)'},...
                'Value',2,'Position',[xpos1+2.5*xw,top-10*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.zcorr_callback);
            obj.guihandles.corrzselect.TooltipString=sprintf('Way of correcting for different z positions (e.g. due to a z shift between data sets or tilted coverslip). \n none: use absolute original positions. \n cross-correlation: use 3D cross-correlation on a volume defined with the parameter: frames used for CC \n shape: for astigmatism only: determine z from the frame whare sigma_x=sigma_y');
            obj.guihandles.corrzt.TooltipString=obj.guihandles.corrzselect.TooltipString;
            
            
            obj.guihandles.zcorrframest=uicontrol('style','text','String','frames to use for CC: ','Position',[xpos1+1.5*xw,top-11*vsep,xw*2,vsep],'FontSize',fontsize,'Visible','on','HorizontalAlignment',ha);
            obj.guihandles.zcorrframes=uicontrol('style','edit','String','50','Position',[xpos1+3.5*xw,top-11*vsep,xw*.5,vsep],'FontSize',fontsize,'Visible','on');
            obj.guihandles.zcorrframes.TooltipString=sprintf('Number of frames around focus used to calculate 3D cross-correlation and thus x, y and z shifts. \n Should correspond to 300-1500 nm (depends on dz). \n Too small value leads to poor z alignment.');
            obj.guihandles.zcorrframest.TooltipString=obj.guihandles.zcorrframes.TooltipString;
            
            obj.guihandles.filtert=uicontrol('style','text','String','Filter size for peak finding','Position',[xpos1,top-12*vsep,xw*2.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.filter=uicontrol('style','edit','String','2','Position',[xpos1+2.5*xw,top-12*vsep,xw*1,vsep],'FontSize',fontsize);
            obj.guihandles.filter.TooltipString=sprintf('Gaussian filter for peak finding (sigma in pixels). For split PSFs (e.g. double-helix) choose larger value to segment centroid positions of the beads, not individual lobes.');
            obj.guihandles.filtert.TooltipString=obj.guihandles.filter.TooltipString;
            
            obj.guihandles.mindistancet=uicontrol('style','text','String','Minimum distance (pixels)','Position',[xpos1,top-13*vsep,xw*2.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.mindistance=uicontrol('style','edit','String','20','Position',[xpos1+2.5*xw,top-13*vsep,xw*1,vsep],'FontSize',fontsize);
            obj.guihandles.mindistance.TooltipString=sprintf('Minimum distance between beads (in pixels). If beads are closer together, they are removed and not used for averaging. Helps eliminating background contaminations from close by beads');
            obj.guihandles.mindistancet.TooltipString=obj.guihandles.mindistance.TooltipString;           
     
            obj.guihandles.csplinet=uicontrol('style','text','String','Cspline parameters: ','Position',[xpos1,top-15*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            obj.guihandles.roisizet=uicontrol('style','text','String','ROI size: X,Y (pixels): ','Position',[xpos1,top-16*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.ROIxy=uicontrol('style','edit','String','27','Position',[xpos1+2*xw,top-16*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.roisizezt=uicontrol('style','text','String','Z (frames): ','Position',[xpos1+2.5*xw,top-16*vsep,xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.ROIz=uicontrol('style','edit','String','','Position',[xpos1+3.5*xw,top-16*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.roisizet.TooltipString=sprintf('Size of the volume for which cspline coefficients are calculated. \n Should be larger than the ROI used for fitting. \n x,y: typically 17-31 pixels, z: number of frames in stack');
            obj.guihandles.ROIxy.TooltipString=obj.guihandles.roisizet.TooltipString;
            obj.guihandles.roisizezt.TooltipString=obj.guihandles.roisizet.TooltipString;
            obj.guihandles.ROIz.TooltipString=obj.guihandles.roisizet.TooltipString;
            obj.guihandles.roisizezt.Visible='off';
            obj.guihandles.ROIz.Visible='off';
            
            
            obj.guihandles.smootht=uicontrol('style','text','String','Smoothing parameter in Z: ','Position',[xpos1,top-17*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.smoothz=uicontrol('style','edit','String','1','Position',[xpos1+2*xw,top-17*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.smoothz.TooltipString=sprintf('Smoothing paramter in z. Too large values lead to a broadened PSF and loss in accuracy, too small value leads to stripe artifacts. Typically 0.3-5');
            obj.guihandles.smootht.TooltipString=obj.guihandles.smoothz.TooltipString;
            
            obj.guihandles.gausst=uicontrol('style','text','String','Gauss fit parameters: ','Position',[xpos1,top-19*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            obj.guihandles.gaussmint=uicontrol('style','text','String','Range (nm). minimum: ','Position',[xpos1,top-20*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.gaussmin=uicontrol('style','edit','String','-700','Position',[xpos1+2*xw,top-20*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.gaussmaxt=uicontrol('style','text','String','maximum: ','Position',[xpos1+2.5*xw,top-20*vsep,xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.gaussmax=uicontrol('style','edit','String','700','Position',[xpos1+3.5*xw,top-20*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.gaussroit=uicontrol('style','text','String','ROI size (pixels): ','Position',[xpos1,top-21*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.gaussroi=uicontrol('style','edit','String','19','Position',[xpos1+2*xw,top-21*vsep,xw*.5,vsep],'FontSize',fontsize);
           
            
            obj.guihandles.gaussmint.TooltipString=sprintf('z-range (nm) around focus for which to caluclate Gaussian z-calibration. \n Use this to exclude parts with systematically bad fits. This also depends on the ROI size.');
            obj.guihandles.gaussmin.TooltipString=obj.guihandles.gaussmint.TooltipString;
            obj.guihandles.gaussmax.TooltipString=obj.guihandles.gaussmint.TooltipString;
            
            obj.guihandles.gaussroit.TooltipString=sprintf('Size of ROI (in pixels) used for Gaussian fitting. Use the same as in later experiments. \n Larger size gives higher accuracy (maximum 21 pixels allowed by fitter). \n For dense data use smaller ROI sizes to avoid contamination with close-by fluorophores.');
            obj.guihandles.gaussroi.TooltipString=obj.guihandles.gaussroit.TooltipString;
         
            obj.guihandles.run=uicontrol('style','pushbutton','String','Calculate bead calibration','Position',[xpos1,2.5*vsep,xw*4,vsep],'FontSize',fontsize,'Callback',@obj.run_callback);
            %obj.guihandles.help=uicontrol('style','pushbutton','String','Help','Position',[xpos1+xw,top-23*vsep,xw*2,vsep],'FontSize',fontsize,'Callback',@obj.help_callback);
                      
            obj.guihandles.status=uicontrol('style','text','String','Status','Position',[xpos1,1*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment','left');
            
            if extended  %called from our propriety fitting software SMAP: extended funtionality. Hidden if called directly          
                obj.guihandles.posfromsmap=uicontrol('style','checkbox','String','SMAP positions','Position',[xpos1,top-3*vsep,xw*1.5,vsep],'FontSize',fontsize,'Value',false);
                obj.guihandles.posfromsmap.TooltipString='Use positions determined by SMAP. This allows for filtering and manual deletion of beads, as well as the use of ROIs.';
                
                obj.guihandles.spatialcalt=uicontrol('style','text','String','Spatially resolved calibration: ','Position',[xpos1,top-23*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
                obj.guihandles.spatialmode=uicontrol('style','popupmenu','String',{'none','horizontal split','vertical split','M x N tiles','coordinates'},'Position',[xpos1,top-24*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha,'Callback',{@spatialselect_callback,obj});
                obj.guihandles.spatial_xtext=uicontrol('style','text','String','xtext','Position',[xpos1,top-25*vsep,xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha,'Visible','off');
                obj.guihandles.spatial_xval=uicontrol('style','edit','String','','Position',[xpos1+xw,top-25*vsep,xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha,'Visible','off');
                obj.guihandles.spatial_ytext=uicontrol('style','text','String','ytext','Position',[xpos1+2*xw,top-25*vsep,xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha,'Visible','off');
                obj.guihandles.spatial_yval=uicontrol('style','edit','String','','Position',[xpos1+3*xw,top-25*vsep,xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha,'Visible','off');   
                obj.guihandles.emgain=uicontrol('style','checkbox','String','EM gain used (mirrored)','Position',[xpos1,top-26*vsep,2*xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha); 
                obj.smappos=varargin{1};
            else
                set(h, 'HandleVisibility', 'off'); %not affected by close all command
            end
            modality_callback(obj,0,0)
        end
        function selectfiles_callback(obj,a,b)
            sf=selectManyFiles;
            sf.guihandles.filelist.String=(obj.guihandles.filelist.String);
            waitfor(sf.handle);
            obj.guihandles.filelist.String=sf.filelist;
            obj.guihandles.filelist.Value=1;
            if isempty(obj.guihandles.outputfile.String)|| strcmp(obj.guihandles.outputfile.String,'bead_3dcal.mat')
                [path,file]=fileparts(sf.filelist{1});
                obj.guihandles.outputfile.String=[path filesep file '_3Dcorr.mat'];
                
            end
        end
        function selectoutputfile_callback(obj,a,b)
            of=obj.guihandles.outputfile.String;
            if isempty(of)
                of='_corr3D.mat';
            end
            [f,p]=uiputfile(of);
            if f
            obj.guihandles.outputfile.String=[p,f];
            end
        end
        function modality_callback(obj,a,b)
            astigonly={'gausst','gaussmin','gaussmint','gaussmax','gaussmaxt','gaussroi','gaussroit'};
            switch obj.guihandles.modality.String{obj.guihandles.modality.Value}
                case 'astigmatism'
                    vis='on';
                    obj.guihandles.corrzselect.String={'none','cross-correlation','shape (astig)'};
                    
                otherwise
                    vis='off';
                    obj.guihandles.corrzselect.String={'none','cross-correlation'};
                    obj.guihandles.corrzselect.Value=min(obj.guihandles.corrzselect.Value,2);
            end
            for k=1:length(astigonly)
                obj.guihandles.(astigonly{k}).Visible=vis;
            end
            zcorr_callback(obj,0,0)
        end
        function zcorr_callback(obj,a,b)
            corrz={'zcorrframest','zcorrframes'};
            if contains(obj.guihandles.corrzselect.String{obj.guihandles.corrzselect.Value},'cross')
                vis='on';
            else
                vis='off';
            end
            for k=1:length(corrz)
                obj.guihandles.(corrz{k}).Visible=vis;
            end
        end
        function out=run_callback(obj,a,b)
            p.filelist=obj.guihandles.filelist.String;
            p.outputfile=obj.guihandles.outputfile.String;
            p.dz=str2double(obj.guihandles.dz.String);
            p.modality=obj.guihandles.modality.String{obj.guihandles.modality.Value};
            p.zcorr=obj.guihandles.corrzselect.String{obj.guihandles.corrzselect.Value};
            p.ROIxy=str2double(obj.guihandles.ROIxy.String);
            p.ROIz=str2double(obj.guihandles.ROIz.String);
%             p.smoothxy=str2double(obj.guihandles.smoothxy.String);
            p.smoothxy=0;
            p.smoothz=str2double(obj.guihandles.smoothz.String);
            p.gaussrange=[str2double(obj.guihandles.gaussmin.String) str2double(obj.guihandles.gaussmax.String)];
            p.filtersize=str2double(obj.guihandles.filter.String);
            p.zcorrframes=str2double(obj.guihandles.zcorrframes.String);
            p.gaussroi=str2double(obj.guihandles.gaussroi.String);
            p.status=obj.guihandles.status;
            p.mindistance=str2double(obj.guihandles.mindistance.String);
            if isempty(p.filelist)
                warndlg('please select image files first')
                return
            end
            
            if ~isempty(obj.smappos) %called from SMAP
                p.smap=true;
                p.imageRoi=obj.smappos.imageROI;
                if obj.guihandles.posfromsmap.Value %use positions passed on from SMAP
                    p.beadpos=obj.smappos.positions;
                end
                if isfield(obj.smappos,'framerangeuse')
                    p.framerangeuse=obj.smappos.framerangeuse;
                end
                %determine xrange, yrange for spatial calibration
                switch obj.guihandles.spatialmode.Value
                    case 1 %all
                        if isfield(obj.smappos,'xrange')
                            p.xrange=obj.smappos.xrange;
                        else
                            p.xrange=[-inf inf];
                        end
                        if isfield(obj.smappos,'yrange')
                            p.yrange=obj.smappos.yrange;
                        else
                            p.yrange=[-inf inf];
                        end                       
                        
                    case 2%split horz
                        midp=str2double(obj.guihandles.spatial_xval.String);
                        p.yrange=[-inf  inf];p.xrange=[-inf midp inf];
                    case 3% split vert
                       midp=str2double(obj.guihandles.spatial_xval.String);
                        p.yrange=[-inf  midp inf];p.xrange=[-inf inf];
                    case 4%MxN
                        mn=str2num(obj.guihandles.spatial_xval.String);
                        imsize=str2num(obj.guihandles.spatial_yval.String);
                        if length(imsize)==1
                            imsize(2)=imsize(1);
                        end
                        if length(mn)==1
                            mn(2)=mn(1);
                        end                       
                        p.xrange=0:floor(imsize(1)/mn(1)):imsize(1);p.xrange(end)=imsize(1);
                        p.yrange=0:floor(imsize(2)/mn(2)):imsize(1);p.yrange(end)=imsize(2);
                    case 5
                        p.xrange=str2num(obj.guihandles.spatial_xval.String);
                        p.yrange=str2num(obj.guihandles.spatial_yval.String);
                end
                p.emgain=obj.guihandles.emgain.Value;
            else
                p.emgain=false;
            end
            calibrate3D(p);
        end    
        function help_callback(obj,a,b)
            helpstring={'write documentation'};
            f=figure;
            h=uicontrol('Parent',f,'Style','text','Units','normalized','Position',[0 0 1 1],'HorizontalAlignment','Left');
            h.String=textwrap(h,helpstring);
            
        end
    end
end


function spatialselect_callback(a,b,obj) %for SMAP extension and spatial calibration option: change visibility of controls
xt='';xv='';yt='';yv='';
xtv='off';xvv='off';ytv='off';yvv='off';
switch obj.guihandles.spatialmode.Value
    case 1 %none
    case {2,3} % split
        xt='split at pixel'; xv='256'; xtv='on';xvv='on';
    case 4 % tiles
        xt='M, N'; xv='3, 3'; xtv='on';xvv='on';
        yt='Image size x,y'; yv='512, 512'; ytv='on';yvv='on';
    case 5 %coordinates
        xt='x coords'; xv='0:256:512'; xtv='on';xvv='on';
        yt='y coords'; yv='0:256:512'; ytv='on';yvv='on';  
end
obj.guihandles.spatial_xtext.String=xt;obj.guihandles.spatial_xtext.Visible=xtv;
obj.guihandles.spatial_ytext.String=yt;obj.guihandles.spatial_ytext.Visible=ytv;
obj.guihandles.spatial_xval.String=xv;obj.guihandles.spatial_xval.Visible=xvv;
obj.guihandles.spatial_yval.String=yv;obj.guihandles.spatial_yval.Visible=yvv;
end
