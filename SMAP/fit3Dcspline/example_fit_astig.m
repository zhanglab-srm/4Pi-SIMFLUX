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

%% example script for fitting of 3D data
%Installation instructions
%  Requires Matlab 2016a or newer including the following tool boxes:
%  image analysis, curve fitting, optimization, statistics
%  GPU support for CUDA 8


%% add path to helper functions
addpath('shared')

%% make bead calibration
%run 3D calibration GUI and make 3D calibration
%e.g. using bead files: in example_data/beadstacks_3D_astig 
% save e.g. as example_data/bead_astig_3dcal.mat

% or generate calibration files programmatically:
if ~exist([pwd filesep 'example_data' filesep 'bead_astig_3dcal.mat'],'file') %only run if no calibration files have been generated, as this takes time.
    example_calibration
end


%% load bead calibration
cal=load([pwd filesep 'example_data' filesep 'bead_astig_3dcal.mat']); %load bead calibration


%% either simulate data or load experimental tiff file
mode = 0;
if mode ==1 % simulate data
    p.offset=0;
    p.conversion=1;
    
    %numlocs: number of simulated molecules. For maximum fitting speeds on  the GPU, 
    %it should be >10^5. For smaller GPUs, adust to avoid memory overflow
    %error to 10^5-10^5. For CPU fitter, use 10^3-10^4.
    numlocs=50000;
    RoiPixelsize=13; %ROI in pixels for simulation
    dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
    z0=cal.cspline.z0; % distance and midpoint of stack in spline PSF, needed to translate into nm coordinates
    dx=floor(RoiPixelsize/2); %distance of center from edge
    ground_truth.z=linspace(-600,600,numlocs)'; %define some coordinates. Alternatively, use rand
    ground_truth.x=linspace(-0.5,0.5,numlocs)';
    ground_truth.y=sin(ground_truth.x*4*pi);
    coordinates=horzcat(ground_truth.x+dx,ground_truth.y+dx,ground_truth.z/dz+z0);
    Intensity=2000; %number of photons / localization
    background=10; %number of background photons per pixel
    imstack = simSplinePSF(RoiPixelsize,cal.cspline.coeff,Intensity,background,coordinates); %simulate images
    
else %experimental data
    
    file=[pwd filesep 'example_data' filesep 'single_bead.tif']; %experimental astgmatic bead stack.
    imstackadu=readfile_tif(file); % Stack of ROIs in photons.
    p.offset=400;p.conversion=.1;
    imstack=(single(imstackadu)-p.offset)/p.conversion;% if necessary, convert ADU into photons.
    ground_truth.z=((1:size(imstack,3))'-size(imstack,3)/2+1)*10; %dz=10 nm; convert frame to nm   
end


%% fit image stacks

% load sCMOS varmap or set sCMOSvarmap to zero
% sCMOSvarmap=ones(size(imstack),'single'); %the variance map, same size as imstack, single format. 
sCMOSvarmap=0; %if scalar : use EMCCD fitter;

numlocs=size(imstack,3); 
imstack=single(imstack); %imstack needs to be in photons. The fitters require the stacks in single-format;

if isfield(cal,'gauss_zfit')&&~isempty(cal.gauss_zfit) %only if Gaussian PSF model was calibrated in the calibrate3D_GUI, we can test the Gaussian fit
    %fit Gaussian model, direct z fit, emCCD mode
    tic
    [P,CRLB]=mleFit_LM(imstack,3,50,single(cal.gauss_zfit),sCMOSvarmap,1);
    tgz=toc;
    disp(['Gauss z: ' num2str(numlocs/tgz) ' fits/s']);
    gaussz.x=P(:,1);gaussz.y=P(:,2); %x,y in pixels 
    gaussz.z=P(:,5)*1000;


    %fit elliptical Gaussian model, extract z from sx, sy. emCCD mode
    tic
    [P,CRLB]=mleFit_LM(imstack,4,50,1,sCMOSvarmap,1);
    tgsxsy=toc;
    disp(['Gaussxy: ' num2str(numlocs/tgsxsy) ' fits/s']);

    gausssxsy.x=P(:,1);gausssxsy.y=P(:,2); %x,y in pixels 
    sx=P(:,5);sy=P(:,6);
    gausssxsy.z=sxsy2z(sx,sy,cal.gauss_sx2_sy2); %z from sx, sy

else
    gaussz.x=zeros(size(imstack,3),1);gaussz.y=zeros(size(imstack,3),1); gaussz.z=zeros(size(imstack,3),1);
    gausssxsy.x=zeros(size(imstack,3),1);gausssxsy.y=zeros(size(imstack,3),1); gausssxsy.z=zeros(size(imstack,3),1);
end

%fit experimental astigmatic PSF, cspline, emCCD mode
tic
[Pcspline,CRLB,LL]=mleFit_LM(imstack,5,50,single(cal.cspline.coeff),sCMOSvarmap,1);
tspline=toc;
disp(['cspline: ' num2str(numlocs/tspline) ' fits/s']);
dx=floor(size(imstack,1)/2);
cspline.x=Pcspline(:,1)-dx;cspline.y=Pcspline(:,2)-dx; %x,y in pixels 
cspline.z=(Pcspline(:,5)-cal.cspline.z0)*cal.cspline.dz;

z_err=sqrt(CRLB(:,5))*cal.cspline.dz; %CRLB is the variance of the parameters, to obtain a localization precision we have to take the square root

% calculate error for all fits. 
zrange=600; %Only take into accoutn central part. Focus +/- zrange
inz=abs(ground_truth.z)<zrange;

%difference between fitted z and ground truth
cspline.zrel=cspline.z-ground_truth.z;
gaussz.zrel=gaussz.z-ground_truth.z;
gaussz.zrel(isnan(gaussz.zrel))=0;
gausssxsy.zrel=gausssxsy.z-ground_truth.z;
gausssxsy.zrel(isnan(gausssxsy.zrel))=0;

%numpoints
numpoints=min(2000,size(imstack,3));
range=1:round(length(ground_truth.z)/numpoints):length(ground_truth.z);
%plot fitted z vs ground-truth z
    figure(101)
    hold off
    plot(ground_truth.z(range),gaussz.zrel(range)-mean(gaussz.zrel(inz)),'r.')
    hold on
    plot(ground_truth.z(range),gausssxsy.zrel(range)-mean(gausssxsy.zrel(inz)),'.','Color',[0.7,0.7,0])
    plot(ground_truth.z(range),cspline.zrel(range)-mean(cspline.zrel(inz)),'b.')
    gs=fit(ground_truth.z(range),double(gaussz.zrel(range)-mean(gaussz.zrel(inz))),'smoothingspline','SmoothingParam',.00001);
    gss=fit(ground_truth.z(range),double(gausssxsy.zrel(range)-mean(gausssxsy.zrel(inz))),'smoothingspline','SmoothingParam',.00001);
    gz=fit(ground_truth.z(range),double(cspline.zrel(range)-mean(cspline.zrel(inz))),'smoothingspline','SmoothingParam',.00001);

    plot(ground_truth.z(range),ground_truth.z(range)*0,'k','LineWidth',3)
    plot(ground_truth.z(range),gs(ground_truth.z(range)),'r','LineWidth',3)
    plot(ground_truth.z(range),gss(ground_truth.z(range)),'Color',[0.7,0.7,0],'LineWidth',3)
    plot(ground_truth.z(range),gz(ground_truth.z(range)),'b','LineWidth',3)

    plot([-zrange -zrange],[-1000 1000],'c')
    plot([zrange zrange],[-1000 1000],'c')

    ss=std(cspline.zrel(inz))*10;
    ylim([-ss ss])
    xlabel('ground truth z (nm)')
    ylabel('z (fitted) - z (ground truth) (nm)')

    title('Accuracy')

    cspline.dz=nanstd(ground_truth.z(inz)-cspline.z(inz));
    zgauss.dz=nanstd(ground_truth.z(inz)-gaussz.z(inz));
    gausssxsy.dz=nanstd(ground_truth.z(inz)-gausssxsy.z(inz));

    legendtxt{4}='smoothing spline';
    legendtxt{3}=['spline fit. error: ' num2str(cspline.dz, '%3.1f') ' nm'];
    legendtxt{1}=['Gaussian z fit. error: ' num2str(zgauss.dz, '%3.1f') ' nm'];
    legendtxt{2}=['Gaussian sx, sy fit. error: ' num2str(gausssxsy.dz, '%3.1f') ' nm'];
    legend(legendtxt)

%calculate lateral error
if mode == 1
    cspline.dx=nanstd(ground_truth.x(inz)-cspline.x(inz)); %in pixels
    cspline.dy=nanstd(ground_truth.y(inz)-cspline.y(inz));
    dx_gaussz=nanstd(ground_truth.x(inz)-gaussz.x(inz));
    dy_gaussz=nanstd(ground_truth.y(inz)-gaussz.y(inz));
    gausssxsy.dx=nanstd(ground_truth.x(inz)-gausssxsy.x(inz));
    gausssxsy.dy=nanstd(ground_truth.y(inz)-gausssxsy.y(inz));
    %plot 3D scatter plot
   
    figure(102);
    hold off
    scatter3(ground_truth.x(range),ground_truth.y(range),ground_truth.z(range),3);
    hold on
    scatter3(cspline.x(range),cspline.y(range),cspline.z(range),5);
    xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
    legend('ground truth','cspline fit')
    title('fitted vs ground truth positions')
else
    cspline.dx=nanstd(cspline.x(inz)); %in pixels
    cspline.dy=nanstd(cspline.y(inz));
    dx_gaussz=nanstd(gaussz.x(inz));
    dy_gaussz=nanstd(gaussz.y(inz));
    gausssxsy.dx=nanstd(gausssxsy.x(inz));
    gausssxsy.dy=nanstd(gausssxsy.y(inz));
end

