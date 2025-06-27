close all
clear
clc

addpath("HelperFun\")

%% hyper parameters for PSF model used for fit
paraSim.NA = 1.35;                                                      % numerical aperture of obj
paraSim.refmed = 1.406;                                                 % refractive index of sample medium
paraSim.refcov = 1.518;                                                 % refractive index of converslip
paraSim.refimm = 1.406;                                                 % refractive index of immersion oil
paraSim.lambda = 600;                                                   % wavelength of emission
paraSim.objStage0 = -0;                                                 % nm, initial objStage0 position,relative to focus at coverslip
paraSim.zemit0 = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0);  % reference emitter z position, nm, distance of molecule to coverslip
paraSim. pixelSizeX = 120;                                              % nm, pixel size of the image
paraSim. pixelSizeY = 120;                                              % nm, pixel size of the image
paraSim.Npupil = 64;                                                    % sampling at the pupil plane

paraSim.aberrations = [2,-2,0.0; 2,2,0.1; 3,-1,0.0; 3,1,0.0; 4,0,0.00; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,-0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];
paraSim.aberrations(:,3) =  paraSim.aberrations(:,3)*paraSim.lambda;

%% Simulate parameters
Nmol = 101;
Npixels = 31;
Nphotons = 4000;
bg = 12;
Nfits = 500;
Iters = 100;

paraSim.Nmol = Nmol;
paraSim.sizeX = Npixels;
paraSim.sizeY = Npixels;

paraSim.xemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.yemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.zemit = linspace(-900,900,Nmol)*1;                             %nm
zlist = linspace(-500,500,Nmol)*1;
paraSim.objStage = linspace(-1000,1000,Nmol)*0;                        %nm

%% generate PSF and Calibration
[PSFs,Waberration] = vectorPSF_Final(paraSim);
PSF_coffs = Spline3D_interp(PSFs);

z_len = Nmol;
zlist_calib = linspace(-800,800,Nmol);
Result = zeros(3, z_len);
dz = paraSim.zemit(2) - paraSim.zemit(1);
z0 = round(z_len / 2) - 1;
zspline = zlist / dz + z0;
zspline_calib = zlist_calib / dz + z0;

Npixels = 21;
gt.x(:,1) = Npixels/2 + 0*rand([z_len 1]);
gt.y(:,1) = Npixels/2 + 0*rand([z_len 1]);
gt.N = Nphotons + 0*rand([z_len 1]);
gt.bg = bg+ 0*rand([z_len 1]);
coordinates = horzcat(gt.x, gt.y, zspline_calib');
imstack = simSplinePSF(Npixels, PSF_coffs, gt.N, gt.bg, coordinates);   %simulate images

[P,CRLB,LL] = mleFit_LM(single(imstack),4,50,2,0,1,0);
wxlist = P(:,5);
wylist = P(:,6);

% with wx and wy
zrange_calib = [min(zlist_calib), max(zlist_calib)];                    %range of z
zcali_wxy_N = 5;
zcali_pz2wx = polyfit (zlist_calib, wxlist, zcali_wxy_N);
zcali_pz2wy = polyfit (zlist_calib, wylist, zcali_wxy_N);

%% Generate all the PSF
Npixels = 15;
ground_truth.x(:,1) = Npixels/2 + 0*rand([Nmol 1]);
ground_truth.y(:,1) = Npixels/2 + 0*rand([Nmol 1]);
ground_truth.N = Nphotons+ 0*rand([Nmol 1]);
ground_truth.bg = bg+ 0*rand([Nmol 1]);

ground_truth.zspline = zspline';
coordinates = horzcat(ground_truth.x, ground_truth.y, ground_truth.zspline);
PSFs = simSplinePSF(Npixels, PSF_coffs, ground_truth.N, ground_truth.bg, coordinates);%simulate images

for iter = 1:Iters
    iter
    for i=1:length(zlist)
        PSF = PSFs(:,:,i);
        imstack = repmat(PSF, [1,1,Nfits]);
        imstack = poissrnd(imstack, Npixels, Npixels, Nfits);

        [P,CRLB,LL] = mleFit_LM(single(imstack),4,50,2,0,1,0);
        x = P(:,1);
        y = P(:,2);
        stdx = std(x)*paraSim.pixelSizeX;
        stdy = std(y)*paraSim.pixelSizeY;
        wxlist = P(:,5);
        wylist = P(:,6);
        retz = fitZwithWXY(double(wxlist), double(wylist), double(zcali_pz2wx), double(zcali_pz2wy), zrange_calib);
        stdz = std(retz);
        Result(1, i) = Result(1, i) + stdx;
        Result(2, i) = Result(2, i) + stdy;
        Result(3, i) = Result(3, i) + stdz;
    end
end

Result = Result / Iters;
plot(zlist, Result(1,:));
hold on
plot(zlist, Result(2,:));
hold on
plot(zlist, Result(3,:));
data_save.z_save = zlist;
data_save.STD_x = Result(1,:);
data_save.STD_y = Result(2,:);
data_save.STD_z = Result(3,:);
save('Fitting_std_single_len.mat','data_save');