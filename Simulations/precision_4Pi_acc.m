close all
clear
clc

addpath("HelperFun\")

%% hyper parameters for PSF model used for fit
paraSim.NA = 1.35;                                                  % numerical aperture of obj
paraSim.refmed = 1.406;                                             % refractive index of sample medium
paraSim.refcov = 1.518;                                             % refractive index of converslip
paraSim.refimm = 1.406;                                             % refractive index of immersion oil
paraSim.lambda = 600;                                               % wavelength of emission
paraSim.objStage0_upper = -0;                                       % nm, initial objStage0 position,relative to focus at coverslip
paraSim.objStage0_lower = -0;                                       % nm, initial objStage0 position,relative to focus at coverslip
paraSim.zemit0_upper = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_upper);  % reference emitter z position, nm, distance of molecule to coverslip
paraSim.zemit0_lower = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_lower);  % reference emitter z position, nm, distance of molecule to coverslip

paraSim. pixelSizeX = 120;                                          % nm, pixel size of the image
paraSim. pixelSizeY = 120;                                          % nm, pixel size of the image
paraSim.Npupil = 64;                                                % sampling at the pupil plane

paraSim.aberrations(:,:,1) = [2,-2,0.0; 2,2,-0.1; 3,-1,0.0; 3,1,0.0; 4,0,0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];
paraSim.aberrations(:,:,2) = [2,-2,0.0; 2,2,0.1; 3,-1,0.0; 3,1,0.0; 4,0,0.0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];

paraSim.aberrations(:,3,:) =  paraSim.aberrations(:,3,:)*paraSim.lambda;

paraSim.offset = [0 0];
paraSim.phaseshift = [0 ,pi/2, pi, 3*pi/2];

%% parameters for molecules for simulation
Nmol = 101;
Npixels = 31;
paraSim.Nmol = Nmol;
paraSim.sizeX = Npixels;
paraSim.sizeY = Npixels;

Nfits = 500;
iter_times = 100;
Npixels = 15;
Nphotons_total = 8000;
bg_total = 12 * 2;

Nphotons_channel = Nphotons_total / 4;
Nphotons = Nphotons_channel * [1 1 1 1];
bg = bg_total / 4 *[1 1 1 1];

zrange = linspace(-500,500,Nmol)*1;

paraSim.xemit = (-200+400*rand(1,Nmol))*0;                              % nm
paraSim.yemit = (-200+400*rand(1,Nmol))*0;                              % nm
paraSim.zemit = linspace(-600,600,Nmol)*1;                              % nm
% paraSim.zemit = zrange;
paraSim.objStage = linspace(-1000,1000,Nmol)*0;                         % nm

[PSFs,PSFsUpper,PSFsLower,WaberrationUpper, WaberrationLower] = vectorPSF_4Pi(paraSim);

%% generate IAB model
ipalm_im  = PSFs;

phaseshift = paraSim.phaseshift;
k = 2 * pi / paraSim.lambda;    %lambda in nm
zcand = paraSim.zemit;          % if move obj stage use paraSim.objStage
dz = zcand(2) - zcand(1);

imsz = paraSim.sizeX;

%make I, A, B and their splines
I = squeeze((ipalm_im(:, :, :, 1) + ipalm_im(:, :, :, 3)) / 2);

kz2 = 2 * k * zcand';
kz2 = permute(repmat(kz2, 1, imsz, imsz), [2, 3, 1]);

F_phi1 = squeeze(ipalm_im(:, :, :, 1)) - I;
F_phi2 = squeeze(ipalm_im(:, :, :, 2)) - I;
phi1 = phaseshift(1);
phi2 = phaseshift(2);

A = (F_phi1 .* sin(kz2 + phi2) - F_phi2 .* sin(kz2 + phi1)) / sin(phi2 - phi1);
B = (-F_phi1 .* cos(kz2 + phi2) + F_phi2 .* cos(kz2 + phi1)) / sin(phi2 - phi1);

%check: restore PSF in the 4th quadrant
phi4 = phaseshift(4);
check_PSF = I + A .* cos(kz2 + phi4) + B .* sin(kz2 + phi4);
Ispline = Spline3D_interp(I);
Aspline = Spline3D_interp(A);
Bspline = Spline3D_interp(B);

%% test fitter
lambdanm = paraSim.lambda;

z0 = round(length(zcand) / 2) - 1;
% k = 2 * pi * paraSim.refimm / lambdanm; %lambda in nm

PSF.Ispline = Ispline;
PSF.Aspline = Aspline;
PSF.Bspline = Bspline;

%% Generate PSF
phi0 = [0, pi/2, pi, 1.5 * pi] + pi / 4;

ground_truth.x(:,1) = Npixels/2 + 0 *rand([Nmol 1]);
ground_truth.y(:,1) = Npixels/2 + 0 *rand([Nmol 1]);
ground_truth.N = repmat(Nphotons, [Nmol 1]);
ground_truth.bg =repmat(bg, [Nmol 1]);
ground_truth.znm = zrange';
ground_truth.zspline = ground_truth.znm / dz + z0;
ground_truth.phi =  wrapTo2Pi(2 * k * ground_truth.znm);

scale = 0;
ground_truth.x(:,2)=ground_truth.x(:,1)+(rand(Nmol, 1))*scale;
ground_truth.y(:,2)=ground_truth.y(:,1)+(rand(Nmol, 1))*scale;

ground_truth.x(:,3)=ground_truth.x(:,1)+(rand(Nmol, 1))*scale;
ground_truth.y(:,3)=ground_truth.y(:,1)+(rand(Nmol, 1))*scale;

ground_truth.x(:,4)=ground_truth.x(:,1)+(rand(Nmol, 1))*scale;
ground_truth.y(:,4)=ground_truth.y(:,1)+(rand(Nmol, 1))*scale;

% coordinates for simulation
coordinates = zeros(Nmol,4,length(phi0));
for kk=1:1:length(phi0)
    coordinates(:,:,kk) = [ground_truth.x(:,kk) ground_truth.y(:,kk) ground_truth.zspline ground_truth.phi];
end

PSFs = simSplinePSF_4pi(Npixels, PSF, ground_truth.N, ground_truth.bg, coordinates, phi0);  %simulate images


RMSEX = zeros(length(zrange),1);
RMSEY = zeros(length(zrange),1);
RMSEZ = zeros(length(zrange),1);


RMSEX_save = zeros(length(zrange),1);
RMSEY_save = zeros(length(zrange),1);
RMSEZ_save = zeros(length(zrange),1);

for iter = 1:iter_times
    iter
    n = 1;
    for i = 1:len(zrange)
        z = zrange(i);
        PSF = PSFs(:,:,i,:);
        imstack = repmat(PSF, [1, 1, Nfits, 1]);
        imstack = poissrnd(imstack, Npixels, Npixels, Nfits, 4);

        [P,CRLB,LL] = mleFit_LM(single(sum(imstack, 4)),4,50,2,0,1,0);
        X=P(:,1);
        Y=P(:,2);
        [P,CRLB, LogL] = GPUmleFit_LM_4Pi_old(imstack,0.9);

        offset = 0;
        phi =P(:,6+offset);
        phi = wrapTo2Pi(phi);

        z_phi = z_from_phi_YL(repmat(z,Nfits,1)./dz+z0, phi, k, z0, dz);
        RMSEX(n,1) = std(X);
        RMSEY(n,1) = std(Y);
        RMSEZ(n,1) = std(z_phi);
        n=n+1;
    end

    RMSEY_save = RMSEY + RMSEY_save;
    RMSEX_save = RMSEX + RMSEX_save;
    RMSEZ_save = RMSEZ + RMSEZ_save;
end

RMSEX = RMSEX_save/iter_times;
RMSEY = RMSEY_save/iter_times;
RMSEZ = RMSEZ_save/iter_times;

plot(zrange, RMSEX*paraSim.pixelSizeX)
hold on
plot(zrange, RMSEY*paraSim.pixelSizeY)
hold on
plot(zrange, RMSEZ)
axis([-600 600 1 10])

data_save.z_save = zrange;
data_save.STD_x = RMSEX*paraSim.pixelSizeX;
data_save.STD_y = RMSEY*paraSim.pixelSizeY;
data_save.STD_z = RMSEZ;
save('Fitting_std_4pi.mat','data_save');