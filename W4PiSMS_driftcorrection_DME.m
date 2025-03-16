function [xout2,yout2,zout2,shifts]=W4PiSMS_driftcorrection_DME(xout,yout,zout,tout,num_images,pixelsz,crlbout,rccstr)

addpath(genpath('./dme-v1.2.1'))
dmestr=insertBefore(rccstr, ".mat", "_DME");
rccstr=insertBefore(rccstr, ".mat", "_RCC");
framenum = tout' + 1;

if exist(dmestr, 'file')
    load(dmestr);
    dmeDrift = shifts';
    if exist(rccstr, 'file')
        load(rccstr);
        rccDrift = drift_xyz';
    end
else
    %% config
    usecuda = true;
    timebins = ceil(max(framenum)/num_images);
    N = length(xout);
    if N>10000000 || timebins>1000
        id = xout>=64 & xout<=192 & yout>=64 & yout<=192 & crlbout(:,1)<0.06 & crlbout(:,2)<0.06;
        xout1=xout(id);
        yout1=yout(id);
        zout1=zout(id);
        tout1=tout(id);
        crlbout=crlbout(id,:);
        framenum1=framenum(id);
        
        imagesz=200;
        im=cHistRecon(imagesz*3,imagesz*3,single(yout1)*3,single(xout1)*3,0);
        gaussim=gaussf(im,[1 1]);
        dipshow(gaussim);
        pause(eps)
    else
        xout1=xout;
        yout1=yout;
        zout1=zout;
        tout1=tout;
        framenum1=framenum;
    end
    coords = single([xout1, yout1, zout1/1000]);
    crlbout(:,3) = 0.005;
    crlb = single(crlbout);

    %% User inputs RCC drift computation
    zoom = 6;
    sigma = 1;
    maxpairs = 10000;

    %% User inputs for drift estimation - !keep data types alive!
    coarse_est = true;                      % Coarse drift estimation (bool)
    precision_est = false;                  % Precision estimation (bool)
    coarse_frames_per_bin = int32(100);     % Number of bins for coarse est. (int32)
    framesperbin = int32(10);               % Number of frames per bin (int32)
    maxneighbors_coarse = int32(1000);      % Max neighbors for coarse and precision est. (int32)
    maxneighbors_regular = int32(1000);     % Max neighbors for regular est. (int32)
    coarseSigma= single([0.04,0.04,0.005]);    % Localization precision for coarse estimation (single/float)
    max_iter_coarse = int32(20000);         % Max iterations coarse est. (int32)
    max_iter = int32(20000);                % Max iterations (int32)
    gradientstep = single(1e-6);            % Gradient (single/float)

    %% RCC computation
    %RCC 3D
    if exist(rccstr, 'file')
        load(rccstr);
        rccDrift = drift_xyz';
    else
        tic;
%         [drift_xyz] = rcc3D(coords, framenum, timebins, zoom, sigma, maxpairs,  usecuda);
        [drift_xyz] = W4PiSMS_driftcorrection_RedunLSv13(xout1*pixelsz,yout1*pixelsz,zout1,tout1,num_images);
        drift_xyz(:,1:2) = drift_xyz(:,1:2)/pixelsz;    % pixel
        drift_xyz(:,3) = drift_xyz(:,3)/1000;           % um
        drift_xyz(:,1) = drift_xyz(:,1)-mean(drift_xyz(:,1));
        drift_xyz(:,2) = drift_xyz(:,2)-mean(drift_xyz(:,2));
        drift_xyz(:,3) = drift_xyz(:,3)-mean(drift_xyz(:,3));
        toc;
        rccDrift = drift_xyz';
        save(rccstr, 'drift_xyz');
    end

    %% Drift estimation
    tic;
    [drift, iter1, iter2] = dme_estimate(coords, framenum1, crlb, drift_xyz, usecuda, coarse_frames_per_bin, ...
        framesperbin, maxneighbors_coarse, maxneighbors_regular, coarseSigma, max_iter_coarse,max_iter, gradientstep, precision_est);
    elapsed = toc;
    disp('Time spent ' + string(toc));
    dmeDrift = drift';

    %% result
    shifts = drift;
    save(dmestr, 'shifts');
end

xout2 = xout.*pixelsz - shifts(framenum,1).*pixelsz;
yout2 = yout.*pixelsz - shifts(framenum,2).*pixelsz;
zout2 = zout - shifts(framenum,3).*1000;
if exist("iter1","var")
    dmestr = insertBefore(dmestr, ".mat", ['_' num2str(iter1) 'it_' num2str(iter2) 'it_' num2str(int32(elapsed)) 's']);
end

%% plot
figure('Position',[50 50 1920 1080]);
subplot(2,1,1);
if exist('rccDrift', 'var')
    plot(rccDrift','LineWidth',0.1);
    title('RCC drift');
    legend('x(px)','y(px)','z(um)');
end

subplot(2,1,2);
plot(dmeDrift','LineWidth',0.1);
title('DME Drift');
legend('x(px)','y(px)','z(um)');
saveas(gcf, strrep(dmestr, '.mat', '.png'));
