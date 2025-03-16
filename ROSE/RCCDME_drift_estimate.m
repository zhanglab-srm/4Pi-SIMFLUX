function [xcorr,ycorr,drift,rc]=RCCDME_drift_estimate(smInfo_x,smInfo_y,smInfo,framebin)


%% Prepare data
localizations(:,1)=smInfo_x;
localizations(:,2)=smInfo_y;

N=length(smInfo_x(:));

crlb(1:N,1:2)=0.02;
framenum=smInfo(:,3);

% if handles.parameter.drift_z
%     localizations(:,3)=smInfo(:,3);
% end

%% User inputs RCC drift computation
zoom = 5;
sigma = 1;
maxpairs = 1000;
timebins = ceil((max(framenum)/framebin));
usecuda = true;

%% RCC computation
[drift_xy] = rcc(localizations, framenum, timebins, zoom, sigma, maxpairs, usecuda);

%% User inputs for drift estimation - !keep data types alive!
coarse_est = false;                     % Coarse drift estimation (bool)
precision_est = false;                  % Precision estimation (bool)
coarse_frames_per_bin = int32(100);     % Number of bins for coarse est. (int32)
framesperbin = int32(10);               % Number of frames per bin (int32)
maxneighbors_coarse = int32(1000);      % Max neighbors for coarse and precision est. (int32)
maxneighbors_regular = int32(1000);     % Max neighbors for regular est. (int32)
coarseSigma= single(mean(crlb));        % Localization precision for coarse estimation (single/float)
max_iter_coarse = int32(1000);          % Max iterations coarse est. (int32)
max_iter = int32(10000);                % Max iterations (int32)
gradientstep = single(1e-6);            % Gradient (single/float)

%% Drift estimation
drift = dme_estimate(localizations, framenum, crlb, drift_xy, usecuda, coarse_frames_per_bin, ...
    framesperbin, maxneighbors_coarse,maxneighbors_regular, coarseSigma, max_iter_coarse, max_iter, gradientstep, precision_est);
%% plot
rc=figure();
for i = 1:size(drift,2)
    subplot(size(drift,2),1,i);
    plot(1:1:max(framenum), drift(:,i), 'LineWidth', 2);
    hold on
    plot(1:1:max(framenum), drift_xy(:,i), 'LineWidth', 2);
    legend('Estimated drift (DME)', 'Estimated drift (RCC)') ;
    %     legend('Estimated drift (RCC)') ;
end

%% Output
drift_per_frame = drift(smInfo(:,3),:);
xcorr = double(localizations(:,1)-drift_per_frame(:,1));
ycorr = double(localizations(:,2)-drift_per_frame(:,2));