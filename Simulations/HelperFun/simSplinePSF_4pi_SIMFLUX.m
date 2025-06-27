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
%  If you modify this Program, or any covered work, by
%  linking or combining it with libraries required for interaction
%  with analysis programs such as Igor Pro or Matlab,
%  the licensors of this Program grant you additional permission
%  to convey the resulting work.
%%
function [out_4pi, out_SIMFLUX] = simSplinePSF_4pi_SIMFLUX(Npixels,PSF,I,bg,cor, phi0, paraSim)
t=tic;
if (nargin <5)
   error('Minimal usage: simSplinePSF(Npixels,coeff,I,bg,cor)');
end

% if (bg == 0)
%     bg = 10^-10;
% end    

Nfits = size(cor,1);
spline_xsize = size(PSF.Ispline,1);
spline_ysize = size(PSF.Ispline,2);
spline_zsize = size(PSF.Ispline,3);
off = floor(((spline_xsize+1)-Npixels)/2);
data = zeros(Npixels,Npixels,Nfits,'single');


for kk = 1:Nfits
        xcenter = cor(kk,1,:);
        ycenter = cor(kk,2,:);
        zcenter = cor(kk,3,:);

        xc = -1*(xcenter - Npixels/2+0.5);
        yc = -1*(ycenter - Npixels/2+0.5);
        zc = zcenter - floor(zcenter);

        xstart = floor(xc);
        xc = xc - xstart;

        ystart = floor(yc);
        yc = yc - ystart;


        zstart = floor(zcenter);

for ss=1:length(phi0)
        [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc(ss)),single(yc(ss)),single(zc(ss)));

        for ii = 0:Npixels-1
            for jj = 0:Npixels-1
                 phi = cor(kk, 4,1);
%                  for p = 1:1:length(phi0)
                     temp = fAt3Dj_4Pi(ii+xstart(ss)+off,jj+ystart(ss)+off,zstart(ss),spline_xsize,spline_ysize,spline_zsize,delta_f,PSF, phi, phi0(ss));
                     model = temp * I(kk,ss);% + bg(kk, ss);
                     data(ii+1, jj+1, kk, ss) = model;
                     if temp < 0
                         disp('Neg');
                     end
%                  end
                   
            end
        end
        if toc(t)>1
            disp(kk/Nfits)
            t=tic;
        end
end
end

data_new = zeros(Npixels,Npixels,Nfits,24,'single');
cor_center_fov = Npixels/2;
for i = 1:4
    for j = 1:6
        PSF = data(:,:,:,i);
        xcenter_nm = (cor(:,1,1)-cor_center_fov)*paraSim.pixelSizeX;
        ycenter_nm = (cor(:,2,1) - cor_center_fov)*paraSim.pixelSizeY;
        if j<4  %横向
            % phase_factor = paraSim.period_k * xcenter_nm + (j-1) * paraSim.phase_shift + paraSim.phase_ini;
            phase_factor = paraSim.period_k * xcenter_nm + (j-2) * paraSim.phase_shift + paraSim.phase_ini;
            
        else % 纵向
            % phase_factor = paraSim.period_k * ycenter_nm + (j-1) * paraSim.phase_shift + paraSim.phase_ini;
            phase_factor = paraSim.period_k * ycenter_nm + (j-2) * paraSim.phase_shift + paraSim.phase_ini;
        end
        % phase_factor
        intensity = 1/6*(1+paraSim.modulation_depth*sin(phase_factor));
        for k = 1:Nfits
            data_new(:,:,k,(i-1)*6+j) = PSF(:,:,k) * intensity(k) + bg(k);
        end
    end
end
% scale = [1 1 1 1];
% for p = 1:1:length(phi0)
%     out_4pi(:, :, :, p) = (poissrnd(data(:, :, :, p)*scale(p),Npixels,Npixels,Nfits)); 
%    %out = data;
% end
out_4pi = zeros(Npixels,Npixels,Nfits,4,'single');
out_SIMFLUX = zeros(Npixels,Npixels,Nfits,6,'single');

for i=1:4
    for j = 1:6
        out_4pi(:,:,:,i) = out_4pi(:,:,:,i) + data_new(:,:,:,(i-1)*6+j);
        out_SIMFLUX(:,:,:,j) = out_SIMFLUX(:,:,:,j) + data_new(:,:,:,(i-1)*6+j);
    end
end
% out_4pi = poissrnd(out_4pi,Npixels,Npixels,Nfits,4);
% out_SIMFLUX = poissrnd(out_SIMFLUX,Npixels,Npixels,Nfits,6);
end